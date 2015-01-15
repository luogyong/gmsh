#include <cmath>
#include "Exception.h"

#include "FormulationOSRCVectorOne.h"
#include "FormulationOSRCVectorTwo.h"
#include "FormulationOSRCVectorThree.h"
#include "FormulationOSRCVectorFour.h"
#include "FormulationOSRCVectorFive.h"
#include "FormulationOSRCVectorSix.h"
#include "FormulationOSRCVectorSeven.h"
#include "FormulationOSRCVectorEight.h"
#include "FormulationOSRCVectorNine.h"

#include "FormulationOSRCHelper.h"
#include "FormulationOSRCVector.h"

using namespace std;

FormulationOSRCVector::FormulationOSRCVector(DDMContextOSRCVector& context){
  // Save DDMContext //
  this->context = &context;

  // Get Domain and Auxiliary FunctionSpaces from DDMContext //
  const GroupOfElement&                     dom = context.getDomain();
  const vector<const FunctionSpaceVector*>& phi = context.getPhiFunctionSpace();
  const vector<const FunctionSpaceScalar*>& rho = context.getRhoFunctionSpace();
  const FunctionSpaceVector&                  R = context.getRFunctionSpace();

  // Save Field FunctionSpace
  field   = &context.getFunctionSpace();  // Unknown Field
  fspaceG = &context.getFunctionSpaceG(); // DDM Field

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = dom.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOSRCVector needs a uniform mesh");

  // Get Basis for field //
  const Basis& basisE = field->getBasis(eType);
  const size_t order  = basisE.getOrder();

  // Get auxiliary bases //
               basisR   = &(    R.getBasis(eType)); // Saved for update
  const Basis& basisPhi = phi[0]->getBasis(eType);
  const Basis& basisRho = rho[0]->getBasis(eType);

  // k, kEpsilon and NPade //
  double  k     = context.getWavenumber();
  Complex kE    = context.getComplexWavenumber();
  int     NPade = context.getNPade();
  double  theta = context.getRotation();

  Complex         R0;
  vector<Complex> A(NPade);
  vector<Complex> B(NPade);

  R0 = FormulationOSRCHelper::padeR0(NPade, theta);

  for(int i = 0; i < NPade; i++)
    A[i] = FormulationOSRCHelper::padeA(i + 1, NPade, theta);

  for(int i = 0; i < NPade; i++)
    B[i] = FormulationOSRCHelper::padeB(i + 1, NPade, theta);

  // Get DDM Dofs from DDMContext //
  const map<Dof, Complex>& ddm = context.getDDMDofs();

  // Gaussian Quadrature //
  gauss = new Quadrature(eType, order, 2); // Saved for update()
  const fullMatrix<double>& gC = gauss->getPoints();

  // Pre-evaluate //
  basisE.preEvaluateFunctions(gC);
  basisE.preEvaluateDerivatives(gC);

  basisR->preEvaluateFunctions(gC);
  basisR->preEvaluateDerivatives(gC);

  basisPhi.preEvaluateFunctions(gC);
  basisPhi.preEvaluateDerivatives(gC);

  basisRho.preEvaluateFunctions(gC);
  basisRho.preEvaluateDerivatives(gC);

  jac = new GroupOfJacobian(dom, gC, "both"); // Saved for update()

  // Local Terms //
  RHS  = new TermProjectionGrad<Complex>(*jac, *basisR  , *gauss, *fspaceG,ddm);
  RE   = new TermGradGrad<double>       (*jac, *basisR  ,  basisE  , *gauss);

  ER   = new TermGradGrad<double>       (*jac,  basisE  , *basisR  , *gauss);
  cEcR = new TermCurlCurl<double>       (*jac,  basisE  , *basisR  , *gauss);
  RR   = new TermGradGrad<double>       (*jac, *basisR  , *basisR  , *gauss);
  PR   = new TermGradGrad<double>       (*jac,  basisPhi, *basisR  , *gauss);

  RP   = new TermGradGrad<double>       (*jac, *basisR  ,  basisPhi, *gauss);
  PP   = new TermGradGrad<double>       (*jac,  basisPhi,  basisPhi, *gauss);
  cPcP = new TermCurlCurl<double>       (*jac,  basisPhi,  basisPhi, *gauss);
  dRoP = new TermGradGrad<double>       (*jac,  basisRho,  basisPhi, *gauss);

  RoRo = new TermFieldField<double>     (*jac,  basisRho,            *gauss);
  PdRo = new TermGradGrad<double>       (*jac,  basisPhi,  basisRho, *gauss);

  // Formulations //
  // NB: FormulationOSRCVector is a friend
  //     of FormulationOSRCVector{One,Two,Three,Four,Five,
  //                              Six,Seven,Eight,Nine} !
  //     So it can instanciate those classes...

  FormulationBlock<Complex>* f[6];

  f[0]= new FormulationOSRCVectorOne  (dom,  R    , *field, k     , *RE);
  f[1]= new FormulationOSRCVectorTwo  (dom, *field,  R    , kE, R0, *ER, *cEcR);
  f[2]= new FormulationOSRCVectorThree(dom,  R    ,  R    ,     R0, *RR, *RHS);

  // Save FormulationOSRCVectorThree for update()
  formulationThree = static_cast<FormulationOSRCVectorThree*>(f[2]);

  // Then push them in list
  fList.push_back(f[0]);
  fList.push_back(f[1]);
  fList.push_back(f[2]);

  // Loop on Pade terms
  for(int i = 0; i < NPade; i++){
    f[0]= new FormulationOSRCVectorFour (dom,*phi[i], R     ,R0,A[i],B[i],*PR);

    f[1]= new FormulationOSRCVectorFive (dom, R     ,*phi[i],            *RP);
    f[2]= new FormulationOSRCVectorSix  (dom,*phi[i],*phi[i],kE,B[i],*PP,*cPcP);
    f[3]= new FormulationOSRCVectorSeven(dom,*rho[i],*phi[i],   B[i],    *dRoP);

    f[4]= new FormulationOSRCVectorEight(dom,*rho[i],*rho[i],kE,         *RoRo);
    f[5]= new FormulationOSRCVectorNine (dom,*phi[i],*rho[i],            *PdRo);

    fList.push_back(f[0]);
    fList.push_back(f[1]);
    fList.push_back(f[2]);
    fList.push_back(f[3]);
    fList.push_back(f[4]);
    fList.push_back(f[5]);
  }
}

FormulationOSRCVector::~FormulationOSRCVector(void){
  // Iterate & Delete Formulations //
  list<const FormulationBlock<Complex>*>::iterator end = fList.end();
  list<const FormulationBlock<Complex>*>::iterator it  = fList.begin();

  for(; it !=end; it++)
    delete *it;

  // Delete terms //
  delete RHS;
  delete RE;
  delete ER;
  delete cEcR;
  delete RR;
  delete PR;
  delete RP;
  delete PP;
  delete cPcP;
  delete dRoP;
  delete RoRo;
  delete PdRo;

  // Delete update stuffs //
  delete jac;
  delete gauss;
}

const list<const FormulationBlock<Complex>*>&
FormulationOSRCVector::getFormulationBlocks(void) const{
  return fList;
}

bool FormulationOSRCVector::isBlock(void) const{
  return false;
}

void FormulationOSRCVector::update(void){
  // Delete RHS
  delete RHS;

  // Get DDM Dofs from DDMContext
  const map<Dof, Complex>& ddm = context->getDDMDofs();

  // Pre-evalution
  const fullMatrix<double>& gC = gauss->getPoints();
  basisR->preEvaluateFunctions(gC);

  // New RHS
  RHS = new TermProjectionGrad<Complex>(*jac, *basisR, *gauss, *fspaceG, ddm);

  // Update FormulationOSRCVectorThree (formulationThree):
  //                                             this FormulationBlock holds RHS
  formulationThree->update(*RHS);
}
