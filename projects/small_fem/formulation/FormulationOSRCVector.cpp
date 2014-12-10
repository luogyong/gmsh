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
  field   = &context.getFunctionSpace();  // Testing Field and Unknown Field
  fspaceG = &context.getFunctionSpaceG(); // DDM Field

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = dom.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOSRCVector needs a uniform mesh");

  // Get Vectorial Basis //
  basisV = &(field->getBasis(eType)); // Saved for update()
  const size_t order = basisV->getOrder();

  // Get scalar Basis
  const Basis& basisS = rho[0]->getBasis(eType);

  // k, kEpsilon and NPade //
  double  k     = context.getWavenumber();
  Complex kE    = context.getComplexWavenumber();
  int     NPade = context.getNPade();
  double theta  = context.getRotation();

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

  // Local Terms //
  basisV->preEvaluateFunctions(gC);
  basisV->preEvaluateDerivatives(gC);

  basisS.preEvaluateFunctions(gC);
  basisS.preEvaluateDerivatives(gC);

  jac = new GroupOfJacobian(dom, gC, "both"); // Saved for update()

  // NB: Since the Formulations share the same basis functions,
  //     the local terms will be the same !
  //     It's the Dof numbering imposed by the function spaces that will differ
  RHS = new TermProjectionGrad<Complex>(*jac, *basisV, *gauss, *fspaceG, ddm);
  GG  = new TermGradGrad<double>       (*jac, *basisV,         *gauss);
  dFG = new TermGradGrad<double>       (*jac,  basisS,*basisV, *gauss);
  GdF = new TermGradGrad<double>       (*jac, *basisV, basisS, *gauss);
  CC  = new TermCurlCurl<double>       (*jac, *basisV,         *gauss);
  FF  = new TermFieldField<double>     (*jac,  basisS,         *gauss);

  // Formulations //
  // NB: FormulationOSRCVector is a friend
  //     of FormulationOSRCVector{One,Two,Three,Four,Five,
  //                              Six,Seven,Eight,Nine} !
  //     So it can instanciate those classes...

  FormulationBlock<Complex>* f[6];

  f[0] = new FormulationOSRCVectorOne  (dom,  R    , *field, k     , *GG);
  f[1] = new FormulationOSRCVectorTwo  (dom, *field,  R    , kE, R0, *GG, *CC);
  f[2] = new FormulationOSRCVectorThree(dom,  R    ,  R    ,     R0, *GG, *RHS);

  // Save FormulationOSRCVectorThree for update()
  formulationThree = static_cast<FormulationOSRCVectorThree*>(f[2]);

  // Then push them in list
  fList.push_back(f[0]);
  fList.push_back(f[1]);
  fList.push_back(f[2]);

  // Loop on Pade terms
  for(int i = 0; i < NPade; i++){
    f[0] = new FormulationOSRCVectorFour (dom,*phi[i], R     ,R0,A[i],B[i],*GG);

    f[1] = new FormulationOSRCVectorFive (dom, R     ,*phi[i],             *GG);
    f[2] = new FormulationOSRCVectorSix  (dom,*phi[i],*phi[i],kE,B[i], *GG,*CC);
    f[3] = new FormulationOSRCVectorSeven(dom,*rho[i],*phi[i],   B[i],    *dFG);

    f[4] = new FormulationOSRCVectorEight(dom,*rho[i],*rho[i],kE,          *FF);
    f[5] = new FormulationOSRCVectorNine (dom,*phi[i],*rho[i],            *GdF);

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
  delete GG;
  delete dFG;
  delete GdF;
  delete CC;
  delete FF;

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
  basisV->preEvaluateFunctions(gC);

  // New RHS
  RHS = new TermProjectionGrad<Complex>(*jac, *basisV, *gauss, *fspaceG, ddm);

  // Update FormulationOSRCVectorThree (formulationThree):
  //                                             this FormulationBlock holds RHS
  formulationThree->update(*RHS);
}
