#include <cmath>

#include "Exception.h"

#include "FormulationOSRCScalarOne.h"
#include "FormulationOSRCScalarTwo.h"
#include "FormulationOSRCScalarThree.h"
#include "FormulationOSRCScalarFour.h"
#include "FormulationOSRCScalar.h"

#include "FormulationOSRCHelper.h"

using namespace std;

FormulationOSRCScalar::FormulationOSRCScalar(DDMContextOSRCScalar& context){
  // Save DDMContext //
  this->context = &context;

  // Get Domain and auxiliary FunctionSpaces from DDMContext //
  const GroupOfElement&                     dom = context.getDomain();
  const vector<const FunctionSpaceScalar*>& aux =
    context.getAuxFunctionSpace();

  // Save field FunctionSpace
  field    = &context.getFunctionSpace();  // Unknown Field
  ffspaceG = &context.getFunctionSpaceG(); // DDM Field

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = dom.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOSRCScalar needs a uniform mesh");

  // Get Basis //
              basisU = &field->getBasis(eType); // Saved for update()
  const size_t order = basisU->getOrder();

  // Get Auxiliary Basis //
  const Basis& basisP = aux[0]->getBasis(eType);

  // k, keps and NPade //
  double k     = context.getWavenumber();
  Complex keps = context.getComplexWavenumber();

  // Pade //
  int    NPade = context.getNPade();
  double theta = context.getRotation();
  vector<Complex> A(NPade);
  vector<Complex> B(NPade);

  Complex C0 = FormulationOSRCHelper::padeC0(NPade, theta);

  for(int j = 0; j < NPade; j++)
    A[j] = FormulationOSRCHelper::padeA(j + 1, NPade, theta);

  for(int j = 0; j < NPade; j++)
    B[j] = FormulationOSRCHelper::padeB(j + 1, NPade, theta);

  // Gaussian Quadrature //
  gaussFF = new Quadrature(eType, order, 2); // Saved for update()
  Quadrature gaussGG(eType, order - 1, 2);

  const fullMatrix<double>& gCFF = gaussFF->getPoints();
  const fullMatrix<double>& gCGG = gaussGG.getPoints();

  // Pre-evaluate //
  basisU->preEvaluateFunctions(gCFF);
  basisU->preEvaluateDerivatives(gCGG);

  basisP.preEvaluateFunctions(gCFF);
  basisP.preEvaluateDerivatives(gCGG);

  jacFF = new GroupOfJacobian(dom, gCFF, "jacobian"); // Saved for update()
  GroupOfJacobian jacGG(dom, gCGG, "invert");

  // Get DDM Dofs from DDMContext //
  const map<Dof, Complex>& ddm = context.getDDMDofs();

  // Local Terms //
  UU   = new TermFieldField<double>(*jacFF, *basisU, *basisU, *gaussFF);
  dPdU = new TermGradGrad<double>  ( jacGG,  basisP, *basisU,  gaussGG);
  PP   = new TermFieldField<double>(*jacFF,  basisP,  basisP, *gaussFF);
  dPdP = new TermGradGrad<double>  ( jacGG,  basisP,  basisP,  gaussGG);
  UP   = new TermFieldField<double>(*jacFF, *basisU,  basisP, *gaussFF);
  RHS  =
    new TermProjectionField<Complex>(*jacFF, *basisU, *gaussFF, *ffspaceG, ddm);

  // Formulations //
  // NB: FormulationOSRCScalar is a friend
  //     of FormulationOSRCScalar{One,Two,Three,Four,} !
  //     So it can instanciate those classes...

  // Save FormulationOSRCScalarOne for update()
  formulationOne = new FormulationOSRCScalarOne(dom,*field,k,C0,*UU,*RHS);

  // Then push it in list
  fList.push_back(formulationOne);

  // Loop on phi[j]
  FormulationBlock<Complex>* f[3];
  for(int j = 0; j < NPade; j++){
    f[0]=new FormulationOSRCScalarTwo  (dom,*aux[j],*field,k,keps,A[j],  *dPdU);
    f[1]=new FormulationOSRCScalarThree(dom,*aux[j],       keps,B[j],*PP,*dPdP);
    f[2]=new FormulationOSRCScalarFour (dom,*field ,*aux[j],               *UP);

    fList.push_back(f[0]);
    fList.push_back(f[1]);
    fList.push_back(f[2]);
  }
}

FormulationOSRCScalar::~FormulationOSRCScalar(void){
  // Iterate & Delete Formulations //
  list<const FormulationBlock<Complex>*>::iterator end = fList.end();
  list<const FormulationBlock<Complex>*>::iterator it  = fList.begin();

  for(; it !=end; it++)
    delete *it;

  // Delete terms //
  delete   UU;
  delete dPdU;
  delete   PP;
  delete dPdP;
  delete   UP;
  delete  RHS;

  // Delete update stuffs //
  delete jacFF;
  delete gaussFF;
}

const list<const FormulationBlock<Complex>*>&
FormulationOSRCScalar::getFormulationBlocks(void) const{
  return fList;
}

bool FormulationOSRCScalar::isBlock(void) const{
  return false;
}

void FormulationOSRCScalar::update(void){
  // Delete RHS
  delete RHS;

  // Get DDM Dofs from DDMContext
  const map<Dof, Complex>& ddm = context->getDDMDofs();

  // Pre-evalution
  const fullMatrix<double>& gCFF = gaussFF->getPoints();
  basisU->preEvaluateFunctions(gCFF);

  // New RHS
  RHS = new TermProjectionField<Complex>(*jacFF,*basisU,*gaussFF,*ffspaceG,ddm);

  // Update FormulationOSRCScalarOne (formulationOne):
  //                                             this FormulationBlock holds RHS
  formulationOne->update(*RHS);
}
