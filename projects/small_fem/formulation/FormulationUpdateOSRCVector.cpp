#include "Exception.h"

#include "FormulationHelper.h"
#include "FormulationOSRCHelper.h"

#include "FormulationOSRCVector.h"
#include "FormulationUpdateOSRCVector.h"

using namespace std;

FormulationUpdateOSRCVector::
FormulationUpdateOSRCVector(DDMContextOSRCVector& context){
  // Save DDMContext //
  this->context = &context;

  // Get Domain and FunctionSpace from DDMContext //
  ffspace = &context.getFunctionSpace();
  ddomain = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateOSRCVector needs a uniform mesh");

  // Wavenumber //
  double  k       = context.getWavenumber();
  this->twoJOverK = Complex(0, 2. / k);

  // Pade //
  NPade = context.getNPade();
  A.resize(NPade);
  B.resize(NPade);

  R0 = FormulationOSRCHelper::padeR0(NPade, M_PI / 4.);

  for(int j = 0; j < NPade; j++)
    A[j] = FormulationOSRCHelper::padeA(j + 1, NPade, M_PI / 4.);

  for(int j = 0; j < NPade; j++)
    B[j] = FormulationOSRCHelper::padeB(j + 1, NPade, M_PI / 4.);

  // Basis //
  basis = &ffspace->getBasis(eType);

  // Gaussian Quadrature //
  gauss = new Quadrature(eType, basis->getOrder(), 2);

  // Pre-evalution //
  const fullMatrix<double>& gC = gauss->getPoints();
  basis->preEvaluateFunctions(gC);

  // Jacobian //
  jac = new GroupOfJacobian(*ddomain, gC, "invert");

  // Auxiliary FunctionSpace //
  fR   = &(context.getRFunctionSpace());
  fPhi = &(context.getPhiFunctionSpace());

  // Init Auxiliary Solution //
  solPhi.resize(NPade);
  FormulationHelper::initDofMap(*fR,   *ddomain, solR);
  FormulationHelper::initDofMap(*fPhi, *ddomain, solPhi);

  // Local Terms //
  lGout = new TermGradGrad<double>(*jac, *basis, *gauss);

  // NB: lGin, lR and lPhi will be computed at update, when vol Sys is available
  lGin = NULL;
  lR   = NULL;
  lPhi = NULL;
}

FormulationUpdateOSRCVector::~FormulationUpdateOSRCVector(void){
  delete lGout;

  if(lGin)
    delete lGin;
  if(lR)
    delete lR;
  if(lPhi)
    delete lPhi;

  delete jac;
  delete gauss;
}

Complex FormulationUpdateOSRCVector::weak(size_t dofI, size_t dofJ,
                                          size_t elementId) const{
  return
    Complex(lGout->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationUpdateOSRCVector::rhs(size_t equationI,
                                         size_t elementId) const{
  // Sum the rest & return
  return
    lGin->getTerm(equationI, 0, elementId)                  -
      lR->getTerm(equationI, 0, elementId) * twoJOverK * R0 +
    lPhi->getTerm(equationI, 0, elementId) * twoJOverK;
}

const FunctionSpace& FormulationUpdateOSRCVector::field(void) const{
  return *ffspace;
}

const FunctionSpace& FormulationUpdateOSRCVector::test(void) const{
  return *ffspace;
}

const GroupOfElement& FormulationUpdateOSRCVector::domain(void) const{
  return *ddomain;
}

bool FormulationUpdateOSRCVector::isBlock(void) const{
  return true;
}

void FormulationUpdateOSRCVector::update(void){
  // Delete RHS (lGin, lR, lPhi)
  if(lGin)
    delete lGin;
  if(lR)
    delete lR;
  if(lPhi)
    delete lPhi;

  // Get DDM Dofs and auxiliary solutions (at border) from DDMContext //
  const map<Dof, Complex>& ddm = context->getDDMDofs(); // ddm
  context->getSystem().getSolution(solR, 0);            // solR

  for(int i = 0; i < NPade; i++)
    context->getSystem().getSolution(solPhi[i], 0);     // solPhi[i]

  // Compute all solPhi[:] contributions and put them into solPhi[0]
  getAllPhi();

  // Get References
  map<Dof, Complex>&           alPhi = solPhi[0];
  const FunctionSpaceVector& phiZero = *(*fPhi)[0];

  // Pre-evalution
  const fullMatrix<double>& gC = gauss->getPoints();
  basis->preEvaluateFunctions(gC);

  // New RHS
  lGin = new TermProjectionGrad<Complex>(*jac, *basis, *gauss, *ffspace, ddm);
  lR   = new TermProjectionGrad<Complex>(*jac, *basis, *gauss, *fR     , solR);
  lPhi = new TermProjectionGrad<Complex>(*jac, *basis, *gauss,  phiZero, alPhi);
}

void FormulationUpdateOSRCVector::getAllPhi(void){
  // Sum all contributions into solPhi[0] //
  map<Dof, Complex>::iterator endZero = solPhi[0].end();
  map<Dof, Complex>::iterator itZero;
  map<Dof, Complex>::iterator itI;

  // Firt multiply by A[0] / B[0]
  itZero  = solPhi[0].begin();

  for(; itZero != endZero; itZero++)
    itZero->second = itZero->second * A[0] / B[0];

  // The sum other contributions
  for(int i = 1; i < NPade; i++){
    itZero  = solPhi[0].begin();
    itI     = solPhi[i].begin();

    for(; itZero != endZero; itZero++, itI++)
      itZero->second = itZero->second + itI->second * A[i] / B[i];
  }
}
