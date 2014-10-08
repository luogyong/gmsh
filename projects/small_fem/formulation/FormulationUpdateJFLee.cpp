#include "Exception.h"
#include "FormulationJFLee.h"
#include "FormulationHelper.h"
#include "FormulationUpdateJFLee.h"

using namespace std;

FormulationUpdateJFLee::FormulationUpdateJFLee(DDMContextJFLee& context){
  // Save DDMContext //
  this->context = &context;

  // Get Domain and FunctionSpace (primary and auxiliary) from DDMContext //
  ddomain = &context.getDomain();
  ffield  = &context.getFunctionSpace();
  fPhi    = &context.getPhi();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateJFLee needs a uniform mesh");

  // Basis //
  basis = &ffield->getBasis(eType);

  // Gaussian Quadrature //
  gauss = new Quadrature(eType, basis->getOrder(), 2);

  // Pre-evalution //
  const fullMatrix<double>& gC = gauss->getPoints();
  basis->preEvaluateFunctions(gC);

  // Jacobian //
  jac = new GroupOfJacobian(*ddomain, gC, "invert");

  // Init Auxiliary Solution //
  FormulationHelper::initDofMap(*fPhi, *ddomain, phi);

  // Local Terms //
  lGout = new TermGradGrad<double>(*jac, *basis, *gauss);

  // NB: lGin and lPhi will be computed at update, when vol Sys is available
  lGin  = NULL;
  lPhi  = NULL;
}

FormulationUpdateJFLee::~FormulationUpdateJFLee(void){
  delete lGout;
  if(lGin)
    delete lGin;
  if(lPhi)
    delete lPhi;

  delete jac;
  delete gauss;
}

Complex FormulationUpdateJFLee::weak(size_t dofI, size_t dofJ,
                                     size_t elementId) const{
  return
    Complex(lGout->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationUpdateJFLee::rhs(size_t equationI, size_t elementId) const{
  return
    Complex(-1, 0) * lGin->getTerm(equationI, 0, elementId) +
    Complex(+2, 0) * lPhi->getTerm(equationI, 0, elementId);
}

const FunctionSpace& FormulationUpdateJFLee::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationUpdateJFLee::test(void) const{
  return *ffield;
}

const GroupOfElement& FormulationUpdateJFLee::domain(void) const{
  return *ddomain;
}

bool FormulationUpdateJFLee::isBlock(void) const{
  return true;
}

void FormulationUpdateJFLee::update(void){
  // Delete RHS (lGin & lPhi)
  if(lGin)
    delete lGin;
  if(lPhi)
    delete lPhi;

  // Get DDM Dofs and auxiliary solutions (at border) from DDMContext //
  const map<Dof, Complex>& ddm = context->getDDMDofs(); // ddm
  context->getSystem().getSolution(phi, 0);             // phi (aux)

  // Pre-evalution
  const fullMatrix<double>& gC = gauss->getPoints();
  basis->preEvaluateFunctions(gC);

  // New RHS
  lGin = new TermProjectionGrad<Complex>(*jac, *basis, *gauss, *ffield, ddm);
  lPhi = new TermProjectionGrad<Complex>(*jac, *basis, *gauss, *fPhi  , phi);
}
