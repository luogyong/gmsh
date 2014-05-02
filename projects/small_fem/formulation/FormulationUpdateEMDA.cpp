#include "Exception.h"
#include "FormulationHelper.h"
#include "FormulationUpdateEMDA.h"

using namespace std;

FormulationUpdateEMDA::FormulationUpdateEMDA(DDMContext& context){
  // Check if EMDA DDMContext //
  if(context.typeDDM != DDMContext::typeEMDA)
    throw Exception("FormulationUpdateEMDA needs a EMDA DDMContext");

  // Save DDMContext //
  this->context = &context;

  // Get Domain and FunctionSpace from DDMContext //
  fspace  = &context.getFunctionSpace();
  ddomain = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateEMDA needs a uniform mesh");

  // Wavenumber & Chi //
  this->k   = context.k;
  this->chi = context.EMDA_Chi;

  // Basis //
  basis = &fspace->getBasis(eType);

  // Gaussian Quadrature //
  gauss = new Quadrature(eType, basis->getOrder(), 2);

  // Pre-evalution //
  const fullMatrix<double>& gC = gauss->getPoints();
  basis->preEvaluateFunctions(gC);

  // Jacobian //
  jac = new GroupOfJacobian(*ddomain, gC, "jacobian");

  // Init Volume Solution //
  FormulationHelper::initDofMap(*fspace, *ddomain, sol);

  // Local Terms //
  lGout = new TermFieldField<double>(*jac, *basis, *gauss);

  // NB: lGin & lU will be computed at update, when volume System is available
  lGin = NULL;
  lU   = NULL;
}

FormulationUpdateEMDA::~FormulationUpdateEMDA(void){
  delete lGout;
  if(lGin)
    delete lGin;
  if(lU)
    delete lU;

  delete jac;
  delete gauss;
}

Complex FormulationUpdateEMDA::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return Complex(lGout->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationUpdateEMDA::rhs(size_t equationI, size_t elementId) const{
  return
    Complex(-1      , 0     ) * lGin->getTerm(0, equationI, elementId) +
    Complex(+2 * chi, -2 * k) *   lU->getTerm(0, equationI, elementId);
}

const FunctionSpace& FormulationUpdateEMDA::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationUpdateEMDA::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationUpdateEMDA::domain(void) const{
  return *ddomain;
}

bool FormulationUpdateEMDA::isBlock(void) const{
  return true;
}

void FormulationUpdateEMDA::update(void){
  // Delete RHS (lGin & lU)
  if(lGin)
    delete lGin;
  if(lU)
    delete lU;

  // Get DDM Dofs & Volume solution (at border) from DDMContext //
  const map<Dof, Complex>& ddm = context->getDDMDofs();
  context->system->getSolution(sol, 0);

  // Pre-evalution
  const fullMatrix<double>& gC = gauss->getPoints();
  basis->preEvaluateFunctions(gC);

  // New RHS
  lGin  = new TermProjectionField<Complex>(*jac, *basis, *gauss, *fspace, ddm);
  lU    = new TermProjectionField<Complex>(*jac, *basis, *gauss, *fspace, sol);
}
