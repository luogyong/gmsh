#include "Exception.h"
#include "FormulationHelper.h"
#include "FormulationUpdateOO2.h"

using namespace std;

FormulationUpdateOO2::FormulationUpdateOO2(DDMContext& context){
  // Check if OO2 DDMContext //
  if(context.typeDDM != DDMContext::typeOO2)
    throw Exception("FormulationUpdateOO2 needs a OO2 DDMContext");

  // Save DDMContext //
  this->context = &context;

  // Get Domain and FunctionSpace from DDMContext //
  fspace  = &context.getFunctionSpace();
  ddomain = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateOO2 needs a uniform mesh");

  // a & b //
  this->a = context.OO2_A;
  this->b = context.OO2_B;

  // Get Basis //
  basis = &fspace->getBasis(eType);

  // Gaussian Quadrature (Field - Field & Grad - Grad) //
  gaussFF = new Quadrature(eType, basis->getOrder()    , 2);
  gaussGG = new Quadrature(eType, basis->getOrder() - 1, 2);

  // Pre-evalution //
  const fullMatrix<double>& gCFF = gaussFF->getPoints();
  const fullMatrix<double>& gCGG = gaussGG->getPoints();

  basis->preEvaluateFunctions(gCFF);
  basis->preEvaluateDerivatives(gCGG);

  // Jacobians //
  jacFF = new GroupOfJacobian(*ddomain, gCFF, "jacobian");
  jacGG = new GroupOfJacobian(*ddomain, gCGG, "invert");

  // Init Volume Solution //
  FormulationHelper::initDofMap(*fspace, *ddomain, sol);

  // Local Terms //
  lGout = new TermFieldField<double>(*jacFF, *basis, *gaussFF);

  // NB: lGin, lU & lDU will be computed at update, when vol System is available
  lGin = NULL;
  lU   = NULL;
  lDU  = NULL;
}

FormulationUpdateOO2::~FormulationUpdateOO2(void){
  delete lGout;
  if(lGin)
    delete lGin;
  if(lU)
    delete lU;
  if(lDU)
    delete lDU;

  delete gaussFF;
  delete gaussGG;
  delete jacFF;
  delete jacGG;
}

Complex FormulationUpdateOO2::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return Complex(lGout->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationUpdateOO2::rhs(size_t equationI, size_t elementId) const{
  return
    Complex(-1, 0)     * lGin->getTerm(0, equationI, elementId) +
    Complex(+2, 0) * a *   lU->getTerm(0, equationI, elementId) +
    Complex(-2, 0) * b *  lDU->getTerm(0, equationI, elementId);
}

const FunctionSpace& FormulationUpdateOO2::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationUpdateOO2::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationUpdateOO2::domain(void) const{
  return *ddomain;
}

bool FormulationUpdateOO2::isBlock(void) const{
  return true;
}

void FormulationUpdateOO2::update(void){
  // Delete RHS (lGin, lU & lDU)
  if(lGin)
    delete lGin;
  if(lU)
    delete lU;
  if(lDU)
    delete lDU;

  // Get DDM Dofs & Volume solution (at border) from DDMContext //
  const map<Dof, Complex>& ddm = context->getDDMDofs();
  context->system->getSolution(sol, 0);

  // Pre-evalution
  const fullMatrix<double>& gCFF = gaussFF->getPoints();
  const fullMatrix<double>& gCGG = gaussGG->getPoints();

  basis->preEvaluateFunctions(gCFF);
  basis->preEvaluateDerivatives(gCGG);

  // New RHS
  lGin =
    new TermProjectionField<Complex>(*jacFF, *basis, *gaussFF, *fspace, ddm);

  lU   =
    new TermProjectionField<Complex>(*jacFF, *basis, *gaussFF, *fspace, sol);

  lDU  =
    new TermProjectionGrad<Complex>(*jacGG, *basis, *gaussGG, *fspace, sol);
}
