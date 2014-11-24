#include "Exception.h"
#include "FormulationHelper.h"
#include "FormulationUpdateOO2.h"

using namespace std;

FormulationUpdateOO2::FormulationUpdateOO2(DDMContextOO2& context){
  // Save DDMContext //
  this->context = &context;

  // Get Domain and FunctionSpace from DDMContext //
  fspaceG = &context.getFunctionSpaceG(); // DDM field: testing and unknwon
  fspace  = &context.getFunctionSpace();  // e field
  ddomain = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateOO2 needs a uniform mesh");

  // a & b //
  this->a = context.getA();
  this->b = context.getB();

  // Get Basis //
  basis = &fspaceG->getBasis(eType);

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
    Complex(-1, 0)     * lGin->getTerm(equationI, 0, elementId) +
    Complex(+2, 0) * a *   lU->getTerm(equationI, 0, elementId) +
    Complex(-2, 0) * b *  lDU->getTerm(equationI, 0, elementId);
}

const FunctionSpace& FormulationUpdateOO2::field(void) const{
  return *fspaceG;
}

const FunctionSpace& FormulationUpdateOO2::test(void) const{
  return *fspaceG;
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
  context->getSystem().getSolution(sol, 0);

  // Pre-evalution
  const fullMatrix<double>& gCFF = gaussFF->getPoints();
  const fullMatrix<double>& gCGG = gaussGG->getPoints();

  basis->preEvaluateFunctions(gCFF);
  basis->preEvaluateDerivatives(gCGG);

  // New RHS
  lGin =
    new TermProjectionField<Complex>(*jacFF, *basis, *gaussFF, *fspaceG, ddm);

  lU   =
    new TermProjectionField<Complex>(*jacFF, *basis, *gaussFF, *fspace, sol);

  lDU  =
    new TermProjectionGrad<Complex>(*jacGG, *basis, *gaussGG, *fspace, sol);
}
