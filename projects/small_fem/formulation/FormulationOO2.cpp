#include "FormulationOO2.h"

using namespace std;

FormulationOO2::FormulationOO2(DDMContextOO2& context){
  // Save DDMContext //
  this->context = &context;

  // Get Domain and FunctionSpace from DDMContext //
  fspace  = &context.getFunctionSpace();
  ddomain = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOO2 needs a uniform mesh");

  // a & b //
  this->a = context.getA();
  this->b = context.getB();

  // Get Basis //
  basis = &fspace->getBasis(eType);

  // Gaussian Quadrature (Field - Field & Field - Projection) //
  gaussFF = new Quadrature(eType, basis->getOrder(), 2); // Saved for update()

  // Gaussian Quadrature (Grad - Grad) //
  Quadrature gaussGG(eType, basis->getOrder() - 1, 2);

  // Pre-evalution //
  const fullMatrix<double>& gCFF = gaussFF->getPoints();
  const fullMatrix<double>& gCGG = gaussGG.getPoints();

  basis->preEvaluateFunctions(gCFF);
  basis->preEvaluateDerivatives(gCGG);

  // Jacobians //
  jacFF = new GroupOfJacobian(*ddomain, gCFF, "jacobian"); // Saved for update()
  GroupOfJacobian jacGG(*ddomain, gCGG, "invert");

  // Get DDM Dofs from DDMContext //
  const map<Dof, Complex>& ddm = context.getDDMDofs();

  // Local Terms //
  localTermsFF = new TermFieldField<double>(*jacFF, *basis, *gaussFF);
  localTermsGG = new TermGradGrad<double>(jacGG, *basis, gaussGG);
  localTermsPr =
    new TermProjectionField<Complex>(*jacFF, *basis, *gaussFF, *fspace, ddm);
}

FormulationOO2::~FormulationOO2(void){
  delete localTermsFF;
  delete localTermsGG;
  delete localTermsPr;
  delete jacFF;
  delete gaussFF;
}

Complex FormulationOO2::weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return
    a * localTermsFF->getTerm(dofI, dofJ, elementId) -
    b * localTermsGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOO2::rhs(size_t equationI, size_t elementId) const{
  return localTermsPr->getTerm(equationI, 0, elementId);
}

const FunctionSpace& FormulationOO2::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationOO2::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationOO2::domain(void) const{
  return *ddomain;
}

bool FormulationOO2::isBlock(void) const{
  return true;
}

void FormulationOO2::update(void){
  // Delete RHS (localTermsPr)
  delete localTermsPr;

  // Get DDM Dofs from DDMContext
  const map<Dof, Complex>& ddm = context->getDDMDofs();

  // Pre-evalution
  const fullMatrix<double>& gCFF = gaussFF->getPoints();
  basis->preEvaluateFunctions(gCFF);

  // New RHS
  localTermsPr =
    new TermProjectionField<Complex>(*jacFF, *basis, *gaussFF, *fspace, ddm);
}
