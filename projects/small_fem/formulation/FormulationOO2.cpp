#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationOO2.h"

using namespace std;

FormulationOO2::FormulationOO2(DDMContext& context){
  // Check if OO2 DDMContext //
  if(context.getType() != DDMContext::typeOO2)
    throw Exception("FormulationOO2 needs a OO2 DDMContext");

  // Get Domain and FunctionSpace from DDMContext //
  fspace  = &context.getFunctionSpace();
  ddomain = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOO2 needs a uniform mesh");

  // a & b //
  this->a = context.OO2_A;
  this->b = context.OO2_B;

  // Get Basis //
  const Basis& basis = fspace->getBasis(eType);

  // Gaussian Quadrature (Field - Field & Field - Projection) //
  Quadrature gaussFF(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gCFF = gaussFF.getPoints();

  // Gaussian Quadrature (Grad - Grad) //
  Quadrature gaussGG(eType, basis.getOrder() - 1, 2);
  const fullMatrix<double>& gCGG = gaussGG.getPoints();

  // Pre-evalution //
  basis.preEvaluateFunctions(gCFF);
  basis.preEvaluateDerivatives(gCGG);

  // Jacobians //
  GroupOfJacobian jacFF(*ddomain, gCFF, "jacobian");
  GroupOfJacobian jacGG(*ddomain, gCGG, "invert");

  // Get DDM Dofs from DDMContext //
  const map<Dof, Complex>& ddm = context.getDDMDofs();

  // Local Terms //
  localTermsFF = new TermFieldField<double>(jacFF, basis, gaussFF);
  localTermsGG = new TermGradGrad<double>(jacGG, basis, gaussGG);
  localTermsPr =
    new TermProjectionField<Complex>(jacFF, basis, gaussFF, *fspace, ddm);
}

FormulationOO2::~FormulationOO2(void){
  delete localTermsFF;
  delete localTermsGG;
  delete localTermsPr;
}

Complex FormulationOO2::weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return
    a * localTermsFF->getTerm(dofI, dofJ, elementId) -
    b * localTermsGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOO2::rhs(size_t equationI, size_t elementId) const{
  return localTermsPr->getTerm(0, equationI, elementId);
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
