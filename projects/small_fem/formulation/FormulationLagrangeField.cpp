#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationLagrangeField.h"

using namespace std;

FormulationLagrangeField::
FormulationLagrangeField(const GroupOfElement& domain,
                         const FunctionSpaceScalar& field,
                         const FunctionSpaceScalar& test){
  // Save Domain //
  goe = &domain;

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationLagrangeField needs a uniform mesh");

  // Save FunctionSpace & Get Basis //
  const Basis& basis = test.getBasis(eType);
  const size_t order = basis.getOrder();
  fsF                = &field;
  fsT                = &test;

  // Gaussian Quadrature //
  Quadrature gaussFF(eType, order, 2);
  const fullMatrix<double>& gC = gaussFF.getPoints();

  // Local Terms //
  basis.preEvaluateFunctions(gC);

  GroupOfJacobian jac(domain, gC, "jacobian");

  localTerms = new TermFieldField<double>(jac, basis, gaussFF);
}

FormulationLagrangeField::~FormulationLagrangeField(void){
  delete localTerms;
}

Complex FormulationLagrangeField::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return Complex(localTerms->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationLagrangeField::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationLagrangeField::field(void) const{
  return *fsF;
}

const FunctionSpace& FormulationLagrangeField::test(void) const{
  return *fsT;
}

const GroupOfElement& FormulationLagrangeField::domain(void) const{
  return *goe;
}
