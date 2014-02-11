#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationLagrangeField.h"

using namespace std;

FormulationLagrangeField::
FormulationLagrangeField(const GroupOfElement& goe,
                         const FunctionSpaceScalar& fsField,
                         const FunctionSpaceScalar& fsTest){
  // Save Domain //
  this->goe = &goe;

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = goe.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationLagrangeField needs a uniform mesh");

  // Save FunctionSpace & Get Basis //
  const Basis& basis = fsTest.getBasis(eType);
  const size_t order = basis.getOrder();
  fsF                = &fsField;
  fsT                = &fsTest;

  // Gaussian Quadrature //
  Quadrature gaussFF(eType, order, 2);
  const fullMatrix<double>& gC = gaussFF.getPoints();
  const fullVector<double>& gW = gaussFF.getWeights();

  // Local Terms //
  basis.preEvaluateFunctions(gC);

  GroupOfJacobian jac(goe, gC, "jacobian");

  localTerms = new TermFieldField(jac, basis, gW);
}

FormulationLagrangeField::~FormulationLagrangeField(void){
  delete localTerms;
}

complex<double> FormulationLagrangeField::weak(size_t dofI, size_t dofJ,
                                               size_t elementId) const{
  return
    complex<double>(localTerms->getTerm(dofI, dofJ, elementId), 0);
}

complex<double> FormulationLagrangeField::rhs(size_t equationI,
                                              size_t elementId) const{
  return complex<double>(0, 0);
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
