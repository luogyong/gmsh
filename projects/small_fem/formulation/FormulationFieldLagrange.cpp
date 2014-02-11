#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationFieldLagrange.h"

using namespace std;

FormulationFieldLagrange::
FormulationFieldLagrange(const GroupOfElement& domain,
                         const FunctionSpaceScalar& fsField,
                         const FunctionSpaceScalar& fsTest,
                         double (*f)(fullVector<double>& xyz)){
  // Save Domain //
  goe = &domain;

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationFieldLagrange needs a uniform mesh");

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

  GroupOfJacobian jac(domain, gC, "jacobian");

  localTerms      = new TermFieldField(jac, basis, gW);
  projectionTerms = new TermProjectionField<double>(jac, basis, gW, gC, f);
}

FormulationFieldLagrange::~FormulationFieldLagrange(void){
  delete localTerms;
  delete projectionTerms;
}

Complex FormulationFieldLagrange::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  return Complex(localTerms->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationFieldLagrange::rhs(size_t equationI, size_t elementId) const{
  return Complex(projectionTerms->getTerm(0, equationI, elementId), 0);
}

const FunctionSpace& FormulationFieldLagrange::field(void) const{
  return *fsF;
}

const FunctionSpace& FormulationFieldLagrange::test(void) const{
  return *fsT;
}

const GroupOfElement& FormulationFieldLagrange::domain(void) const{
  return *goe;
}
