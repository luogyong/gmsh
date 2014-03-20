#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationSommerfeld.h"

using namespace std;

FormulationSommerfeld::FormulationSommerfeld(const GroupOfElement& domain,
                                             const FunctionSpaceScalar& fs,
                                             double k){
  // Save Domain //
  goe = &domain;

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationSommerfeld needs a uniform mesh");

  // Wavenumber //
  this->k = k;

  // Save FunctionSpace & Get Basis //
  const Basis& basis = fs.getBasis(eType);
  const size_t order = basis.getOrder();
  fspace             = &fs;

  // Gaussian Quadrature //
  Quadrature gaussFF(eType, order, 2);
  const fullMatrix<double>& gC = gaussFF.getPoints();

  // Local Terms //
  basis.preEvaluateFunctions(gC);

  GroupOfJacobian jac(domain, gC, "jacobian");

  localTerms = new TermFieldField<double>(jac, basis, gaussFF);
}

FormulationSommerfeld::~FormulationSommerfeld(void){
  delete localTerms;
}

Complex FormulationSommerfeld::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return Complex(0, -1 * k * localTerms->getTerm(dofI, dofJ, elementId));
}

Complex FormulationSommerfeld::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationSommerfeld::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationSommerfeld::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationSommerfeld::domain(void) const{
  return *goe;
}
