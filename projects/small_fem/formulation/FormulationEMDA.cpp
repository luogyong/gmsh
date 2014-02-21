#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationEMDA.h"

using namespace std;

FormulationEMDA::FormulationEMDA(const GroupOfElement& domain,
                                 const FunctionSpaceScalar& fs,
                                 double k,
                                 double chi,
                                 const std::map<Dof, Complex>& ddmDof){

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationEMDA needs a uniform mesh");

  // Wavenumber & Chi //
  this->k   = k;
  this->chi = chi;

  // Save FunctionSpace & Domain //
  fspace = &fs;
  goe    = &domain;

  // Basis //
  const Basis& basis = fs.getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Pre-evalution //
  basis.preEvaluateFunctions(gC);

  // Jacobian //
  GroupOfJacobian jac(domain, gC, "jacobian");

  // Local Terms //
  localLHS = new TermFieldField<double>(jac, basis, gauss);
  localRHS = new TermProjectionField<Complex>(jac, basis, gauss, fs, ddmDof);
}

FormulationEMDA::~FormulationEMDA(void){
  delete localLHS;
  delete localRHS;
}

Complex FormulationEMDA::weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return Complex(chi, -k) * localLHS->getTerm(dofI, dofJ, elementId);
}

Complex FormulationEMDA::rhs(size_t equationI, size_t elementId) const{
  return localRHS->getTerm(0, equationI, elementId);
}

const FunctionSpace& FormulationEMDA::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationEMDA::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationEMDA::domain(void) const{
  return *goe;
}
