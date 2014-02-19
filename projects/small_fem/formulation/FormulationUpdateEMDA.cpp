#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "Exception.h"
#include "FormulationUpdateEMDA.h"

using namespace std;

FormulationUpdateEMDA::
FormulationUpdateEMDA(const GroupOfElement& domain,
                      const FunctionSpaceScalar& fs,
                      double k,
                      double chi,
                      const std::map<Dof, Complex>& sol,
                      const std::map<Dof, Complex>& oldG){

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
  const fullVector<double>& gW = gauss.getWeights();

  // Pre-evalution //
  basis.preEvaluateFunctions(gC);

  // Jacobian //
  GroupOfJacobian jac(domain, gC, "jacobian");

  // Local Terms //
  lGout = new TermFieldField<double>(jac, basis, gW);
  lGin  = new TermProjectionField<Complex>(jac, basis, gW, gC, fs, oldG);
  lU    = new TermProjectionField<Complex>(jac, basis, gW, gC, fs, sol);
}

FormulationUpdateEMDA::~FormulationUpdateEMDA(void){
  delete lGout;
  delete lGin;
  delete lU;
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
  return *goe;
}
