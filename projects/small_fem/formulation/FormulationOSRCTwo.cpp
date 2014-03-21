#include <cmath>

#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCTwo.h"

using namespace std;

FormulationOSRCTwo::FormulationOSRCTwo(const GroupOfElement& domain,
                                       const FunctionSpaceScalar& fField,
                                       const FunctionSpaceScalar& fTest,
                                       double  k,
                                       Complex keps){

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOSRCTwo needs a uniform mesh");

  // Wavenumber (normal and complexified) //
  this->k    = k;
  this->keps = keps;

  // Pade A1 //
  A1 = padeAj(1, 1, M_PI / 4.);

  // FunctionSpace (test and field) and Domain //
  ffField = &fField;
  ffTest  = &fTest;
  ddomain = &domain;

  // Basis (for test functions) //
  const Basis& basis = fTest.getBasis(eType);

  // Quadrature //
  Quadrature gauss(eType, basis.getOrder() - 1, 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Basis derivatives pre evaluation //
  basis.preEvaluateDerivatives(gC);

  // Jacobians //
  GroupOfJacobian jac(domain, gC, "invert");

  // Local Term //
  local = new TermGradGrad<double>(jac, basis, gauss);
}

FormulationOSRCTwo::~FormulationOSRCTwo(void){
  delete local;
}

Complex FormulationOSRCTwo::weak(size_t dofI, size_t dofJ,
                                 size_t elementId) const{
  return
    Complex(0, k) * A1 / (keps * keps) * local->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCTwo::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCTwo::field(void) const{
  return *ffField;
}

const FunctionSpace& FormulationOSRCTwo::test(void) const{
  return *ffTest;
}

const GroupOfElement& FormulationOSRCTwo::domain(void) const{
  return *ddomain;
}

double FormulationOSRCTwo::pade_aj(int j, int N){
  return 2. / (2. * N + 1.) * sqrt(sin((double)j * M_PI / (2. * N + 1.)));
}

double FormulationOSRCTwo::pade_bj(int j, int N){
  return sqrt(cos((double)j * M_PI / (2. *N + 1.)));
}

Complex FormulationOSRCTwo::padeAj(int j, int N, double theta){
  Complex one = Complex(1, 0);
  Complex res;
  Complex z;

  z   = Complex(cos(-theta / 2.), sin(-theta / 2.));
  res = z * pade_aj(j, N);

  z   = Complex(cos(-theta) - 1., sin(-theta));
  res = res / ((one + z * pade_bj(j, N)) * (one + z * pade_bj(j, N)));

  return res;
}
