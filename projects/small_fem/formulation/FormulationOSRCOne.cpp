#include <cmath>

#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCOne.h"

using namespace std;

FormulationOSRCOne::FormulationOSRCOne(const GroupOfElement& domain,
                                       const FunctionSpaceScalar& fs,
                                       double k,
                                       const map<Dof, Complex>& ddmDof){

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOSRCOne needs a uniform mesh");

  // Wavenumber //
  this->k = k;

  // Pade C0
  C0 = padeC0(1, M_PI / 4.);

  // Function Space & Domain //
  fspace  = &fs;
  ddomain = &domain;

  // Basis //
  const Basis& basis = fs.getBasis(eType);

  // Quadrature //
  Quadrature gauss(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Basis pre evaluation //
  basis.preEvaluateFunctions(gC);

  // Jacobians //
  GroupOfJacobian jac(domain, gC, "jacobian");

  // Local Terms //
  localLHS = new TermFieldField<double>(jac, basis, gauss);
  localRHS = new TermProjectionField<Complex>(jac, basis, gauss, fs, ddmDof);
}

FormulationOSRCOne::~FormulationOSRCOne(void){
  delete localLHS;
  delete localRHS;
}

Complex FormulationOSRCOne::weak(size_t dofI, size_t dofJ,
                                 size_t elementId) const{

  return Complex(0, -k) * C0 * localLHS->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCOne::rhs(size_t equationI, size_t elementId) const{
  return localRHS->getTerm(0, equationI, elementId);
}

const FunctionSpace& FormulationOSRCOne::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationOSRCOne::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationOSRCOne::domain(void) const{
  return *ddomain;
}

double FormulationOSRCOne::pade_aj(int j, int N){
  return 2. / (2. * N + 1.) * sqrt(sin((double)j * M_PI / (2. * N + 1.)));
}

double FormulationOSRCOne::pade_bj(int j, int N){
  return sqrt(cos((double)j * M_PI / (2. *N + 1.)));
}

Complex FormulationOSRCOne::padeC0(int N, double theta){
  Complex sum = Complex(1, 0);
  Complex one = Complex(1, 0);
  Complex z   = Complex(cos(-theta) - 1,  sin(-theta));

  for(int j = 1; j <= N; j++)
    sum += (z * pade_aj(j, N)) / (one + z * pade_bj(j, N));

  z = Complex(cos(theta / 2.), sin(theta / 2.));

  return sum * z;
}
