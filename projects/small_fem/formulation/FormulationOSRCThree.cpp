#include <cmath>

#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCThree.h"

using namespace std;

FormulationOSRCThree::FormulationOSRCThree(const GroupOfElement& domain,
                                           const FunctionSpaceScalar& fspace,
                                           Complex keps){

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOSRCThree needs a uniform mesh");

  // Complexified wavenumber //
  this->keps = keps;

  // Pade B1 //
  B1 = padeBj(1, 1, M_PI / 4.);

  // Function Space & Domain //
  ffspace = &fspace;
  ddomain = &domain;

  // Basis //
  const Basis& basis = fspace.getBasis(eType);

  // Quadrature //
  Quadrature gaussFF(eType, basis.getOrder(), 2);
  Quadrature gaussGG(eType, basis.getOrder() - 1, 2);

  const fullMatrix<double>& gCFF = gaussFF.getPoints();
  const fullMatrix<double>& gCGG = gaussGG.getPoints();

  // Basis pre evaluation //
  basis.preEvaluateFunctions(gCFF);
  basis.preEvaluateDerivatives(gCGG);

  // Jacobians //
  GroupOfJacobian jacFF(domain, gCFF, "jacobian");
  GroupOfJacobian jacGG(domain, gCGG, "invert");

  // Local Terms //
  localFF = new TermFieldField<double>(jacFF, basis, gaussFF);
  localGG = new TermGradGrad<double>  (jacGG, basis, gaussGG);
}

FormulationOSRCThree::~FormulationOSRCThree(void){
  delete localFF;
  delete localGG;
}

Complex FormulationOSRCThree::weak(size_t dofI, size_t dofJ,
                                   size_t elementId) const{
  return
     localFF->getTerm(dofI, dofJ, elementId)
    -localGG->getTerm(dofI, dofJ, elementId) * B1 / (keps * keps);
}

Complex FormulationOSRCThree::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCThree::field(void) const{
  return *ffspace;
}

const FunctionSpace& FormulationOSRCThree::test(void) const{
  return *ffspace;
}

const GroupOfElement& FormulationOSRCThree::domain(void) const{
  return *ddomain;
}

double FormulationOSRCThree::pade_aj(int j, int N){
  return 2. / (2. * N + 1.) * sqrt(sin((double)j * M_PI / (2. * N + 1.)));
}

double FormulationOSRCThree::pade_bj(int j, int N){
  return sqrt(cos((double)j * M_PI / (2. *N + 1.)));
}

Complex FormulationOSRCThree::padeBj(int j, int N, double theta){
  Complex one = Complex(1, 0);
  Complex res;
  Complex z;

  z   = Complex(cos(-theta), sin(-theta));
  res = z * pade_bj(j, N);

  z   = Complex(cos(-theta) - 1., sin(-theta));
  res = res / (one + z * pade_bj(j, N));

  return res;
}
