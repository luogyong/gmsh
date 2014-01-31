#include <complex>

#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationSteadyWaveScalar.h"

using namespace std;

template<>
FormulationSteadyWaveScalar<complex<double> >::
FormulationSteadyWaveScalar(const GroupOfElement& goe,
                            const FunctionSpaceScalar& fs,
                            double k){
  // Save Domain //
  this->goe = &goe;

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = goe.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationSteadyWaveScalar needs a uniform mesh");

  // Wave Squared //
  kSquare = k * k;

  // Save FunctionSpace & Get Basis //
  const Basis& basis = fs.getBasis(eType);
  const size_t order = basis.getOrder();
  fspace             = &fs;

  // Gaussian Quadrature //
  Quadrature gaussGradGrad(eType, order - 1, 2);
  Quadrature gaussFF(eType, order, 2);

  const fullMatrix<double>& gC1 = gaussGradGrad.getPoints();
  const fullVector<double>& gW1 = gaussGradGrad.getWeights();

  const fullMatrix<double>& gC2 = gaussFF.getPoints();
  const fullVector<double>& gW2 = gaussFF.getWeights();

  // Local Terms //
  basis.preEvaluateDerivatives(gC1);
  basis.preEvaluateFunctions(gC2);

  GroupOfJacobian jac1(goe, gC1, "invert");
  GroupOfJacobian jac2(goe, gC2, "jacobian");

  localTerms1 = new TermGradGrad(jac1, basis, gW1);
  localTerms2 = new TermFieldField(jac2, basis, gW2);
}

template<>
FormulationSteadyWaveScalar<complex<double> >::
~FormulationSteadyWaveScalar(void){
  delete localTerms1;
  delete localTerms2;
}

template<>
complex<double> FormulationSteadyWaveScalar<complex<double> >::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return
    complex<double>(localTerms1->getTerm(dofI, dofJ, elementId) -
                    localTerms2->getTerm(dofI, dofJ, elementId) * kSquare, 0);
}

template<>
complex<double> FormulationSteadyWaveScalar<complex<double> >::
rhs(size_t equationI, size_t elementId) const{
  return complex<double>(0, 0);
}

template<>
bool FormulationSteadyWaveScalar<complex<double> >::isGeneral(void) const{
  return false;
}

template<>
complex<double> FormulationSteadyWaveScalar<complex<double> >::
weakB(size_t dofI, size_t dofJ, size_t elementId) const{
  return complex<double>(0, 0);
}

template<>
const FunctionSpace& FormulationSteadyWaveScalar<complex<double> >::
fs(void) const{
  return *fspace;
}

template<>
const GroupOfElement& FormulationSteadyWaveScalar<complex<double> >::
domain(void) const{
  return *goe;
}
