#include "FormulationPML.h"

FormulationPML::FormulationPML(const GroupOfElement& domain,
                               const FunctionSpace& fs,
                               double k){
  // Wavenumber //
  this->k = k;

  // Steady wave formulations //
  wave = new FormulationSteadyWave<Complex>(domain, fs, k);
}

FormulationPML::~FormulationPML(void){
  delete wave;
}

Complex FormulationPML::weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return Complex(0, -1 / k) * wave->weak(dofI, dofJ, elementId);
}

Complex FormulationPML::rhs(size_t equationI, size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationPML::field(void) const{
  return wave->field();
}

const FunctionSpace& FormulationPML::test(void) const{
  return wave->test();
}

const GroupOfElement& FormulationPML::domain(void) const{
  return wave->domain();
}
