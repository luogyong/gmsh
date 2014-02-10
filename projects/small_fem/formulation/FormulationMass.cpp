#include <complex>
#include "FormulationMass.h"

using namespace std;

// Real Version //
// ------------ //
template<>
double FormulationMass<double>::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  return localTerms->getTerm(dofI, dofJ, elementId);
}

template<>
double FormulationMass<double>::
weakB(size_t dofI, size_t dofJ, size_t elementId) const{

  return 0;
}

template<>
double FormulationMass<double>::
rhs(size_t equationI, size_t elementId) const{

  return 0;
}

// Complex Version //
// --------------- //
template<>
complex<double> FormulationMass<complex<double> >::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  return complex<double>(localTerms->getTerm(dofI, dofJ, elementId), 0);
}

template<>
complex<double> FormulationMass<complex<double> >::
weakB(size_t dofI, size_t dofJ, size_t elementId) const{

  return complex<double>(0, 0);
}

template<>
complex<double> FormulationMass<complex<double> >::
rhs(size_t equationI, size_t elementId) const{

  return complex<double>(0, 0);
}
