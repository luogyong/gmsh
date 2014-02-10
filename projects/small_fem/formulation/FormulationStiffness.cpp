#include <complex>
#include "FormulationStiffness.h"

using namespace std;

// Real Version //
// ------------ //
template<>
double FormulationStiffness<double>::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  return localTerms->getTerm(dofI, dofJ, elementId);
}

template<>
double FormulationStiffness<double>::
weakB(size_t dofI, size_t dofJ, size_t elementId) const{

  return 0;
}

template<>
double FormulationStiffness<double>::
rhs(size_t equationI, size_t elementId) const{

  return 0;
}

// Complex Version //
// --------------- //
template<>
complex<double> FormulationStiffness<complex<double> >::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  return complex<double>(localTerms->getTerm(dofI, dofJ, elementId), 0);
}

template<>
complex<double> FormulationStiffness<complex<double> >::
weakB(size_t dofI, size_t dofJ, size_t elementId) const{

  return complex<double>(0, 0);
}

template<>
complex<double> FormulationStiffness<complex<double> >::
rhs(size_t equationI, size_t elementId) const{

  return complex<double>(0, 0);
}
