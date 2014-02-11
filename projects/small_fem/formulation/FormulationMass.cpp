#include "SmallFem.h"
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
rhs(size_t equationI, size_t elementId) const{

  return 0;
}

// Complex Version //
// --------------- //
template<>
Complex FormulationMass<Complex>::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  return Complex(localTerms->getTerm(dofI, dofJ, elementId), 0);
}

template<>
Complex FormulationMass<Complex>::
rhs(size_t equationI, size_t elementId) const{

  return Complex(0, 0);
}
