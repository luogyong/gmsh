#include "SmallFem.h"
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
rhs(size_t equationI, size_t elementId) const{

  return 0;
}

// Complex Version //
// --------------- //
template<>
Complex FormulationStiffness<Complex>::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  return Complex(localTerms->getTerm(dofI, dofJ, elementId), 0);
}

template<>
Complex FormulationStiffness<Complex>::
rhs(size_t equationI, size_t elementId) const{

  return Complex(0, 0);
}
