#include "FormulationProjectionScalar.h"

// Real Version //
// ------------ //
template<>
double FormulationProjectionScalar<double>::weak(size_t dofI, size_t dofJ,
                                                 size_t elementId) const{

  return localTerms1->getTerm(dofI, dofJ, elementId);
}

template<>
double FormulationProjectionScalar<double>::rhs(size_t equationI,
                                                size_t elementId) const{

  return localTerms2->getTerm(0, equationI, elementId);
}

// Complex Version //
// --------------- //
template<>
Complex FormulationProjectionScalar<Complex>::weak(size_t dofI, size_t dofJ,
                                                   size_t elementId) const{

  return Complex(localTerms1->getTerm(dofI, dofJ, elementId), 0);
}

template<>
Complex FormulationProjectionScalar<Complex>::rhs(size_t equationI,
                                                  size_t elementId) const{

  return localTerms2->getTerm(0, equationI, elementId);
}
