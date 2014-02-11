#include "FormulationProjectionVector.h"

// Real Version //
// ------------ //
template<>
double FormulationProjectionVector<double>::weak(size_t dofI, size_t dofJ,
                                                 size_t elementId) const{

  return localTerms1->getTerm(dofI, dofJ, elementId);
}

template<>
double FormulationProjectionVector<double>::rhs(size_t equationI,
                                                size_t elementId) const{

  return localTerms2->getTerm(0, equationI, elementId);
}

// Complex Version //
// --------------- //
template<>
Complex FormulationProjectionVector<Complex>::weak(size_t dofI, size_t dofJ,
                                                   size_t elementId) const{

  return Complex(localTerms1->getTerm(dofI, dofJ, elementId), 0);
}

template<>
Complex FormulationProjectionVector<Complex>::rhs(size_t equationI,
                                                  size_t elementId) const{

  return localTerms2->getTerm(0, equationI, elementId);
}
