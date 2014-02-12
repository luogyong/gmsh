#include "SmallFem.h"
#include "FormulationProjection.h"

// Real Version //
// ------------ //
template<>
double FormulationProjection<double>::weak(size_t dofI, size_t dofJ,
                                           size_t elementId) const{

  return localTerms1->getTerm(dofI, dofJ, elementId);
}

template<>
double FormulationProjection<double>::rhs(size_t equationI,
                                          size_t elementId) const{

  return localTerms2->getTerm(0, equationI, elementId);
}

// Complex Version //
// --------------- //
template<>
Complex FormulationProjection<Complex>::weak(size_t dofI, size_t dofJ,
                                             size_t elementId) const{

  return Complex(localTerms1->getTerm(dofI, dofJ, elementId), 0);
}

template<>
Complex FormulationProjection<Complex>::rhs(size_t equationI,
                                            size_t elementId) const{

  return localTerms2->getTerm(0, equationI, elementId);
}
