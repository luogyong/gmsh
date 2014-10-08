#include "FormulationJFLeeOne.h"

using namespace std;

FormulationJFLeeOne::FormulationJFLeeOne(void){
}

FormulationJFLeeOne::
FormulationJFLeeOne(const GroupOfElement& domain,
                    const FunctionSpace& field,
                    double k,
                    const TermProjectionGrad<Complex>& localRHS){
  // Save Data //
  this->k        = k;
  this->ffield   = &field;
  this->ddomain  = &domain;
  this->localRHS = &localRHS;
}

FormulationJFLeeOne::~FormulationJFLeeOne(void){
}

Complex FormulationJFLeeOne::weak(size_t dofI, size_t dofJ,
                                  size_t elementId) const{
  return 0;
}

Complex FormulationJFLeeOne::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, -k) * localRHS->getTerm(equationI, 0, elementId);
}

const FunctionSpace& FormulationJFLeeOne::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationJFLeeOne::test(void) const{
  return *ffield;
}

const GroupOfElement& FormulationJFLeeOne::domain(void) const{
  return *ddomain;
}

bool FormulationJFLeeOne::isBlock(void) const{
  return true;
}

void FormulationJFLeeOne::update(TermProjectionGrad<Complex>& localRHS){
  this->localRHS = &localRHS;
}
