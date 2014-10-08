#include "FormulationJFLeeSeven.h"

using namespace std;

FormulationJFLeeSeven::FormulationJFLeeSeven(void){
}

FormulationJFLeeSeven::
FormulationJFLeeSeven(const GroupOfElement& domain,
                      const FunctionSpaceScalar& field,
                      const TermFieldField<double>& local){
  // Save Data //
  this->ffield  = &field;
  this->ddomain = &domain;
  this->local   = &local;
}

FormulationJFLeeSeven::~FormulationJFLeeSeven(void){
}

Complex FormulationJFLeeSeven::weak(size_t dofI, size_t dofJ,
                                    size_t elementId) const{
  return local->getTerm(dofI, dofJ, elementId);
}

Complex FormulationJFLeeSeven::rhs(size_t equationI, size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationJFLeeSeven::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationJFLeeSeven::test(void) const{
  return *ffield;
}

const GroupOfElement& FormulationJFLeeSeven::domain(void) const{
  return *ddomain;
}

bool FormulationJFLeeSeven::isBlock(void) const{
  return true;
}
