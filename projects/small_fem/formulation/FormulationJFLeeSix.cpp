#include "FormulationJFLeeSix.h"

using namespace std;

FormulationJFLeeSix::FormulationJFLeeSix(void){
}

FormulationJFLeeSix::
FormulationJFLeeSix(const GroupOfElement& domain,
                    const FunctionSpaceVector& field,
                    const FunctionSpaceVector& test,
                    Complex C1,
                    const TermCurlCurl<double>& local){
  // Save Data //
  this->C1      = C1;
  this->ffield  = &field;
  this->ttest   = &test;
  this->ddomain = &domain;
  this->local   = &local;
}

FormulationJFLeeSix::~FormulationJFLeeSix(void){
}

Complex FormulationJFLeeSix::weak(size_t dofI, size_t dofJ,
                                    size_t elementId) const{
  return C1 * local->getTerm(dofI, dofJ, elementId);
}

Complex FormulationJFLeeSix::rhs(size_t equationI, size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationJFLeeSix::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationJFLeeSix::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationJFLeeSix::domain(void) const{
  return *ddomain;
}

bool FormulationJFLeeSix::isBlock(void) const{
  return true;
}
