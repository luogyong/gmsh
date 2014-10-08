#include "FormulationJFLeeTwo.h"

using namespace std;

FormulationJFLeeTwo::FormulationJFLeeTwo(void){
}

FormulationJFLeeTwo::
FormulationJFLeeTwo(const GroupOfElement& domain,
                    const FunctionSpaceVector& field,
                    const FunctionSpaceVector& test,
                    double k,
                    const TermGradGrad<double>& local){
  // Save Data //
  this->k       = k;
  this->ffield  = &field;
  this->ttest   = &test;
  this->ddomain = &domain;
  this->local   = &local;
}

FormulationJFLeeTwo::~FormulationJFLeeTwo(void){
}

Complex FormulationJFLeeTwo::weak(size_t dofI, size_t dofJ,
                                  size_t elementId) const{
  return Complex(0, -k) * local->getTerm(dofI, dofJ, elementId);
}

Complex FormulationJFLeeTwo::rhs(size_t equationI, size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationJFLeeTwo::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationJFLeeTwo::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationJFLeeTwo::domain(void) const{
  return *ddomain;
}

bool FormulationJFLeeTwo::isBlock(void) const{
  return true;
}
