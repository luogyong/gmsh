#include "FormulationJFLeeEight.h"

using namespace std;

FormulationJFLeeEight::FormulationJFLeeEight(void){
}

FormulationJFLeeEight::
FormulationJFLeeEight(const GroupOfElement& domain,
                      const FunctionSpaceVector& field,
                      const FunctionSpaceScalar& test,
                      const TermGradGrad<double>& local){
  // Save Data //
  this->ffield  = &field;
  this->ttest   = &test;
  this->ddomain = &domain;
  this->local   = &local;
}

FormulationJFLeeEight::~FormulationJFLeeEight(void){
}

Complex FormulationJFLeeEight::weak(size_t dofI, size_t dofJ,
                                    size_t elementId) const{
  return local->getTerm(dofI, dofJ, elementId);
}

Complex FormulationJFLeeEight::rhs(size_t equationI, size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationJFLeeEight::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationJFLeeEight::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationJFLeeEight::domain(void) const{
  return *ddomain;
}

bool FormulationJFLeeEight::isBlock(void) const{
  return true;
}
