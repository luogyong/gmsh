#include "FormulationJFLeeFive.h"

using namespace std;

FormulationJFLeeFive::FormulationJFLeeFive(void){
}

FormulationJFLeeFive::
FormulationJFLeeFive(const GroupOfElement& domain,
                     const FunctionSpaceVector& field,
                     const FunctionSpaceVector& test,
                     double k,
                     const TermGradGrad<double>& local){
  // Save Data //
  this->kSquare = k * k;
  this->ffield  = &field;
  this->ttest   = &test;
  this->ddomain = &domain;
  this->local   = &local;
}

FormulationJFLeeFive::~FormulationJFLeeFive(void){
}

Complex FormulationJFLeeFive::weak(size_t dofI, size_t dofJ,
                                    size_t elementId) const{
  return Complex(-kSquare * local->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationJFLeeFive::rhs(size_t equationI, size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationJFLeeFive::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationJFLeeFive::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationJFLeeFive::domain(void) const{
  return *ddomain;
}

bool FormulationJFLeeFive::isBlock(void) const{
  return true;
}
