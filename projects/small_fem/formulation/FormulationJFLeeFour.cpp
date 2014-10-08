#include "FormulationJFLeeFour.h"

using namespace std;

FormulationJFLeeFour::FormulationJFLeeFour(void){
}

FormulationJFLeeFour::
FormulationJFLeeFour(const GroupOfElement& domain,
                      const FunctionSpaceVector& field,
                      double k,
                      const TermGradGrad<double>& local){
  // Save Data //
  this->kSquare = k * k;
  this->ffield  = &field;
  this->ddomain = &domain;
  this->local   = &local;
}

FormulationJFLeeFour::~FormulationJFLeeFour(void){
}

Complex FormulationJFLeeFour::weak(size_t dofI, size_t dofJ,
                                    size_t elementId) const{
  return Complex(kSquare * local->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationJFLeeFour::rhs(size_t equationI, size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationJFLeeFour::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationJFLeeFour::test(void) const{
  return *ffield;
}

const GroupOfElement& FormulationJFLeeFour::domain(void) const{
  return *ddomain;
}

bool FormulationJFLeeFour::isBlock(void) const{
  return true;
}
