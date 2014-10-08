#include "FormulationJFLeeThree.h"

using namespace std;

FormulationJFLeeThree::FormulationJFLeeThree(void){
}

FormulationJFLeeThree::
FormulationJFLeeThree(const GroupOfElement& domain,
                      const FunctionSpaceScalar& field,
                      const FunctionSpaceVector& test,
                      Complex C2,
                      const TermGradGrad<double>& local){
  // Save Data //
  this->C2      = C2;
  this->ffield  = &field;
  this->ttest   = &test;
  this->ddomain = &domain;
  this->local   = &local;
}

FormulationJFLeeThree::~FormulationJFLeeThree(void){
}

Complex FormulationJFLeeThree::weak(size_t dofI, size_t dofJ,
                                    size_t elementId) const{
  return Complex(-1, 0) * C2 * local->getTerm(dofI, dofJ, elementId);
}

Complex FormulationJFLeeThree::rhs(size_t equationI, size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationJFLeeThree::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationJFLeeThree::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationJFLeeThree::domain(void) const{
  return *ddomain;
}

bool FormulationJFLeeThree::isBlock(void) const{
  return true;
}
