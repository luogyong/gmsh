#include "Exception.h"
#include "FormulationOSRCVectorThree.h"

using namespace std;

FormulationOSRCVectorThree::FormulationOSRCVectorThree(void){
}

FormulationOSRCVectorThree::
FormulationOSRCVectorThree(const GroupOfElement& domain,
                           const FunctionSpace& field,
                           const FunctionSpace& test,
                           const TermGradGrad<double>& localGG){
  // Save Data //
  this->ffield  = &field;
  this->ttest   = &test;
  this->ddomain = &domain;
  this->localGG = &localGG;
}

FormulationOSRCVectorThree::~FormulationOSRCVectorThree(void){
}

Complex FormulationOSRCVectorThree::weak(size_t dofI, size_t dofJ,
                                         size_t elementId) const{

  return Complex(-1, 0) * localGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCVectorThree::rhs(size_t equationI,
                                        size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCVectorThree::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorThree::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationOSRCVectorThree::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorThree::isBlock(void) const{
  return true;
}
