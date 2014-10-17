#include "Exception.h"
#include "FormulationOSRCVectorNine.h"

using namespace std;

FormulationOSRCVectorNine::FormulationOSRCVectorNine(void){
}

FormulationOSRCVectorNine::
FormulationOSRCVectorNine(const GroupOfElement& domain,
                          const FunctionSpace& field,
                          const FunctionSpace& test,
                          const TermGradGrad<double>& localGG){
  // Save Data //
  this->ffield  = &field;
  this->ttest   = &test;
  this->ddomain = &domain;
  this->localGG = &localGG;
}

FormulationOSRCVectorNine::~FormulationOSRCVectorNine(void){
}

Complex FormulationOSRCVectorNine::weak(size_t dofI, size_t dofJ,
                                        size_t elementId) const{

  return localGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCVectorNine::rhs(size_t equationI,
                                       size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCVectorNine::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorNine::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationOSRCVectorNine::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorNine::isBlock(void) const{
  return true;
}
