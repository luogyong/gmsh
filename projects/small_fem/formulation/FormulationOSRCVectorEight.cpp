#include "Exception.h"
#include "FormulationOSRCVectorEight.h"

using namespace std;

FormulationOSRCVectorEight::FormulationOSRCVectorEight(void){
}

FormulationOSRCVectorEight::
FormulationOSRCVectorEight(const GroupOfElement& domain,
                           const FunctionSpace& field,
                           const FunctionSpace& test,
                           const TermGradGrad<double>& localGG){
  // Save Data //
  this->ffield  = &field;
  this->ttest   = &test;
  this->ddomain = &domain;
  this->localGG = &localGG;
}

FormulationOSRCVectorEight::~FormulationOSRCVectorEight(void){
}

Complex FormulationOSRCVectorEight::weak(size_t dofI, size_t dofJ,
                                         size_t elementId) const{

  return localGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCVectorEight::rhs(size_t equationI,
                                        size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCVectorEight::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorEight::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationOSRCVectorEight::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorEight::isBlock(void) const{
  return true;
}
