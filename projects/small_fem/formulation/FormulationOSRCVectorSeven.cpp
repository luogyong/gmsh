#include "Exception.h"
#include "FormulationOSRCVectorSeven.h"

using namespace std;

FormulationOSRCVectorSeven::FormulationOSRCVectorSeven(void){
}

FormulationOSRCVectorSeven::
FormulationOSRCVectorSeven(const GroupOfElement& domain,
                           const FunctionSpace& field,
                           const FunctionSpace& test,
                           Complex Bi,
                           const TermGradGrad<double>& localGG){
  // Save Data //
  this->minusBi = Complex(-1, 0) * Bi;
  this->ffield  = &field;
  this->ttest   = &test;
  this->ddomain = &domain;
  this->localGG = &localGG;
}

FormulationOSRCVectorSeven::~FormulationOSRCVectorSeven(void){
}

Complex FormulationOSRCVectorSeven::weak(size_t dofI, size_t dofJ,
                                         size_t elementId) const{

  return minusBi * localGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCVectorSeven::rhs(size_t equationI,
                                        size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCVectorSeven::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorSeven::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationOSRCVectorSeven::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorSeven::isBlock(void) const{
  return true;
}
