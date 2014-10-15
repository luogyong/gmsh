#include "Exception.h"
#include "FormulationOSRCVectorSeven.h"

using namespace std;

FormulationOSRCVectorSeven::FormulationOSRCVectorSeven(void){
}

FormulationOSRCVectorSeven::
FormulationOSRCVectorSeven(const GroupOfElement& domain,
                           const FunctionSpace& field,
                           const TermFieldField<double>& localFF){
  // Save Data //
  this->ffield  = &field;
  this->ddomain = &domain;
  this->localFF = &localFF;
}

FormulationOSRCVectorSeven::~FormulationOSRCVectorSeven(void){
}

Complex FormulationOSRCVectorSeven::weak(size_t dofI, size_t dofJ,
                                       size_t elementId) const{

  return localFF->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCVectorSeven::rhs(size_t equationI,
                                      size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCVectorSeven::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorSeven::test(void) const{
  return *ffield;
}

const GroupOfElement& FormulationOSRCVectorSeven::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorSeven::isBlock(void) const{
  return true;
}
