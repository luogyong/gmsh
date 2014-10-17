#include "Exception.h"
#include "FormulationOSRCHelper.h"
#include "FormulationOSRCVectorOne.h"

using namespace std;

FormulationOSRCVectorOne::FormulationOSRCVectorOne(void){
}

FormulationOSRCVectorOne::
FormulationOSRCVectorOne(const GroupOfElement& domain,
                         const FunctionSpace& field,
                         const FunctionSpace& test,
                         double  k,
                         const TermGradGrad<double>& localGG){
  // Save Data //
  this->jK      = Complex(0, k);
  this->ffield  = &field;
  this->ttest   = &test;
  this->ddomain = &domain;
  this->localGG = &localGG;
}

FormulationOSRCVectorOne::~FormulationOSRCVectorOne(void){
}

Complex FormulationOSRCVectorOne::weak(size_t dofI, size_t dofJ,
                                       size_t elementId) const{

  return jK * localGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCVectorOne::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCVectorOne::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorOne::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationOSRCVectorOne::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorOne::isBlock(void) const{
  return true;
}
