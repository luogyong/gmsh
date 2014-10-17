#include "Exception.h"
#include "FormulationOSRCVectorFour.h"

using namespace std;

FormulationOSRCVectorFour::FormulationOSRCVectorFour(void){
}

FormulationOSRCVectorFour::
FormulationOSRCVectorFour(const GroupOfElement& domain,
                          const FunctionSpace& field,
                          const FunctionSpace& test,
                          Complex R0,
                          Complex Ai,
                          Complex Bi,
                          const TermGradGrad<double>& localGG){
  // Save Data //
  this->alpha   = Complex(-1, 0) * Ai / (Bi * R0);
  this->ffield  = &field;
  this->ttest   = &test;
  this->ddomain = &domain;
  this->localGG = &localGG;
}

FormulationOSRCVectorFour::~FormulationOSRCVectorFour(void){
}

Complex FormulationOSRCVectorFour::weak(size_t dofI, size_t dofJ,
                                        size_t elementId) const{

  return alpha * localGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCVectorFour::rhs(size_t equationI,
                                       size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCVectorFour::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorFour::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationOSRCVectorFour::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorFour::isBlock(void) const{
  return true;
}
