#include "Exception.h"
#include "FormulationOSRCVectorNine.h"

using namespace std;

FormulationOSRCVectorNine::FormulationOSRCVectorNine(void){
}

FormulationOSRCVectorNine::
FormulationOSRCVectorNine(const GroupOfElement& domain,
                          const FunctionSpace& field,
                          const FunctionSpace& test,
                          Complex kEps,
                          Complex Ai,
                          double  k,
                          const TermCurlCurl<double>& localCC){
  // Save Data //
  this->jOverK                 = Complex( 0, 1. / k);
  this->minusOneOverKEpsSquare = Complex(-1, 0) / (kEps * kEps);
  this->Ai                     = Ai;
  this->ffield                 = &field;
  this->ttest                  = &test;
  this->ddomain                = &domain;
  this->localCC                = &localCC;
}

FormulationOSRCVectorNine::~FormulationOSRCVectorNine(void){
}

Complex FormulationOSRCVectorNine::weak(size_t dofI, size_t dofJ,
                                        size_t elementId) const{

  return
    jOverK * Ai * minusOneOverKEpsSquare
                * localCC->getTerm(dofI, dofJ, elementId);
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
