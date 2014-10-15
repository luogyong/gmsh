#include "Exception.h"
#include "FormulationOSRCVectorTwo.h"

using namespace std;

FormulationOSRCVectorTwo::FormulationOSRCVectorTwo(void){
}

FormulationOSRCVectorTwo::
FormulationOSRCVectorTwo(const GroupOfElement& domain,
                         const FunctionSpace& field,
                         const FunctionSpace& test,
                         Complex kEps,
                         const TermGradGrad<double>& localGG,
                         const TermCurlCurl<double>& localCC){
  // Save Data //
  this->minusOneOverKEpsSquare = Complex(-1, 0) / (kEps * kEps);
  this->ffield                 = &field;
  this->ttest                  = &test;
  this->ddomain                = &domain;
  this->localGG                = &localGG;
  this->localCC                = &localCC;
}

FormulationOSRCVectorTwo::~FormulationOSRCVectorTwo(void){
}

Complex FormulationOSRCVectorTwo::weak(size_t dofI, size_t dofJ,
                                       size_t elementId) const{

  return
    localGG->getTerm(dofI, dofJ, elementId) +
    localCC->getTerm(dofI, dofJ, elementId) * minusOneOverKEpsSquare;
}

Complex FormulationOSRCVectorTwo::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCVectorTwo::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorTwo::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationOSRCVectorTwo::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorTwo::isBlock(void) const{
  return true;
}
