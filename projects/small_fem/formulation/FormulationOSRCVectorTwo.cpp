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
                         Complex R0,
                         const TermGradGrad<double>& localGG,
                         const TermCurlCurl<double>& localCC){
  // Save Data //
  this->oneOverKEpsSquare = Complex(1, 0) / (kEps * kEps);
  this->oneOverR0         = Complex(1, 0) / R0;
  this->ffield            = &field;
  this->ttest             = &test;
  this->ddomain           = &domain;
  this->localGG           = &localGG;
  this->localCC           = &localCC;
}

FormulationOSRCVectorTwo::~FormulationOSRCVectorTwo(void){
}

Complex FormulationOSRCVectorTwo::weak(size_t dofI, size_t dofJ,
                                       size_t elementId) const{

  return
    oneOverR0                     * localGG->getTerm(dofI, dofJ, elementId) -
    oneOverR0 * oneOverKEpsSquare * localCC->getTerm(dofI, dofJ, elementId);
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
