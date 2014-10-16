#include "Exception.h"
#include "FormulationOSRCVectorTen.h"

using namespace std;

FormulationOSRCVectorTen::FormulationOSRCVectorTen(void){
}

FormulationOSRCVectorTen::
FormulationOSRCVectorTen(const GroupOfElement& domain,
                         const FunctionSpace& field,
                         const FunctionSpace& test,
                         Complex kEps,
                         Complex Ai,
                         double  k,
                         const TermGradGrad<double>& localGG){
  // Save Data //
  this->jOverK                = Complex(0, 1. / k);
  this->plusOneOverKEpsSquare = Complex(1, 0) / (kEps * kEps);
  this->Ai                    = Ai;
  this->ffield                = &field;
  this->ttest                 = &test;
  this->ddomain               = &domain;
  this->localGG               = &localGG;
}

FormulationOSRCVectorTen::~FormulationOSRCVectorTen(void){
}

Complex FormulationOSRCVectorTen::weak(size_t dofI, size_t dofJ,
                                       size_t elementId) const{

  return
    jOverK * Ai * plusOneOverKEpsSquare
                * localGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCVectorTen::rhs(size_t equationI,
                                      size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCVectorTen::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorTen::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationOSRCVectorTen::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorTen::isBlock(void) const{
  return true;
}
