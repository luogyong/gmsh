#include "Exception.h"
#include "FormulationOSRCVectorFour.h"

using namespace std;

FormulationOSRCVectorFour::FormulationOSRCVectorFour(void){
}

FormulationOSRCVectorFour::
FormulationOSRCVectorFour(const GroupOfElement& domain,
                          const FunctionSpace& field,
                          Complex kEps,
                          Complex Bin,
                          const TermGradGrad<double>& localGG,
                          const TermCurlCurl<double>& localCC){
  // Save Data //
  this->minusOneOverKEpsSquare = Complex(-1, 0) / (kEps * kEps);
  this->Bi                     = Bi;
  this->ffield                 = &field;
  this->ddomain                = &domain;
  this->localGG                = &localGG;
  this->localCC                = &localCC;
}

FormulationOSRCVectorFour::~FormulationOSRCVectorFour(void){
}

Complex FormulationOSRCVectorFour::weak(size_t dofI, size_t dofJ,
                                        size_t elementId) const{
  return
    localCC->getTerm(dofI, dofJ, elementId) * Bi * minusOneOverKEpsSquare +
    localGG->getTerm(dofI, dofJ, elementId);

}

Complex FormulationOSRCVectorFour::rhs(size_t equationI,
                                       size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCVectorFour::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorFour::test(void) const{
  return *ffield;
}

const GroupOfElement& FormulationOSRCVectorFour::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorFour::isBlock(void) const{
  return true;
}
