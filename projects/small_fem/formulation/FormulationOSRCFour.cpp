#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCFour.h"

using namespace std;

FormulationOSRCFour::FormulationOSRCFour(void){
}

FormulationOSRCFour::
FormulationOSRCFour(const GroupOfElement& domain,
                    const FunctionSpace& field,
                    const FunctionSpace& auxiliary,
                    const TermFieldField<double>& localTerm){
  // Save Data //
  this->ffield    = &field;
  this->faux      = &auxiliary;
  this->ddomain   = &domain;
  this->localTerm = &localTerm;
}

FormulationOSRCFour::~FormulationOSRCFour(void){
}

Complex FormulationOSRCFour::weak(size_t dofI, size_t dofJ,
                                  size_t elementId) const{
  return
    Complex(-1, 0) * localTerm->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCFour::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCFour::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCFour::test(void) const{
  return *faux;
}

const GroupOfElement& FormulationOSRCFour::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCFour::isBlock(void) const{
  return true;
}
