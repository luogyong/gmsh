#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCScalarFour.h"

using namespace std;

FormulationOSRCScalarFour::FormulationOSRCScalarFour(void){
}

FormulationOSRCScalarFour::
FormulationOSRCScalarFour(const GroupOfElement& domain,
                          const FunctionSpace& field,
                          const FunctionSpace& auxiliary,
                          const TermFieldField<double>& localTerm){
  // Save Data //
  this->ffield    = &field;
  this->faux      = &auxiliary;
  this->ddomain   = &domain;
  this->localTerm = &localTerm;
}

FormulationOSRCScalarFour::~FormulationOSRCScalarFour(void){
}

Complex FormulationOSRCScalarFour::weak(size_t dofI, size_t dofJ,
                                        size_t elementId) const{
  return
    Complex(-1, 0) * localTerm->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCScalarFour::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCScalarFour::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCScalarFour::test(void) const{
  return *faux;
}

const GroupOfElement& FormulationOSRCScalarFour::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCScalarFour::isBlock(void) const{
  return true;
}
