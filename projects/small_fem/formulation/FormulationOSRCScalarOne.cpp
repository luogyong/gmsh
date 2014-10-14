#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCScalarOne.h"
#include "FormulationOSRCHelper.h"

using namespace std;

FormulationOSRCScalarOne::FormulationOSRCScalarOne(void){
}

FormulationOSRCScalarOne::
FormulationOSRCScalarOne(const GroupOfElement& domain,
                         const FunctionSpace& field,
                         double k,
                         int NPade,
                         const TermFieldField<double>& localLHS,
                         const TermProjectionField<Complex>& localRHS){
  // Save Data //
  this->k        = k;
  this->ffield   = &field;
  this->ddomain  = &domain;
  this->localLHS = &localLHS;
  this->localRHS = &localRHS;

  // Pade C0
  C0 = FormulationOSRCHelper::padeC0(NPade, M_PI / 4.);
}

FormulationOSRCScalarOne::~FormulationOSRCScalarOne(void){
}

Complex FormulationOSRCScalarOne::weak(size_t dofI, size_t dofJ,
                                       size_t elementId) const{

  return Complex(0, -k) * C0 * localLHS->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCScalarOne::rhs(size_t equationI, size_t elementId) const{
  return localRHS->getTerm(equationI, 0, elementId);
}

const FunctionSpace& FormulationOSRCScalarOne::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCScalarOne::test(void) const{
  return *ffield;
}

const GroupOfElement& FormulationOSRCScalarOne::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCScalarOne::isBlock(void) const{
  return true;
}

void FormulationOSRCScalarOne::update(TermProjectionField<Complex>& localRHS){
  this->localRHS = &localRHS;
}
