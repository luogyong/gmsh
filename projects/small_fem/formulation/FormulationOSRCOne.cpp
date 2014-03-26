#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCOne.h"

using namespace std;

FormulationOSRCOne::FormulationOSRCOne(void){
}

FormulationOSRCOne::
FormulationOSRCOne(const GroupOfElement& domain,
                   const FunctionSpaceScalar& field,
                   double k,
                   const TermFieldField<double>& localLHS,
                   const TermProjectionField<Complex>& localRHS){
  // Save Data //
  this->k        = k;
  this->ffield   = &field;
  this->ddomain  = &domain;
  this->localLHS = &localLHS;
  this->localRHS = &localRHS;

  // Pade C0
  C0 = FormulationOSRC::padeC0(1, M_PI / 4.);
}

FormulationOSRCOne::~FormulationOSRCOne(void){
}

Complex FormulationOSRCOne::weak(size_t dofI, size_t dofJ,
                                 size_t elementId) const{

  return Complex(0, -k) * C0 * localLHS->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCOne::rhs(size_t equationI, size_t elementId) const{
  return localRHS->getTerm(0, equationI, elementId);
}

const FunctionSpace& FormulationOSRCOne::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCOne::test(void) const{
  return *ffield;
}

const GroupOfElement& FormulationOSRCOne::domain(void) const{
  return *ddomain;
}
