#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationLagrangeTwo.h"

using namespace std;

FormulationLagrangeTwo::FormulationLagrangeTwo(void){
}

FormulationLagrangeTwo::
FormulationLagrangeTwo(const GroupOfElement& domain,
                       const FunctionSpaceScalar& field,
                       const FunctionSpaceScalar& test,
                       const TermFieldField<double>& localTerm){
  // Save Data //
  this->ddomain   = &domain;
  this->ffield    = &field;
  this->ttest     = &test;
  this->localTerm = &localTerm;
}

FormulationLagrangeTwo::~FormulationLagrangeTwo(void){
}

Complex FormulationLagrangeTwo::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  return Complex(localTerm->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationLagrangeTwo::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationLagrangeTwo::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationLagrangeTwo::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationLagrangeTwo::domain(void) const{
  return *ddomain;
}
