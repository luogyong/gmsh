#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationLagrangeOne.h"

using namespace std;

FormulationLagrangeOne::FormulationLagrangeOne(void){
}

FormulationLagrangeOne::
FormulationLagrangeOne(const GroupOfElement& domain,
                       const FunctionSpaceScalar& field,
                       const FunctionSpaceScalar& test,
                       const TermFieldField<double>& localTerm,
                       const TermProjectionField<double>& projectionTerm){
  // Save Data //
  this->ddomain        = &domain;
  this->ffield         = &field;
  this->ttest          = &test;
  this->localTerm      = &localTerm;
  this->projectionTerm = &projectionTerm;
}

FormulationLagrangeOne::~FormulationLagrangeOne(void){
}

Complex FormulationLagrangeOne::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  return Complex(localTerm->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationLagrangeOne::rhs(size_t equationI, size_t elementId) const{
  return Complex(projectionTerm->getTerm(equationI, 0, elementId), 0);
}

const FunctionSpace& FormulationLagrangeOne::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationLagrangeOne::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationLagrangeOne::domain(void) const{
  return *ddomain;
}

bool FormulationLagrangeOne::isBlock(void) const{
  return true;
}
