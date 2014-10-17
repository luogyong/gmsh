#include "Exception.h"
#include "FormulationOSRCVectorThree.h"

using namespace std;

FormulationOSRCVectorThree::FormulationOSRCVectorThree(void){
}

FormulationOSRCVectorThree::
FormulationOSRCVectorThree(const GroupOfElement& domain,
                           const FunctionSpace& field,
                           const FunctionSpace& test,
                           Complex R0,
                           const TermGradGrad<double>& localGG,
                           const TermProjectionGrad<Complex>& localRHS){
  // Save Data //
  this->oneOverR0 = Complex(1, 0) / R0;
  this->ffield    = &field;
  this->ttest     = &test;
  this->ddomain   = &domain;
  this->localGG   = &localGG;
  this->localRHS  = &localRHS;
}

FormulationOSRCVectorThree::~FormulationOSRCVectorThree(void){
}

Complex FormulationOSRCVectorThree::weak(size_t dofI, size_t dofJ,
                                         size_t elementId) const{

  return localGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCVectorThree::rhs(size_t equationI,
                                        size_t elementId) const{

  return oneOverR0 * localRHS->getTerm(equationI, 0, elementId);
}

const FunctionSpace& FormulationOSRCVectorThree::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorThree::test(void) const{
  return *ttest;
}

const GroupOfElement& FormulationOSRCVectorThree::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorThree::isBlock(void) const{
  return true;
}

void FormulationOSRCVectorThree::update(TermProjectionGrad<Complex>& localRHS){
  this->localRHS = &localRHS;
}
