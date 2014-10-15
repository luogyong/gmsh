#include "Exception.h"
#include "FormulationOSRCHelper.h"
#include "FormulationOSRCVectorOne.h"

using namespace std;

FormulationOSRCVectorOne::FormulationOSRCVectorOne(void){
}

FormulationOSRCVectorOne::
FormulationOSRCVectorOne(const GroupOfElement& domain,
                         const FunctionSpace& field,
                         double  k,
                         Complex C0 ,
                         const TermGradGrad<double>& localLHS,
                         const TermProjectionGrad<Complex>& localRHS){
  // Save Data //
  this->jOverK   = Complex(0, 1. / k);
  this ->C0      = C0;
  this->ffield   = &field;
  this->ddomain  = &domain;
  this->localLHS = &localLHS;
  this->localRHS = &localRHS;
}

FormulationOSRCVectorOne::~FormulationOSRCVectorOne(void){
}

Complex FormulationOSRCVectorOne::weak(size_t dofI, size_t dofJ,
                                       size_t elementId) const{

  return jOverK * C0 * localLHS->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCVectorOne::rhs(size_t equationI, size_t elementId) const{
  return localRHS->getTerm(equationI, 0, elementId);
}

const FunctionSpace& FormulationOSRCVectorOne::field(void) const{
  return *ffield;
}

const FunctionSpace& FormulationOSRCVectorOne::test(void) const{
  return *ffield;
}

const GroupOfElement& FormulationOSRCVectorOne::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCVectorOne::isBlock(void) const{
  return true;
}

void FormulationOSRCVectorOne::update(TermProjectionGrad<Complex>& localRHS){
  this->localRHS = &localRHS;
}
