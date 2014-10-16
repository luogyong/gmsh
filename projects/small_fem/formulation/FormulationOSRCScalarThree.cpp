#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCScalarThree.h"
#include "FormulationOSRCHelper.h"

using namespace std;

FormulationOSRCScalarThree::FormulationOSRCScalarThree(void){
}

FormulationOSRCScalarThree::
FormulationOSRCScalarThree(const GroupOfElement& domain,
                           const FunctionSpace& auxiliary,
                           Complex keps,
                           int NPade,
                           int jPade,
                           const TermFieldField<double>& localFF,
                           const TermGradGrad<double>& localGG){
  // Save Data //
  this->keps    = keps;
  this->faux    = &auxiliary;
  this->ddomain = &domain;
  this->localFF = &localFF;
  this->localGG = &localGG;

  // Pade Bj //
  Bj = FormulationOSRCHelper::padeB(jPade, NPade, M_PI / 4.);
}

FormulationOSRCScalarThree::~FormulationOSRCScalarThree(void){
}

Complex FormulationOSRCScalarThree::weak(size_t dofI, size_t dofJ,
                                         size_t elementId) const{
  return
     localFF->getTerm(dofI, dofJ, elementId)
    -localGG->getTerm(dofI, dofJ, elementId) * Bj / (keps * keps);
}

Complex FormulationOSRCScalarThree::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCScalarThree::field(void) const{
  return *faux;
}

const FunctionSpace& FormulationOSRCScalarThree::test(void) const{
  return *faux;
}

const GroupOfElement& FormulationOSRCScalarThree::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCScalarThree::isBlock(void) const{
  return true;
}
