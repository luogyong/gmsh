#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCThree.h"

using namespace std;

FormulationOSRCThree::FormulationOSRCThree(void){
}

FormulationOSRCThree::
FormulationOSRCThree(const GroupOfElement& domain,
                     const FunctionSpaceScalar& auxiliary,
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
  Bj = FormulationOSRC::padeBj(jPade, NPade, M_PI / 4.);
}

FormulationOSRCThree::~FormulationOSRCThree(void){
}

Complex FormulationOSRCThree::weak(size_t dofI, size_t dofJ,
                                   size_t elementId) const{
  return
     localFF->getTerm(dofI, dofJ, elementId)
    -localGG->getTerm(dofI, dofJ, elementId) * Bj / (keps * keps);
}

Complex FormulationOSRCThree::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCThree::field(void) const{
  return *faux;
}

const FunctionSpace& FormulationOSRCThree::test(void) const{
  return *faux;
}

const GroupOfElement& FormulationOSRCThree::domain(void) const{
  return *ddomain;
}
