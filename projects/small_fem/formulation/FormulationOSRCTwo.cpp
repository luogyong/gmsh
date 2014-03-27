#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCTwo.h"

using namespace std;

FormulationOSRCTwo::FormulationOSRCTwo(void){
}

FormulationOSRCTwo::FormulationOSRCTwo(const GroupOfElement& domain,
                                       const FunctionSpaceScalar& auxiliary,
                                       const FunctionSpaceScalar& field,
                                       double  k,
                                       Complex keps,
                                       int NPade,
                                       int jPade,
                                       const TermGradGrad<double>& localTerm){
  // Save Data //
  this->k         = k;
  this->keps      = keps;
  this->ffield    = &field;
  this->faux      = &auxiliary;
  this->ddomain   = &domain;
  this->localTerm = &localTerm;

  // Pade Aj //
  Aj = FormulationOSRC::padeAj(jPade, NPade, M_PI / 4.);
}

FormulationOSRCTwo::~FormulationOSRCTwo(void){
}

Complex FormulationOSRCTwo::weak(size_t dofI, size_t dofJ,
                                 size_t elementId) const{
  return
    Complex(0, k) * Aj / (keps * keps) *
    localTerm->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCTwo::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCTwo::field(void) const{
  return *faux;
}

const FunctionSpace& FormulationOSRCTwo::test(void) const{
  return *ffield;
}

const GroupOfElement& FormulationOSRCTwo::domain(void) const{
  return *ddomain;
}
