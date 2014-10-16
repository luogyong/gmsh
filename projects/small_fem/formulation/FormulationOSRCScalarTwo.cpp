#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCScalarTwo.h"
#include "FormulationOSRCHelper.h"

using namespace std;

FormulationOSRCScalarTwo::FormulationOSRCScalarTwo(void){
}

FormulationOSRCScalarTwo::
FormulationOSRCScalarTwo(const GroupOfElement& domain,
                         const FunctionSpace& auxiliary,
                         const FunctionSpace& field,
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
  Aj = FormulationOSRCHelper::padeA(jPade, NPade, M_PI / 4.);
}

FormulationOSRCScalarTwo::~FormulationOSRCScalarTwo(void){
}

Complex FormulationOSRCScalarTwo::weak(size_t dofI, size_t dofJ,
                                       size_t elementId) const{
  return
    Complex(0, k) * Aj / (keps * keps) *
    localTerm->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCScalarTwo::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCScalarTwo::field(void) const{
  return *faux;
}

const FunctionSpace& FormulationOSRCScalarTwo::test(void) const{
  return *ffield;
}

const GroupOfElement& FormulationOSRCScalarTwo::domain(void) const{
  return *ddomain;
}

bool FormulationOSRCScalarTwo::isBlock(void) const{
  return true;
}
