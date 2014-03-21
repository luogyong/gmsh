#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
#include "FormulationOSRCFour.h"

using namespace std;

FormulationOSRCFour::FormulationOSRCFour(const GroupOfElement& domain,
                                         const FunctionSpaceScalar& fField,
                                         const FunctionSpaceScalar& fTest){

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOSRCFour needs a uniform mesh");

  // FunctionSpace (test and field) and Domain //
  ffField = &fField;
  ffTest  = &fTest;
  ddomain = &domain;

  // Basis (for test functions) //
  const Basis& basis = fTest.getBasis(eType);

  // Quadrature //
  Quadrature gauss(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Basis pre evaluation //
  basis.preEvaluateFunctions(gC);

  // Jacobians //
  GroupOfJacobian jac(domain, gC, "jacobian");

  // Local Term //
  local = new TermFieldField<double>(jac, basis, gauss);
}

FormulationOSRCFour::~FormulationOSRCFour(void){
  delete local;
}

Complex FormulationOSRCFour::weak(size_t dofI, size_t dofJ,
                                  size_t elementId) const{
  return
    Complex(-1, 0) * local->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOSRCFour::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationOSRCFour::field(void) const{
  return *ffField;
}

const FunctionSpace& FormulationOSRCFour::test(void) const{
  return *ffTest;
}

const GroupOfElement& FormulationOSRCFour::domain(void) const{
  return *ddomain;
}
