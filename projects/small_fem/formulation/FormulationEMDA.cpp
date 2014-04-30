#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationEMDA.h"

using namespace std;

FormulationEMDA::FormulationEMDA(DDMContext& context){
  // Check if EMDA DDMContext //
  if(context.getType() != DDMContext::typeEMDA)
    throw Exception("FormulationEMDA needs a EMDA DDMContext");

  // Get Domain and FunctionSpace from DDMContext //
  fspace  = &context.getFunctionSpace();
  ddomain = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationEMDA needs a uniform mesh");

  // Wavenumber & Chi //
  this->k   = context.k;
  this->chi = context.EMDA_Chi;

  // Basis //
  const Basis& basis = fspace->getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Pre-evalution //
  basis.preEvaluateFunctions(gC);

  // Jacobian //
  GroupOfJacobian jac(*ddomain, gC, "jacobian");

  // Get DDM Dofs from DDMContext //
  const map<Dof, Complex>& ddm = context.getDDMDofs();

  // Local Terms //
  localLHS = new TermFieldField<double>(jac, basis, gauss);
  localRHS = new TermProjectionField<Complex>(jac, basis, gauss, *fspace, ddm);
}

FormulationEMDA::~FormulationEMDA(void){
  delete localLHS;
  delete localRHS;
}

Complex FormulationEMDA::weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return Complex(chi, -k) * localLHS->getTerm(dofI, dofJ, elementId);
}

Complex FormulationEMDA::rhs(size_t equationI, size_t elementId) const{
  return localRHS->getTerm(0, equationI, elementId);
}

const FunctionSpace& FormulationEMDA::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationEMDA::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationEMDA::domain(void) const{
  return *ddomain;
}

bool FormulationEMDA::isBlock(void) const{
  return true;
}
