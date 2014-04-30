#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "Exception.h"
#include "FormulationUpdateEMDA.h"

using namespace std;

#include <set>
static
void initMap(const FunctionSpace& fs,
             const GroupOfElement& goe, map<Dof, Complex>& data){

  set<Dof> dSet;
  fs.getKeys(goe, dSet);

  set<Dof>::iterator it  = dSet.begin();
  set<Dof>::iterator end = dSet.end();

  for(; it != end; it++)
    data.insert(pair<Dof, Complex>(*it, 0));
}

FormulationUpdateEMDA::FormulationUpdateEMDA(DDMContext& context){
  // Check if EMDA DDMContext //
  if(context.typeDDM != DDMContext::typeEMDA)
    throw Exception("FormulationUpdateEMDA needs a EMDA DDMContext");

  // Get Domain and FunctionSpace from DDMContext //
  fspace  = &context.getFunctionSpace();
  ddomain = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateEMDA needs a uniform mesh");

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

  // Get Volume Solution //
  map<Dof, Complex> sol;
  initMap(*fspace, *ddomain, sol);
  context.system->getSolution(sol, 0);

  // Get DDM Dofs from DDMContext //
  const map<Dof, Complex>& ddm = context.getDDMDofs();

  // Local Terms //
  lGout = new TermFieldField<double>(jac, basis, gauss);
  lGin  = new TermProjectionField<Complex>(jac, basis, gauss, *fspace, ddm);
  lU    = new TermProjectionField<Complex>(jac, basis, gauss, *fspace, sol);
}

FormulationUpdateEMDA::~FormulationUpdateEMDA(void){
  delete lGout;
  delete lGin;
  delete lU;
}

Complex FormulationUpdateEMDA::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return Complex(lGout->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationUpdateEMDA::rhs(size_t equationI, size_t elementId) const{
  return
    Complex(-1      , 0     ) * lGin->getTerm(0, equationI, elementId) +
    Complex(+2 * chi, -2 * k) *   lU->getTerm(0, equationI, elementId);
}

const FunctionSpace& FormulationUpdateEMDA::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationUpdateEMDA::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationUpdateEMDA::domain(void) const{
  return *ddomain;
}

bool FormulationUpdateEMDA::isBlock(void) const{
  return true;
}
