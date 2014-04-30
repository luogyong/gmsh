#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "Exception.h"
#include "FormulationUpdateOO2.h"

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

FormulationUpdateOO2::FormulationUpdateOO2(DDMContext& context){
  // Check if OO2 DDMContext //
  if(context.typeDDM != DDMContext::typeOO2)
    throw Exception("FormulationUpdateOO2 needs a OO2 DDMContext");

  // Get Domain and FunctionSpace from DDMContext //
  fspace  = &context.getFunctionSpace();
  ddomain = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateOO2 needs a uniform mesh");

  // a & b //
  this->a = context.OO2_A;
  this->b = context.OO2_B;

  // Get Basis //
  const Basis& basis = fspace->getBasis(eType);

  // Gaussian Quadrature (Field - Field) //
  Quadrature gaussFF(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gCFF = gaussFF.getPoints();

  // Gaussian Quadrature (Grad - Grad) //
  Quadrature gaussGG(eType, basis.getOrder() - 1, 2);
  const fullMatrix<double>& gCGG = gaussGG.getPoints();

  // Pre-evalution //
  basis.preEvaluateFunctions(gCFF);
  basis.preEvaluateDerivatives(gCGG);

  // Jacobians //
  GroupOfJacobian jacFF(*ddomain, gCFF, "jacobian");
  GroupOfJacobian jacGG(*ddomain, gCGG, "invert");

  // Get Volume Solution //
  map<Dof, Complex> sol;
  initMap(*fspace, *ddomain, sol);
  context.system->getSolution(sol, 0);

 // Get DDM Dofs from DDMContext //
  const map<Dof, Complex>& ddm = context.getDDMDofs();

  // Local Terms //
  lGout = new TermFieldField<double>(jacFF, basis, gaussFF);
  lGin  = new TermProjectionField<Complex>(jacFF, basis, gaussFF, *fspace, ddm);
  lU    = new TermProjectionField<Complex>(jacFF, basis, gaussFF, *fspace, sol);
  lDU   = new  TermProjectionGrad<Complex>(jacGG, basis, gaussGG, *fspace, sol);
}

FormulationUpdateOO2::~FormulationUpdateOO2(void){
  delete lGout;
  delete lGin;
  delete lU;
  delete lDU;
}

Complex FormulationUpdateOO2::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return Complex(lGout->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationUpdateOO2::rhs(size_t equationI, size_t elementId) const{
  return
    Complex(-1, 0)     * lGin->getTerm(0, equationI, elementId) +
    Complex(+2, 0) * a *   lU->getTerm(0, equationI, elementId) +
    Complex(-2, 0) * b *  lDU->getTerm(0, equationI, elementId);
}


const FunctionSpace& FormulationUpdateOO2::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationUpdateOO2::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationUpdateOO2::domain(void) const{
  return *ddomain;
}

bool FormulationUpdateOO2::isBlock(void) const{
  return true;
}
