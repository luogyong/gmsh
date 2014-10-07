#include "Exception.h"
#include "FormulationOSRC.h"
#include "FormulationHelper.h"
#include "FormulationUpdateOSRC.h"

using namespace std;

FormulationUpdateOSRC::FormulationUpdateOSRC(DDMContextOSRC& context){
  // Save DDMContext //
  this->context = &context;

  // Get Domain and FunctionSpace from DDMContext //
  ffspace = &context.getFunctionSpace();
  ddomain = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateOSRC needs a uniform mesh");

  // Wavenumber //
  this->k = context.getWavenumber();

  // Pade //
  NPade = context.getNPade();
  A.resize(NPade);
  B.resize(NPade);

  C0 = FormulationOSRC::padeC0(NPade, M_PI / 4.);

  for(int j = 0; j < NPade; j++)
    A[j] = FormulationOSRC::padeAj(j + 1, NPade, M_PI / 4.);
  for(int j = 0; j < NPade; j++)
    B[j] = FormulationOSRC::padeBj(j + 1, NPade, M_PI / 4.);

  // Basis //
  basis = &ffspace->getBasis(eType);

  // Gaussian Quadrature //
  gauss = new Quadrature(eType, basis->getOrder(), 2);

  // Pre-evalution //
  const fullMatrix<double>& gC = gauss->getPoints();
  basis->preEvaluateFunctions(gC);

  // Jacobian //
  jac = new GroupOfJacobian(*ddomain, gC, "jacobian");

  // Init Volume & Auxiliary Solution //
  solPhi.resize(NPade);
  FormulationHelper::initDofMap(*ffspace, *ddomain, solU);
  FormulationHelper::initDofMap(context.getAuxFunctionSpace(),
                                *ddomain, solPhi);

  // Init UPhi //
  UPhi = solU;
  resetUPhi();

  // Local Terms //
  lGout = new TermFieldField<double>(*jac, *basis, *gauss);

  // NB: lGin, lC0 and lAB will be computed at update, when vol Sys is available
  lGin  = NULL;
  lC0   = NULL;
  lAB   = NULL;
}

FormulationUpdateOSRC::~FormulationUpdateOSRC(void){
  delete lGout;
  if(lGin)
    delete lGin;
  if(lC0)
    delete lC0;
  if(lAB)
    delete lAB;

  delete jac;
  delete gauss;
}

Complex FormulationUpdateOSRC::weak(size_t dofI, size_t dofJ,
                                    size_t elementId) const{
  return
    Complex(lGout->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationUpdateOSRC::rhs(size_t equationI, size_t elementId) const{
  return
    Complex(-1,  0    ) *     lGin->getTerm(equationI, 0, elementId) +
    Complex( 0, -2 * k) * C0 * lC0->getTerm(equationI, 0, elementId) +
    Complex( 0, -2 * k) *      lAB->getTerm(equationI, 0, elementId);
}

const FunctionSpace& FormulationUpdateOSRC::field(void) const{
  return *ffspace;
}

const FunctionSpace& FormulationUpdateOSRC::test(void) const{
  return *ffspace;
}

const GroupOfElement& FormulationUpdateOSRC::domain(void) const{
  return *ddomain;
}

bool FormulationUpdateOSRC::isBlock(void) const{
  return true;
}

void FormulationUpdateOSRC::update(void){
  // Delete RHS (lGin, lC0 & lAB)
  if(lGin)
    delete lGin;
  if(lC0)
    delete lC0;
  if(lAB)
    delete lAB;

  // Get DDM Dofs, Volume and auxiliary solutions (at border) from DDMContext //
  const map<Dof, Complex>& ddm = context->getDDMDofs(); // ddm
  context->getSystem().getSolution(solU, 0);            // solU (field)

  for(int i = 0; i < NPade; i++)
    context->getSystem().getSolution(solPhi[i], 0);     // solPhi[] (aux)

  // UPhi[d] = sum_j (solU[d] - solPhi[j][d]) * A[j] / B[j] //
  resetUPhi();
  getUPhi();

  // Pre-evalution
  const fullMatrix<double>& gC = gauss->getPoints();
  basis->preEvaluateFunctions(gC);

  // New RHS
  lGin = new TermProjectionField<Complex>(*jac, *basis, *gauss, *ffspace, ddm);
  lC0  = new TermProjectionField<Complex>(*jac, *basis, *gauss, *ffspace, solU);
  lAB  = new TermProjectionField<Complex>(*jac, *basis, *gauss, *ffspace, UPhi);
}

void FormulationUpdateOSRC::resetUPhi(void){
  map<Dof, Complex>:: iterator UPhiEnd = UPhi.end();
  map<Dof, Complex>:: iterator UPhiIt;

  for(UPhiIt = UPhi.begin(); UPhiIt != UPhiEnd; UPhiIt++)
    UPhiIt->second = Complex(0, 0);
}

void FormulationUpdateOSRC::getUPhi(void){
  // UPhi[d] = sum_j (solU[d] - solPhi[j][d]) * A[j] / B[j] //

  // Iterator on solU and solPhi[j]
  map<Dof, Complex>:: iterator UPhiEnd = UPhi.end();
  map<Dof, Complex>:: iterator UPhiIt;
  map<Dof, Complex>::const_iterator solUIt;
  map<Dof, Complex>::const_iterator solPhiJIt;

  // Loop on j (Pade terms) and Degrees of freedom (iterators)
  for(int j = 0; j < NPade; j++)
    for(UPhiIt = UPhi.begin(),
          solUIt = solU.begin(), solPhiJIt = solPhi[j].begin();
        UPhiIt != UPhiEnd; UPhiIt++, solUIt++, solPhiJIt++)

      UPhiIt->second += (solUIt->second - solPhiJIt->second) * A[j] / B[j];
}
