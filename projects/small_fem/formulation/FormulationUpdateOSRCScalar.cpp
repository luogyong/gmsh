#include "Exception.h"

#include "FormulationHelper.h"
#include "FormulationOSRCHelper.h"

#include "FormulationOSRCScalar.h"
#include "FormulationUpdateOSRCScalar.h"

using namespace std;

FormulationUpdateOSRCScalar::
FormulationUpdateOSRCScalar(DDMContextOSRCScalar& context){
  // Save DDMContext //
  this->context = &context;

  // Get Domain and FunctionSpace from DDMContext //
  fspaceG   = &context.getFunctionSpaceG();   // DDM field: testing and unknwon
  fspace    = &context.getFunctionSpace();    // e field
  fspaceAux = &context.getAuxFunctionSpace(); // Auxiliary fields
  ddomain   = &context.getDomain();

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = ddomain->isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateOSRCScalar needs a uniform mesh");

  // Wavenumber //
  this->k = context.getWavenumber();

  // Pade //
  double theta = context.getRotation();
  NPade        = context.getNPade();
  A.resize(NPade);
  B.resize(NPade);

  C0 = FormulationOSRCHelper::padeC0(NPade, theta);

  for(int j = 0; j < NPade; j++)
    A[j] = FormulationOSRCHelper::padeA(j + 1, NPade, theta);
  for(int j = 0; j < NPade; j++)
    B[j] = FormulationOSRCHelper::padeB(j + 1, NPade, theta);

  // Basis //
  basis = &fspaceG->getBasis(eType);

  // Gaussian Quadrature //
  gauss = new Quadrature(eType, basis->getOrder(), 2);

  // Pre-evalution //
  const fullMatrix<double>& gC = gauss->getPoints();
  basis->preEvaluateFunctions(gC);

  // Jacobian //
  jac = new GroupOfJacobian(*ddomain, gC, "jacobian");

  // Init Volume & Auxiliary Solution //
  solPhi.resize(NPade);
  FormulationHelper::initDofMap(*fspace,    *ddomain, solU);
  FormulationHelper::initDofMap(*fspaceAux, *ddomain, solPhi);

  // Local Terms //
  lGout = new TermFieldField<double>(*jac, *basis, *gauss);

  // NB: lGin, lU & lPhi will be computed at update, when vol Sys is available
  lGin = NULL;
  lU   = NULL;
  lPhi = new TermProjectionField<Complex>*[NPade];

  for(int i = 0; i < NPade; i++)
    lPhi[i] = NULL;
}

FormulationUpdateOSRCScalar::~FormulationUpdateOSRCScalar(void){
  if(lGin)
    delete lGin;

  if(lU)
    delete lU;

  for(int i = 0; i < NPade; i++)
    if(lPhi[i])
      delete lPhi[i];

  delete[] lPhi;
  delete   lGout;
  delete   jac;
  delete   gauss;
}

Complex FormulationUpdateOSRCScalar::weak(size_t dofI, size_t dofJ,
                                          size_t elementId) const{
  return
    Complex(lGout->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationUpdateOSRCScalar::rhs(size_t equationI,
                                         size_t elementId) const{
  // Pade Sum //
  Complex sum = Complex(0, 0);

  for(int i = 0; i < NPade; i++)
    sum += A[i] / B[i] * (     lU->getTerm(equationI, 0, elementId) -
                          lPhi[i]->getTerm(equationI, 0, elementId));
  // FEM Term //
  return
    Complex(-1,  0    ) *      lGin->getTerm(equationI, 0, elementId) +
    Complex( 0, -2 * k) * C0  *  lU->getTerm(equationI, 0, elementId) +
    Complex( 0, -2 * k) * sum;
}

const FunctionSpace& FormulationUpdateOSRCScalar::field(void) const{
  return *fspaceG;
}

const FunctionSpace& FormulationUpdateOSRCScalar::test(void) const{
  return *fspaceG;
}

const GroupOfElement& FormulationUpdateOSRCScalar::domain(void) const{
  return *ddomain;
}

bool FormulationUpdateOSRCScalar::isBlock(void) const{
  return true;
}

void FormulationUpdateOSRCScalar::update(void){
  // Delete RHS (lGin, lU and lPhi)
  if(lGin)
    delete lGin;
  if(lU)
    delete lU;

  for(int i = 0; i < NPade; i++)
    if(lPhi[i])
      delete lPhi[i];

  // Get DDM Dofs, Volume and auxiliary solutions (at border) from DDMContext //
  const map<Dof, Complex>& ddm = context->getDDMDofs(); // ddm
  context->getSystem().getSolution(solU, 0);            // solU (field)

  for(int i = 0; i < NPade; i++)
    context->getSystem().getSolution(solPhi[i], 0);     // solPhi[] (aux)

  // Pre-evalution
  const fullMatrix<double>& gC = gauss->getPoints();
  basis->preEvaluateFunctions(gC);

  // New RHS
  lGin = new TermProjectionField<Complex>(*jac, *basis, *gauss, *fspaceG, ddm);
  lU   = new TermProjectionField<Complex>(*jac, *basis, *gauss, *fspace, solU);

  for(int i = 0; i < NPade; i++)
    lPhi[i] = new TermProjectionField<Complex>(*jac, *basis, *gauss,
                                               *(*fspaceAux)[i], solPhi[i]);
}
