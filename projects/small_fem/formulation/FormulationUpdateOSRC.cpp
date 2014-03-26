#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"

#include "FormulationOSRC.h"
#include "FormulationUpdateOSRC.h"

using namespace std;

FormulationUpdateOSRC::FormulationUpdateOSRC(const GroupOfElement& domain,
                                             const FunctionSpaceScalar& fspace,
                                             double k,
                                             const map<Dof, Complex>& solU,
                                             const map<Dof, Complex>& solPhi,
                                             const map<Dof, Complex>& oldG){

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateOSRC needs a uniform mesh");

  // Wavenumber //
  this->k = k;

  // FunctionSpace and Domain //
  ffspace = &fspace;
  ddomain = &domain;

  // Pade //
  C0 = FormulationOSRC::padeC0(1, M_PI / 4.);
  A1 = FormulationOSRC::padeAj(1, 1, M_PI / 4.);
  B1 = FormulationOSRC::padeBj(1, 1, M_PI / 4.);

  // Basis //
  const Basis& basis = fspace.getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Jacobian //
  GroupOfJacobian jac(domain, gC, "jacobian");

  // Difference between solU and solPhi //
  if(solU.size() != solPhi.size())
    throw
      Exception("FormulationUpdateOSRC: solU and solPhi must have same size");

  map<Dof, Complex> UPhi = solU;

  map<Dof, Complex>::iterator       end    = UPhi.end();
  map<Dof, Complex>::iterator       itUPhi = UPhi.begin();
  map<Dof, Complex>::const_iterator itPhi  = solPhi.begin();

  for(; itUPhi != end; itUPhi++, itPhi++)
    itUPhi->second = itUPhi->second - itPhi->second;

  // Local Terms //
  lGout = new TermFieldField<double>(jac, basis, gauss);
  lGin  = new TermProjectionField<Complex>(jac, basis, gauss, fspace, oldG);
  lC0   = new TermProjectionField<Complex>(jac, basis, gauss, fspace, solU);
  lAB   = new TermProjectionField<Complex>(jac, basis, gauss, fspace, UPhi);
}

FormulationUpdateOSRC::~FormulationUpdateOSRC(void){
  delete lGout;
  delete lGin;
  delete lC0;
  delete lAB;
}

Complex FormulationUpdateOSRC::weak(size_t dofI, size_t dofJ,
                                    size_t elementId) const{
  return
    Complex(lGout->getTerm(dofI, dofJ, elementId), 0);
}

Complex FormulationUpdateOSRC::rhs(size_t equationI, size_t elementId) const{
  return
    Complex(-1,  0    ) *           lGin->getTerm(0, equationI, elementId) +
    Complex( 0, -2 * k) * C0 *       lC0->getTerm(0, equationI, elementId) +
    Complex( 0, -2 * k) * A1 / B1 *  lAB->getTerm(0, equationI, elementId);
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
