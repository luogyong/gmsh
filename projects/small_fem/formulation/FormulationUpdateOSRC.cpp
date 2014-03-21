#include <cmath>

#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"
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
  C0 = padeC0(1, M_PI / 4.);
  A1 = padeAj(1, 1, M_PI / 4.);
  B1 = padeBj(1, 1, M_PI / 4.);

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

double FormulationUpdateOSRC::pade_aj(int j, int N){
  return 2. / (2. * N + 1.) * sqrt(sin((double)j * M_PI / (2. * N + 1.)));
}

double FormulationUpdateOSRC::pade_bj(int j, int N){
  return sqrt(cos((double)j * M_PI / (2. *N + 1.)));
}

Complex FormulationUpdateOSRC::padeC0(int N, double theta){
  Complex sum = Complex(1, 0);
  Complex one = Complex(1, 0);
  Complex z   = Complex(cos(-theta) - 1,  sin(-theta));

  for(int j = 1; j <= N; j++)
    sum += (z * pade_aj(j, N)) / (one + z * pade_bj(j, N));

  z = Complex(cos(theta / 2.), sin(theta / 2.));

  return sum * z;
}

Complex FormulationUpdateOSRC::padeAj(int j, int N, double theta){
  Complex one = Complex(1, 0);
  Complex res;
  Complex z;

  z   = Complex(cos(-theta / 2.), sin(-theta / 2.));
  res = z * pade_aj(j, N);

  z   = Complex(cos(-theta) - 1., sin(-theta));
  res = res / ((one + z * pade_bj(j, N)) * (one + z * pade_bj(j, N)));

  return res;
}

Complex FormulationUpdateOSRC::padeBj(int j, int N, double theta){
  Complex one = Complex(1, 0);
  Complex res;
  Complex z;

  z   = Complex(cos(-theta), sin(-theta));
  res = z * pade_bj(j, N);

  z   = Complex(cos(-theta) - 1., sin(-theta));
  res = res / (one + z * pade_bj(j, N));

  return res;
}
