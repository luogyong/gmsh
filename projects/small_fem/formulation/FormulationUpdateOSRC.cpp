#include "Exception.h"
#include "Quadrature.h"
#include "GroupOfJacobian.h"

#include "FormulationOSRC.h"
#include "FormulationUpdateOSRC.h"

using namespace std;

FormulationUpdateOSRC::
FormulationUpdateOSRC(const GroupOfElement& domain,
                      const FunctionSpaceScalar& fspace,
                      double k,
                      int NPade,
                      const map<Dof, Complex>& solU,
                      const vector<map<Dof, Complex> >& solPhi,
                      const map<Dof, Complex>& oldG){

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateOSRC needs a uniform mesh");

  // Sanity //
  if(solPhi.size() != (size_t)(NPade))
    throw Exception
      ("FormulationUpdateOSRC: size of solPhi should be equal to NPade");

  for(int j = 0; j < NPade; j++)
    if(solPhi[j].size() != solU.size())
      throw Exception
        ("FormulationUpdateOSRC: solPhi[j] and solU should have the same size");

  if(oldG.size() != solU.size())
    throw Exception
      ("FormulationUpdateOSRC: oldG and solU should have the same size");

  // Wavenumber //
  this->k = k;

  // FunctionSpace and Domain //
  ffspace = &fspace;
  ddomain = &domain;

  // Pade //
  A.resize(NPade);
  B.resize(NPade);

  C0 = FormulationOSRC::padeC0(NPade, M_PI / 4.);

  for(int j = 0; j < NPade; j++)
    A[j] = FormulationOSRC::padeAj(j + 1, NPade, M_PI / 4.);
  for(int j = 0; j < NPade; j++)
    B[j] = FormulationOSRC::padeBj(j + 1, NPade, M_PI / 4.);

  // Basis //
  const Basis& basis = fspace.getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Jacobian //
  GroupOfJacobian jac(domain, gC, "jacobian");

  // UPhi[d] = sum_j (solU[d] - solPhi[j][d]) * A[j] / B[j] //

  // Init UPhi
  map<Dof, Complex> UPhi = solU;

  map<Dof, Complex>:: iterator UPhiEnd = UPhi.end();
  map<Dof, Complex>:: iterator UPhiIt;

  for(UPhiIt = UPhi.begin(); UPhiIt != UPhiEnd; UPhiIt++)
    UPhiIt->second = Complex(0, 0);

  // Iterator on solU and solPhi[j]
  map<Dof, Complex>::const_iterator solUIt;
  map<Dof, Complex>::const_iterator solPhiJIt;

  // Loop on j (Pade terms) and Degrees of freedom (iterators)
  for(int j = 0; j < NPade; j++)
    for(UPhiIt = UPhi.begin(),
          solUIt = solU.begin(), solPhiJIt = solPhi[j].begin();
        UPhiIt != UPhiEnd; UPhiIt++, solUIt++, solPhiJIt++)

      UPhiIt->second += (solUIt->second - solPhiJIt->second) * A[j] / B[j];

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
    Complex(-1,  0    ) *      lGin->getTerm(0, equationI, elementId) +
    Complex( 0, -2 * k) * C0 * lC0->getTerm(0, equationI, elementId)  +
    Complex( 0, -2 * k) *      lAB->getTerm(0, equationI, elementId);
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
