#include "FormulationPML.h"

#include "GroupOfJacobian.h"
#include "Quadrature.h"

using namespace std;

FormulationPML::
FormulationPML(const GroupOfElement& domain,
               const FunctionSpace& fs,
               double k,
               void   (*fS)(fullVector<double>& xyz, fullMatrix<Complex>& T),
               Complex (*fM)(fullVector<double>& xyz)){

  // Check domain stats: uniform mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationPML needs a uniform mesh");

  // Save Data //
  ddomain = &domain;
  ffs     = &fs;

  // Wavenumber squared //
  this->kSquare = k * k;

  // Get Basis //
  const Basis& basis = fs.getBasis(eType);
  const size_t order = basis.getOrder();

  // Gaussian Quadrature //
  Quadrature gaussStif(eType, order - 1, 2);
  Quadrature gaussMass(eType, order,     2);

  const fullMatrix<double>& gCS = gaussStif.getPoints();
  const fullMatrix<double>& gCM = gaussMass.getPoints();

  // Functions //
  basis.preEvaluateDerivatives(gCS);
  basis.preEvaluateFunctions(gCM);

  // Local Terms //
  GroupOfJacobian jacS(domain, gCS, "invert");
  GroupOfJacobian jacM(domain, gCM, "jacobian");

  stif = new TermGradGrad<Complex>  (jacS, basis, gaussStif, fS);
  mass = new TermFieldField<Complex>(jacM, basis, gaussMass, fM);
}

FormulationPML::~FormulationPML(void){
  delete stif;
  delete mass;
}

Complex FormulationPML::weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return
                stif->getTerm(dofI, dofJ, elementId)
    - kSquare * mass->getTerm(dofI, dofJ, elementId);
}

Complex FormulationPML::rhs(size_t equationI, size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationPML::field(void) const{
  return *ffs;
}

const FunctionSpace& FormulationPML::test(void) const{
  return *ffs;
}

const GroupOfElement& FormulationPML::domain(void) const{
  return *ddomain;
}
