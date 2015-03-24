#include "SmallFem.h"
#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"

#include "TermFieldField.h"
#include "TermGradGrad.h"

#include "FormulationImpedance.h"

using namespace std;


FormulationImpedance::FormulationImpedance(const GroupOfElement& domain,
                                           const FunctionSpace& fs,
                                           double k, Complex epsr, Complex mur){
  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationImpedance needs a uniform mesh");

  // Save Domain //
  goe = &domain;

  // Wavenumber, epsilon_r and mu_r //
  this->k    = k;
  this->epsr = epsr;
  this->mur  = mur;

  // Save FunctionSpace & Get Basis //
  const Basis& basis = fs.getBasis(eType);
  const size_t order = basis.getOrder();
  fspace             = &fs;

  // Gaussian Quadrature //
  Quadrature gauss(eType, order, 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Local Terms //
  basis.preEvaluateFunctions(gC);

  if(fs.isScalar()){
    GroupOfJacobian jac(domain, gC, "jacobian");
    localTerms = new TermFieldField<double>(jac, basis, gauss);
  }

  else{
    GroupOfJacobian jac(domain, gC, "invert");
    localTerms = new TermGradGrad<double>(jac, basis, gauss);
  }
}

FormulationImpedance::~FormulationImpedance(void){
  delete localTerms;
}

Complex FormulationImpedance::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return Complex(0, -1 * k) * std::sqrt(epsr * mur)
                            * localTerms->getTerm(dofI, dofJ, elementId);
}

Complex FormulationImpedance::rhs(size_t equationI, size_t elementId) const{
  return Complex(0, 0);
}

const FunctionSpace& FormulationImpedance::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationImpedance::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationImpedance::domain(void) const{
  return *goe;
}

bool FormulationImpedance::isBlock(void) const{
  return true;
}
