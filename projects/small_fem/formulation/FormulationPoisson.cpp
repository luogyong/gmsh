#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationPoisson.h"

using namespace std;

// Poisson //
FormulationPoisson::FormulationPoisson(const GroupOfElement& domain,
                                       const FunctionSpaceScalar& fs,
                                       double (*f)(fullVector<double>& xyz)){
  // Save Domain //
  goe = &domain;

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationPoisson needs a uniform mesh");

  // Save FunctionSpace & Get Basis //
  const Basis& basis = fs.getBasis(eType);
  const size_t order = basis.getOrder();
  fspace             = &fs;

  // Source Term //
  fSource = f;

  // Gaussian Quadrature //
  Quadrature gaussGradGrad(eType, order - 1, 2);
  Quadrature gaussFF(eType, order, 2);

  const fullMatrix<double>& gCL = gaussGradGrad.getPoints();
  const fullVector<double>& gWL = gaussGradGrad.getWeights();

  const fullMatrix<double>& gCR = gaussFF.getPoints();
  const fullVector<double>& gWR = gaussFF.getWeights();

  // Local Terms //
  basis.preEvaluateDerivatives(gCL);
  basis.preEvaluateFunctions(gCR);

  GroupOfJacobian jacL(domain, gCL, "invert");
  GroupOfJacobian jacR(domain, gCR, "jacobian");

  localTermsL = new TermGradGrad(jacL, basis, gWL);
  localTermsR = new TermProjectionField<double>(jacR, basis, gWR, gCR, fSource);
}

FormulationPoisson::~FormulationPoisson(void){
  delete localTermsL;
  delete localTermsR;
}

double FormulationPoisson::weak(size_t dofI, size_t dofJ,
                                size_t elementId) const{

  return localTermsL->getTerm(dofI, dofJ, elementId);
}

double FormulationPoisson::rhs(size_t equationI,
                               size_t elementId) const{

  return localTermsR->getTerm(0, equationI, elementId);
}

const FunctionSpace& FormulationPoisson::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationPoisson::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationPoisson::domain(void) const{
  return *goe;
}
