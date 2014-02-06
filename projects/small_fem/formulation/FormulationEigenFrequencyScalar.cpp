#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationEigenFrequencyScalar.h"

using namespace std;

FormulationEigenFrequencyScalar::
FormulationEigenFrequencyScalar(const GroupOfElement& goe,
                                const FunctionSpaceScalar& fs){
  // Save Domain //
  this->goe = &goe;

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = goe.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationEigenFrequencyScalar needs a uniform mesh");

  // Save FunctionSpace & Get Basis //
  const Basis& basis = fs.getBasis(eType);
  const size_t order = basis.getOrder();
  fspace             = &fs;

  // Gaussian Quadrature //
  Quadrature gaussGradGrad(eType, order - 1, 2);
  Quadrature gaussFF(eType, order, 2);

  const fullMatrix<double>& gC1 = gaussGradGrad.getPoints();
  const fullVector<double>& gW1 = gaussGradGrad.getWeights();

  const fullMatrix<double>& gC2 = gaussFF.getPoints();
  const fullVector<double>& gW2 = gaussFF.getWeights();

  // Local Terms //
  basis.preEvaluateDerivatives(gC1);
  basis.preEvaluateFunctions(gC2);

  GroupOfJacobian jac1(goe, gC1, "invert");
  GroupOfJacobian jac2(goe, gC2, "jacobian");

  localTerms1 = new TermGradGrad(jac1, basis, gW1);
  localTerms2 = new TermFieldField(jac2, basis, gW2);
}

FormulationEigenFrequencyScalar::~FormulationEigenFrequencyScalar(void){
  delete localTerms1;
  delete localTerms2;
}

std::complex<double>
FormulationEigenFrequencyScalar::weak(size_t dofI, size_t dofJ,
                                      size_t elementId) const{

  return std::complex<double>(localTerms1->getTerm(dofI, dofJ, elementId), 0);
}


std::complex<double>
FormulationEigenFrequencyScalar::weakB(size_t dofI, size_t dofJ,
                                       size_t elementId) const{

  return std::complex<double>(localTerms2->getTerm(dofI, dofJ, elementId), 0);
}

std::complex<double>
FormulationEigenFrequencyScalar::rhs(size_t equationI, size_t elementId) const{
  return std::complex<double>(0, 0);
}

bool FormulationEigenFrequencyScalar::isGeneral(void) const{
  return true;
}

const FunctionSpace& FormulationEigenFrequencyScalar::fsField(void) const{
  return *fspace;
}

const FunctionSpace& FormulationEigenFrequencyScalar::fsTest(void) const{
  return *fspace;
}

const GroupOfElement& FormulationEigenFrequencyScalar::domain(void) const{
  return *goe;
}
