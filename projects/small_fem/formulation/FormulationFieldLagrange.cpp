#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationFieldLagrange.h"

using namespace std;

FormulationFieldLagrange::
FormulationFieldLagrange(const GroupOfElement& goe,
                         const FunctionSpaceScalar& fsField,
                         const FunctionSpaceScalar& fsTest,
                         double (*f)(fullVector<double>& xyz)){
  // Save Domain //
  this->goe = &goe;

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = goe.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationFieldLagrange needs a uniform mesh");

  // Save FunctionSpace & Get Basis //
  const Basis& basis = fsTest.getBasis(eType);
  const size_t order = basis.getOrder();
  fsF                = &fsField;
  fsT                = &fsTest;

  // Gaussian Quadrature //
  Quadrature gaussFF(eType, order, 2);
  const fullMatrix<double>& gC = gaussFF.getPoints();
  const fullVector<double>& gW = gaussFF.getWeights();

  // Local Terms //
  basis.preEvaluateFunctions(gC);

  GroupOfJacobian jac(goe, gC, "jacobian");

  localTerms      = new TermFieldField(jac, basis, gW);
  projectionTerms = new TermProjectionField(jac, basis, gW, gC, f);
}

FormulationFieldLagrange::~FormulationFieldLagrange(void){
  delete localTerms;
  delete projectionTerms;
}

complex<double> FormulationFieldLagrange::weak(size_t dofI, size_t dofJ,
                                               size_t elementId) const{
  return
    complex<double>(localTerms->getTerm(dofI, dofJ, elementId), 0);
}

complex<double> FormulationFieldLagrange::rhs(size_t equationI,
                                              size_t elementId) const{
  return complex<double>(projectionTerms->getTerm(0, equationI, elementId), 0);
}

bool FormulationFieldLagrange::isGeneral(void) const{
  return false;
}

complex<double> FormulationFieldLagrange::weakB(size_t dofI, size_t dofJ,
                                                size_t elementId) const{
  return complex<double>(0, 0);
}


const FunctionSpace& FormulationFieldLagrange::fsField(void) const{
  return *fsF;
}

const FunctionSpace& FormulationFieldLagrange::fsTest(void) const{
  return *fsT;
}

const GroupOfElement& FormulationFieldLagrange::domain(void) const{
  return *goe;
}
