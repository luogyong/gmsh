#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationEigenFrequencyVector.h"

using namespace std;

FormulationEigenFrequencyVector::
FormulationEigenFrequencyVector(GroupOfElement& goe,
                                const FunctionSpaceVector& fs){

  // Check GroupOfElement Stats: Uniform Mesh //
  const vector<size_t>& gType = goe.getTypeStats();
  const size_t nGType = gType.size();
  size_t eType = (size_t)(-1);

  for(size_t i = 0; i < nGType; i++)
    if((eType == (size_t)(-1)) && (gType[i] != 0))
      eType = i;
    else if((eType != (size_t)(-1)) && (gType[i] != 0))
      throw Exception("FormulationEigenFrequencyVector needs a uniform mesh");

  // Save FunctionSpace & Get Basis //
  const Basis& basis = fs.getBasis(eType);
  const size_t order = basis.getOrder();
  fspace             = &fs;

  // Gaussian Quadrature //
  Quadrature gaussCurlCurl(eType, order - 1, 2);
  Quadrature gaussFF(eType, order, 2);

  const fullMatrix<double>& gC1 = gaussCurlCurl.getPoints();
  const fullVector<double>& gW1 = gaussCurlCurl.getWeights();

  const fullMatrix<double>& gC2 = gaussFF.getPoints();
  const fullVector<double>& gW2 = gaussFF.getWeights();

  // Local Terms //
  basis.preEvaluateDerivatives(gC1);
  basis.preEvaluateFunctions(gC2);

  GroupOfJacobian jac1(goe, gC1, "jacobian");
  GroupOfJacobian jac2(goe, gC2, "invert");

  localTerms1 = new TermCurlCurl(jac1, basis, gW1);
  localTerms2 = new TermGradGrad(jac2, basis, gW2);
}

FormulationEigenFrequencyVector::~FormulationEigenFrequencyVector(void){
  delete localTerms1;
  delete localTerms2;
}

std::complex<double>
FormulationEigenFrequencyVector::weak(size_t dofI, size_t dofJ,
                                      size_t elementId) const{

  return std::complex<double>(localTerms1->getTerm(dofI, dofJ, elementId), 0);
}

std::complex<double>
FormulationEigenFrequencyVector::weakB(size_t dofI, size_t dofJ,
                                       size_t elementId) const{

  return std::complex<double>(localTerms2->getTerm(dofI, dofJ, elementId), 0);
}

std::complex<double>
FormulationEigenFrequencyVector::rhs(size_t dofI, size_t elementId) const{
  return std::complex<double>(0, 0);
}

bool FormulationEigenFrequencyVector::isGeneral(void) const{
  return true;
}

const FunctionSpace& FormulationEigenFrequencyVector::fs(void) const{
  return *fspace;
}
