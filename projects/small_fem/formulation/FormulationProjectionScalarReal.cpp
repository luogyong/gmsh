#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "GroupOfElement.h"
#include "Quadrature.h"

#include "FormulationProjectionScalar.h"

using namespace std;

template<>
FormulationProjectionScalar<double>::
FormulationProjectionScalar(const GroupOfElement& goe,
                            const FunctionSpaceScalar& fs,
                            double (*f)(fullVector<double>& xyz)){
  // Save Domain //
  this->goe = &goe;

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = goe.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw ("FormulationProjectionScalar<real> needs a uniform mesh");

  // Save fspace //
  fspace = &fs;
  basis  = &fs.getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis->getOrder(), 2);

  const fullMatrix<double>& gC = gauss.getPoints();
  const fullVector<double>& gW = gauss.getWeights();

  // Local Terms //
  basis->preEvaluateFunctions(gC);
  GroupOfJacobian jac(goe, gC, "jacobian");

  localTerms1 = new TermFieldField(jac, *basis, gW);
  localTerms2 = new TermProjectionField(jac, *basis, gW, gC, f);
}

template<>
FormulationProjectionScalar<double>::~FormulationProjectionScalar(void){
  delete localTerms2;
  delete localTerms1;
}

template<>
double FormulationProjectionScalar<double>::weak(size_t dofI, size_t dofJ,
                                                 size_t elementId) const{

  return localTerms1->getTerm(dofI, dofJ, elementId);
}

template<>
double FormulationProjectionScalar<double>::rhs(size_t equationI,
                                                size_t elementId) const{

  return localTerms2->getTerm(0, equationI, elementId);
}

template<>
bool FormulationProjectionScalar<double>::isGeneral(void) const{
  return false;
}

template<>
double FormulationProjectionScalar<double>::weakB(size_t dofI, size_t dofJ,
                                                  size_t elementId) const{
  return 0;
}

template<>
const FunctionSpace& FormulationProjectionScalar<double>::fs(void) const{
  return *fspace;
}

template<>
const GroupOfElement& FormulationProjectionScalar<double>::domain(void) const{
  return *goe;
}
