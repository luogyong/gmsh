////////////////////////////////////////////////////////////////
// Templates Implementations for FormulationProjectionVector: //
// Inclusion compilation model                                //
//                                                            //
// Damn you gcc: we want 'export' !                           //
////////////////////////////////////////////////////////////////

#include "GroupOfJacobian.h"
#include "Quadrature.h"

template<typename scalar>
FormulationProjectionVector<scalar>::
FormulationProjectionVector(const GroupOfElement& domain,
                            const FunctionSpaceVector& fs,
                            fullVector<scalar> (*f)(fullVector<double>& xyz)){
  // Save Domain //
  ddomain = &domain;

  // Check GroupOfElement Stats: Uniform Mesh //
  std::pair<bool, size_t> uniform = domain.isUniform();
  size_t                    eType = uniform.second;

  if(!uniform.first)
    throw ("FormulationProjectionVector needs a uniform mesh");

  // Save fspace //
  fspace = &fs;
  const Basis& basis = fs.getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis.getOrder(), 2);

  const fullMatrix<double>& gC = gauss.getPoints();
  const fullVector<double>& gW = gauss.getWeights();

  // Local Terms //
  basis.preEvaluateFunctions(gC);
  GroupOfJacobian jac(domain, gC, "invert");

  localTerms1 = new TermGradGrad(jac, basis, gW);
  localTerms2 = new TermProjectionGrad<scalar>(jac, basis, gW, gC, f);
}

template<typename scalar>
FormulationProjectionVector<scalar>::~FormulationProjectionVector(void){
  delete localTerms1;
  delete localTerms2;
}

template<typename scalar>
const FunctionSpace& FormulationProjectionVector<scalar>::field(void) const{
  return *fspace;
}

template<typename scalar>
const FunctionSpace& FormulationProjectionVector<scalar>::test(void) const{
  return *fspace;
}

template<typename scalar>
const GroupOfElement& FormulationProjectionVector<scalar>::domain(void) const{
  return *ddomain;
}
