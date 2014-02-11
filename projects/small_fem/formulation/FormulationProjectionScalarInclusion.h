////////////////////////////////////////////////////////////////
// Templates Implementations for FormulationProjectionScalar: //
// Inclusion compilation model                                //
//                                                            //
// Damn you gcc: we want 'export' !                           //
////////////////////////////////////////////////////////////////

#include "GroupOfJacobian.h"
#include "Quadrature.h"

template<typename scalar>
FormulationProjectionScalar<scalar>::
FormulationProjectionScalar(const GroupOfElement& domain,
                            const FunctionSpaceScalar& fs,
                            scalar (*f)(fullVector<double>& xyz)){
  // Save Domain //
  ddomain = &domain;

  // Check GroupOfElement Stats: Uniform Mesh //
  std::pair<bool, size_t> uniform = domain.isUniform();
  size_t                    eType = uniform.second;

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
  GroupOfJacobian jac(domain, gC, "jacobian");

  localTerms1 = new TermFieldField(jac, *basis, gW);
  localTerms2 = new TermProjectionField<scalar>(jac, *basis, gW, gC, f);
}

template<typename scalar>
FormulationProjectionScalar<scalar>::~FormulationProjectionScalar(void){
  delete localTerms2;
  delete localTerms1;
}

template<typename scalar>
const FunctionSpace& FormulationProjectionScalar<scalar>::field(void) const{
  return *fspace;
}

template<typename scalar>
const FunctionSpace& FormulationProjectionScalar<scalar>::test(void) const{
  return *fspace;
}

template<typename scalar>
const GroupOfElement& FormulationProjectionScalar<scalar>::domain(void) const{
  return *ddomain;
}
