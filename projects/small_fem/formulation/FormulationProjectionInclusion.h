//////////////////////////////////////////////////////////
// Templates Implementations for FormulationProjection: //
// Inclusion compilation model                          //
//                                                      //
// Damn you gcc: we want 'export' !                     //
//////////////////////////////////////////////////////////

#include "GroupOfJacobian.h"
#include "Quadrature.h"
#include "Exception.h"

#include "TermProjectionField.h"
#include "TermProjectionGrad.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"

template<typename scalar>
FormulationProjection<scalar>::
FormulationProjection(const GroupOfElement& domain,
                      const FunctionSpace& fs,
                      scalar (*f)(fullVector<double>& xyz)){

  // Check FunctionSpace: we need a 0 form //
  if(fs.getForm() != 0)
    throw Exception("FormulationProjection: %s %s",
                    "projection of a scalar function requires",
                    "a 0 form function space");
  // Init //
  fullMatrix<double> gC;
  fullVector<double> gW;
  const Basis& basis = initCommon(domain, fs, gC, gW);

  // Local Terms (0form) //
  basis.preEvaluateFunctions(gC);
  GroupOfJacobian jac(domain, gC, "jacobian");

  localTerms1 = new TermFieldField(jac, basis, gW);
  localTerms2 = new TermProjectionField<scalar>(jac, basis, gW, gC, f);
}

template<typename scalar>
FormulationProjection<scalar>::
FormulationProjection(const GroupOfElement& domain,
                      const FunctionSpace& fs,
                      fullVector<scalar> (*g)(fullVector<double>& xyz)){

  // Check FunctionSpace: we need a 1 form //
  if(fs.getForm() != 1)
    throw Exception("FormulationProjection: %s %s",
                    "projection of a vector function requires",
                    "a 1 form function space");
  // Init //
  fullMatrix<double> gC;
  fullVector<double> gW;
  const Basis& basis = initCommon(domain, fs, gC, gW);

  // Local Terms (1form) //
  basis.preEvaluateFunctions(gC);
  GroupOfJacobian jac(domain, gC, "invert");

  localTerms1 = new TermGradGrad(jac, basis, gW);
  localTerms2 = new TermProjectionGrad<scalar>(jac, basis, gW, gC, g);
}

template<typename scalar>
const Basis& FormulationProjection<scalar>::
initCommon(const GroupOfElement& domain,
           const FunctionSpace& fs,
           fullMatrix<double>& gC,
           fullVector<double>& gW){

  // Save Domain //
  ddomain = &domain;

  // Check GroupOfElement Stats: Uniform Mesh //
  std::pair<bool, size_t> uniform = domain.isUniform();
  size_t                    eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationProjection needs a uniform mesh");

  // Save fspace //
  fspace = &fs;
  const Basis& basis = fs.getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis.getOrder(), 2);

  gC = gauss.getPoints();
  gW = gauss.getWeights();

  // Return Basis //
  return basis;
}

template<typename scalar>
FormulationProjection<scalar>::~FormulationProjection(void){
  delete localTerms2;
  delete localTerms1;
}

template<typename scalar>
const FunctionSpace& FormulationProjection<scalar>::field(void) const{
  return *fspace;
}

template<typename scalar>
const FunctionSpace& FormulationProjection<scalar>::test(void) const{
  return *fspace;
}

template<typename scalar>
const GroupOfElement& FormulationProjection<scalar>::domain(void) const{
  return *ddomain;
}
