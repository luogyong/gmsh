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

  // Init and type //
  int eType = initCommon(domain, fs);

  // Basis
  const Basis& basis = fs.getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Local Terms (0form) //
  basis.preEvaluateFunctions(gC);
  GroupOfJacobian jac(domain, gC, "jacobian");

  localTerms1 = new TermFieldField<scalar>(jac, basis, gauss);
  localTerms2 = new TermProjectionField<scalar>(jac, basis, gauss, f);
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

  // Init and type //
  int eType = initCommon(domain, fs);

  // Basis
  const Basis& basis = fs.getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Local Terms (1form) //
  basis.preEvaluateFunctions(gC);
  GroupOfJacobian jac(domain, gC, "invert");

  localTerms1 = new TermGradGrad<scalar>(jac, basis, gauss);
  localTerms2 = new TermProjectionGrad<scalar>(jac, basis, gauss, g);
}

template<typename scalar>
int FormulationProjection<scalar>::
initCommon(const GroupOfElement& domain,
           const FunctionSpace& fs){

  // Check GroupOfElement Stats: Uniform Mesh //
  std::pair<bool, size_t> uniform = domain.isUniform();
  size_t                    eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationProjection needs a uniform mesh");

  // Save Domain //
  ddomain = &domain;

  // Save fspace //
  fspace = &fs;

  // Return Type //
  return eType;
}

template<typename scalar>
FormulationProjection<scalar>::~FormulationProjection(void){
  delete localTerms2;
  delete localTerms1;
}

template<typename scalar>
scalar FormulationProjection<scalar>::weak(size_t dofI, size_t dofJ,
                                           size_t elementId) const{

  return localTerms1->getTerm(dofI, dofJ, elementId);
}

template<typename scalar>
scalar FormulationProjection<scalar>::rhs(size_t equationI,
                                          size_t elementId) const{

  return localTerms2->getTerm(equationI, 0, elementId);
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

template<typename scalar>
bool FormulationProjection<scalar>::isBlock(void) const{
  return true;
}
