//////////////////////////////////////////////////////////
// Templates Implementations for FormulationSteadyWave: //
// Inclusion compilation model                          //
//                                                      //
// Damn you gcc: we want 'export' !                     //
//////////////////////////////////////////////////////////

#include "GroupOfJacobian.h"
#include "Quadrature.h"
#include "Exception.h"

#include "TermProjectionGrad.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"
#include "TermCurlCurl.h"

template<typename scalar>
FormulationSteadyWave<scalar>::
FormulationSteadyWave(const GroupOfElement& domain,
                      const FunctionSpace& fs,
                      double k){

  // Check domain stats: uniform mesh //
  std::pair<bool, size_t> uniform = domain.isUniform();
  size_t                    eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationSteadyWave needs a uniform mesh");

  // Save Data //
  ddomain = &domain;
  ffs     = &fs;

  // Wave Squared //
  kSquare = k * k;

  // Get Basis //
  const size_t form  = fs.getForm();
  const Basis& basis = fs.getBasis(eType);
  const size_t order = basis.getOrder();

  // Gaussian Quadrature //
  Quadrature gaussStif(eType, order - 1, 2);
  Quadrature gaussMass(eType, order,     2);

  const fullMatrix<double>& gCS = gaussStif.getPoints();
  const fullMatrix<double>& gCM = gaussMass.getPoints();

  // Functions //
  basis.preEvaluateDerivatives(gCS);
  basis.preEvaluateFunctions(gCM);

  // Local Terms //
  switch(form){
  case 0:{
    // Gradiends
    GroupOfJacobian jacM(domain, gCM, "jacobian");
    GroupOfJacobian jacS(domain, gCS, "invert");
    mass  = new TermFieldField<scalar>(jacM, basis, gaussMass);
    stif = new TermGradGrad<scalar>(  jacS, basis, gaussStif);
    break;
  }

  case 1:{
    // Curls //
    GroupOfJacobian jacM(domain, gCM, "invert");
    GroupOfJacobian jacS(domain, gCS, "jacobian");
    mass  = new TermGradGrad<scalar>(jacM, basis, gaussMass);
    stif = new TermCurlCurl<scalar>(jacS, basis, gaussStif);
    break;
  }

  default:
    // Else //
    throw Exception("FormulationSteadyWave does 0 and 1 forms only");
  }

  // No source //
  src = NULL;
}

template<typename scalar>
FormulationSteadyWave<scalar>::
FormulationSteadyWave(const GroupOfElement& domain,
                      const FunctionSpace& fs,
                      double k,
                      void (*nu)(fullVector<double>&, fullMatrix<scalar>&),
                      void (*eps)(fullVector<double>&, fullMatrix<scalar>&),
                      fullVector<scalar> (*source)(fullVector<double>&)){

  // Check domain stats: uniform mesh //
  std::pair<bool, size_t> uniform = domain.isUniform();
  size_t                    eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationSteadyWave needs a uniform mesh");

  // Save Data //
  ddomain = &domain;
  ffs     = &fs;

  // Wave Squared //
  kSquare = k * k;

  // Get Basis //
  const size_t form  = fs.getForm();
  const Basis& basis = fs.getBasis(eType);
  const size_t order = basis.getOrder();

  // Gaussian Quadrature //
  Quadrature gaussStif(eType, order - 1, 2);
  Quadrature gaussMass(eType, order,     2);

  const fullMatrix<double>& gCS = gaussStif.getPoints();
  const fullMatrix<double>& gCM = gaussMass.getPoints();

  // Functions //
  basis.preEvaluateDerivatives(gCS);
  basis.preEvaluateFunctions(gCM);

  // Local Terms //
  if(form != 1)
    throw Exception("FormulationSteadyWave with vector source uses 1 forms");

  GroupOfJacobian jacM(domain, gCM, "invert");
  GroupOfJacobian jacS(domain, gCS, "jacobian");
  stif = new TermCurlCurl<scalar>(jacS, basis, gaussStif, nu);
  mass = new TermGradGrad<scalar>(jacM, basis, gaussMass, eps);
  src  = new TermProjectionGrad<scalar>(jacM, basis, gaussMass, source);
}

template<typename scalar>
FormulationSteadyWave<scalar>::~FormulationSteadyWave(void){
  delete stif;
  delete mass;

  if(src)
    delete src;
}

template<typename scalar>
scalar FormulationSteadyWave<scalar>::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return
    stif->getTerm(dofI, dofJ, elementId) -
    mass->getTerm(dofI, dofJ, elementId) * kSquare;
}

template<typename scalar>
scalar FormulationSteadyWave<scalar>::
rhs(size_t equationI, size_t elementId) const{
  if(src)
    return src->getTerm(equationI, 0, elementId);
  else
    return 0;
}

template<typename scalar>
const FunctionSpace& FormulationSteadyWave<scalar>::field(void) const{
  return *ffs;
}

template<typename scalar>
const FunctionSpace& FormulationSteadyWave<scalar>::test(void) const{
  return *ffs;
}

template<typename scalar>
const GroupOfElement& FormulationSteadyWave<scalar>::domain(void) const{
  return *ddomain;
}

template<typename scalar>
bool FormulationSteadyWave<scalar>::isBlock(void) const{
  return true;
}
