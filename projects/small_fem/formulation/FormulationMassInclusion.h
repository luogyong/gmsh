////////////////////////////////////////////////////
// Templates Implementations for FormulationMass: //
// Inclusion compilation model                    //
//                                                //
// Damn you gcc: we want 'export' !               //
////////////////////////////////////////////////////

#include "GroupOfJacobian.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"
#include "Quadrature.h"

template<typename scalar>
FormulationMass<scalar>::FormulationMass(const GroupOfElement& domain,
                                         const FunctionSpace& field,
                                         const FunctionSpace& test){
  // Save Data //
  ddomain = &domain;
  ffield  = &field;
  ttest   = &test;

  // Check domain stats: uniform mesh //
  std::pair<bool, size_t> uniform = domain.isUniform();
  size_t                    eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationMass needs a uniform mesh");

  // Check if same forms //
  const size_t form = field.getForm();
  if(form != test.getForm())
    throw Exception
      ("FormulationMass needs function spaces of the same form");

  // Get Basis //
  const Basis& basis = test.getBasis(eType);
  const size_t order = basis.getOrder();

  // Gaussian Quadrature //
  Quadrature gauss(eType, order, 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Functions //
  basis.preEvaluateFunctions(gC);

  // Local Terms //
  switch(form){
  // Gradiends
  case 0:{
    GroupOfJacobian jac(domain, gC, "jacobian");
    localTerms = new TermFieldField<scalar>(jac, basis, gauss);
    break;
  }

  // Curls //
  case 1:{
    GroupOfJacobian jac(domain, gC, "invert");
    localTerms = new TermGradGrad<scalar>(jac, basis, gauss);
    break;
  }

  // Else //
  default:
    throw Exception("FormulationMass does 0 and 1 forms only");
  }
}

template<typename scalar>
FormulationMass<scalar>::
FormulationMass(const GroupOfElement& domain,
                const FunctionSpace& field,
                const FunctionSpace& test,
                void (*m)(fullVector<double>&, fullMatrix<scalar>&)){
  // Save Data //
  ddomain = &domain;
  ffield  = &field;
  ttest   = &test;

  // Check domain stats: uniform mesh //
  std::pair<bool, size_t> uniform = domain.isUniform();
  size_t                    eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationMass needs a uniform mesh");

  // Check if same forms //
  const size_t form = field.getForm();
  if(form != test.getForm())
    throw Exception
      ("FormulationMass needs function spaces of the same form");

  if(form != 1)
    throw Exception
      ("FormulationMass needs a 1 form FunctionSpace with tensorial materials");

  // Get Basis //
  const Basis& basis = test.getBasis(eType);
  const size_t order = basis.getOrder();

  // Gaussian Quadrature //
  Quadrature gauss(eType, order, 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Functions //
  basis.preEvaluateFunctions(gC);

  // Local Terms //
  GroupOfJacobian jac(domain, gC, "invert");
  localTerms = new TermGradGrad<scalar>(jac, basis, gauss, m);
}

template<typename scalar>
FormulationMass<scalar>::
FormulationMass(const GroupOfElement& domain,
                const FunctionSpace& field,
                const FunctionSpace& test,
                scalar (*m)(fullVector<double>&)){
  // Save Data //
  ddomain = &domain;
  ffield  = &field;
  ttest   = &test;

  // Check domain stats: uniform mesh //
  std::pair<bool, size_t> uniform = domain.isUniform();
  size_t                    eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationMass needs a uniform mesh");

  // Check if same forms //
  const size_t form = field.getForm();
  if(form != test.getForm())
    throw Exception
      ("FormulationMass needs function spaces of the same form");

  if(form != 0)
    throw Exception
      ("FormulationMass needs a 0 form FunctionSpace with scalar materials");

  // Get Basis //
  const Basis& basis = test.getBasis(eType);
  const size_t order = basis.getOrder();

  // Gaussian Quadrature //
  Quadrature gauss(eType, order, 2);
  const fullMatrix<double>& gC = gauss.getPoints();

  // Functions //
  basis.preEvaluateFunctions(gC);

  // Local Terms //
  GroupOfJacobian jac(domain, gC, "jacobian");
  localTerms = new TermFieldField<scalar>(jac, basis, gauss, m);
}

template<typename scalar>
FormulationMass<scalar>::~FormulationMass(void){
  delete localTerms;
}

template<typename scalar>
scalar FormulationMass<scalar>::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  return localTerms->getTerm(dofI, dofJ, elementId);
}

template<typename scalar>
scalar FormulationMass<scalar>::
rhs(size_t equationI, size_t elementId) const{

  return 0;
}

template<typename scalar>
const FunctionSpace& FormulationMass<scalar>::field(void) const{
  return *ffield;
}

template<typename scalar>
const FunctionSpace& FormulationMass<scalar>::test(void) const{
  return *ttest;
}

template<typename scalar>
const GroupOfElement& FormulationMass<scalar>::domain(void) const{
  return *ddomain;
}

template<typename scalar>
bool FormulationMass<scalar>::isBlock(void) const{
  return true;
}
