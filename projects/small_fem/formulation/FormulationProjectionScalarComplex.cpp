#include "ReferenceSpaceManager.h"
#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "GroupOfElement.h"
#include "Quadrature.h"

#include "FormulationProjectionScalar.h"

using namespace std;

template<>
FormulationProjectionScalar<complex<double> >::
FormulationProjectionScalar(const GroupOfElement& goe,
                            const FunctionSpaceScalar& fs,
                            std::complex<double> (*f)(fullVector<double>& xyz)){
  // Save Domain //
  this->goe = &goe;

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = goe.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception
      ("FormulationProjectionScalar<complex> needs a uniform mesh");

  // Save fspace //
  fspace = &fs;
  basis  = &fs.getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis->getOrder(), 2);

  gC = new fullMatrix<double>(gauss.getPoints());
  gW = new fullVector<double>(gauss.getWeights());

  // Pre-evalution //
  basis->preEvaluateFunctions(*gC);
  jac = new GroupOfJacobian(goe, *gC, "jacobian");

  // f //
  this->f = f;
}

template<>
FormulationProjectionScalar<complex<double> >::
~FormulationProjectionScalar(void){
  delete gC;
  delete gW;
  delete jac;
}

template<>
complex<double> FormulationProjectionScalar<complex<double> >::
weak(size_t dofI, size_t dofJ,size_t elementId) const{

  // Init //
  double phiI;
  double phiJ;

  double det;
  double integral = 0;

  // Get Element //
  const MElement& element = goe->get(elementId);

  // Get Basis Functions //
  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Get Jacobians //
  const vector<const pair<const fullMatrix<double>*, double>*>& allJac =
    jac->getJacobian(elementId).getJacobianMatrix();

  // Loop over Integration Point //
  const size_t G = gW->size();

  for(size_t g = 0; g < G; g++){
    det   = allJac[g]->second;

    phiI = eFun(dofI, g);
    phiJ = eFun(dofJ, g);

    integral += phiI * phiJ * fabs(det) * (*gW)(g);
  }

  return complex<double>(integral, 0);
}

template<>
complex<double> FormulationProjectionScalar<complex<double> >::
rhs(size_t equationI, size_t elementId) const{

  // Init //
  double phi;
  double det;

  fullVector<double> xyz(3);
  double             pxyz[3];
  complex<double>    fxyz;

  complex<double>    integral = complex<double>(0, 0);

  // Get Element //
  const MElement& element = goe->get(elementId);

  // Get Basis Functions //
  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Get Jacobians //
  const vector<const pair<const fullMatrix<double>*, double>*>& allJac =
    jac->getJacobian(elementId).getJacobianMatrix();

  // Loop over Integration Point //
  const size_t G = gW->size();

  for(size_t g = 0; g < G; g++){
    // Compute phi
    det = allJac[g]->second;
    phi = eFun(equationI, g);

    // Compute f in the *physical* coordinate
    ReferenceSpaceManager::mapFromABCtoXYZ(element,
                                           (*gC)(g, 0),
                                           (*gC)(g, 1),
                                           (*gC)(g, 2),
                                           pxyz);
    xyz(0) = pxyz[0];
    xyz(1) = pxyz[1];
    xyz(2) = pxyz[2];

    fxyz = f(xyz);

    // Integrate
    integral += fxyz * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}

template<>
const FunctionSpace& FormulationProjectionScalar<complex<double> >::
field(void) const{
  return *fspace;
}

template<>
const FunctionSpace& FormulationProjectionScalar<complex<double> >::
test(void) const{
  return *fspace;
}

template<>
const GroupOfElement& FormulationProjectionScalar<complex<double> >::
domain(void) const{
  return *goe;
}
