#include <cmath>

#include "Quadrature.h"
#include "Mapper.h"

#include "ReferenceSpaceManager.h"
#include "FormulationSteadySlow.h"

using namespace std;

FormulationSteadySlow::FormulationSteadySlow(const GroupOfElement& domain,
                                             const FunctionSpaceVector& fs,
                                             double k){
  // Dummy materials //
  nu  = dummy;
  eps = dummy;

  // Init //
  init(domain, fs, k);
}

FormulationSteadySlow::
FormulationSteadySlow(const GroupOfElement& domain,
                      const FunctionSpaceVector& fs,
                      double k,
                      void (*nu)(fullVector<double>&, fullMatrix<double>&),
                      void (*eps)(fullVector<double>&, fullMatrix<double>&)){
  // Materials //
  this->nu  = nu;
  this->eps = eps;

  // Init //
  init(domain, fs, k);
}

void FormulationSteadySlow::init(const GroupOfElement& domain,
                                 const FunctionSpaceVector& fs,
                                 double k){

  // Check domain stats: uniform mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationSteadySlow needs a uniform mesh");

  // Wave Squared //
  kSquare = k * k;

  // Domain //
  ddomain = &domain;

  // Save FunctionSpace & Get Basis //
  fspace = &fs;
  basis  = &fs.getBasis(eType);

  // Gaussian Quadrature //
  const size_t order = basis->getOrder();

  Quadrature gaussCurlCurl(eType, order - 1, 2);
  Quadrature gaussFF(eType, order, 2);

  gC1 = new fullMatrix<double>(gaussCurlCurl.getPoints());
  gW1 = new fullVector<double>(gaussCurlCurl.getWeights());

  gC2 = new fullMatrix<double>(gaussFF.getPoints());
  gW2 = new fullVector<double>(gaussFF.getWeights());

  G1 = gW1->size();
  G2 = gW2->size();

  // PreEvaluate
  basis->preEvaluateDerivatives(*gC1);
  basis->preEvaluateFunctions(*gC2);

  jac1 = new GroupOfJacobian(domain, *gC1, "jacobian");
  jac2 = new GroupOfJacobian(domain, *gC2, "invert");
}

FormulationSteadySlow::~FormulationSteadySlow(void){
  delete gC1;
  delete gW1;
  delete gC2;
  delete gW2;
  delete jac1;
  delete jac2;
}

double FormulationSteadySlow::
weak(size_t dofI, size_t dofJ, size_t elementId) const{

  // Init Some Stuff //
  fullVector<double> phiI(3);
  fullVector<double> phiK(3);
  fullVector<double> phiJ(3);

  fullMatrix<double> T(3, 3);

  fullVector<double> xyz(3);
  double pxyz[3];

  const fullMatrix<double>* jac;

  double integral1 = 0;
  double integral2 = 0;
  double det;

  // Get Element //
  const MElement& element = ddomain->get(elementId);

  // Get Basis Functions //
  const fullMatrix<double>& eCurlFun =
    basis->getPreEvaluatedDerivatives(element);

  const fullMatrix<double>& eFun =
    basis->getPreEvaluatedFunctions(element);

  // Get Jacobians //
  const vector<const pair<const fullMatrix<double>*, double>*>& MJac =
    jac1->getJacobian(elementId).getJacobianMatrix();

  const vector<const pair<const fullMatrix<double>*, double>*>& invJac =
    jac2->getJacobian(elementId).getInvertJacobianMatrix();

  // Loop over Integration Point (Term 1) //
  for(int g = 0; g < G1; g++){
    det = MJac[g]->second;
    jac = MJac[g]->first;

    Mapper::hDiv(eCurlFun, dofI, g, *jac, det, phiI);
    Mapper::hDiv(eCurlFun, dofJ, g, *jac, det, phiJ);

    // Compute Nu in the *physical* coordinate
    ReferenceSpaceManager::
      mapFromABCtoXYZ(element, (*gC1)(g, 0), (*gC1)(g, 1), (*gC1)(g, 2), pxyz);

    xyz(0) = pxyz[0];
    xyz(1) = pxyz[1];
    xyz(2) = pxyz[2];

    nu(xyz, T); // T = Nu(x, y, z);

    // phiK = Nu * phiI
    T.mult(phiI, phiK);

    integral1 +=
      ((phiK * phiJ)) * fabs(det) * (*gW1)(g);
  }


  // Loop over Integration Point (Term 2) //
  for(int g = 0; g < G2; g++){
    det = invJac[g]->second;
    jac = invJac[g]->first;

    Mapper::hCurl(eFun, dofI, g, *jac, phiI);
    Mapper::hCurl(eFun, dofJ, g, *jac, phiJ);

    // Compute Epsilon in the *physical* coordinate
    ReferenceSpaceManager::
      mapFromABCtoXYZ(element, (*gC2)(g, 0), (*gC2)(g, 1), (*gC2)(g, 2), pxyz);

    xyz(0) = pxyz[0];
    xyz(1) = pxyz[1];
    xyz(2) = pxyz[2];

    eps(xyz, T); // T = Eps(x, y, z);

    // phiK = Eps * phiI
    T.mult(phiI, phiK);

    integral2 +=
      ((phiK * phiJ) * kSquare) * fabs(det) * (*gW2)(g);
  }

  return integral1 - integral2;
}

double FormulationSteadySlow::rhs(size_t equationI, size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationSteadySlow::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationSteadySlow::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationSteadySlow::domain(void) const{
  return *ddomain;
}

bool FormulationSteadySlow::isBlock(void) const{
  return true;
}

void FormulationSteadySlow::dummy(fullVector<double>& xyz,
                                  fullMatrix<double>& tensor){
  tensor.scale(0);
  tensor(0, 0) = 1;
  tensor(1, 1) = 1;
  tensor(2, 2) = 1;
}
