#include "ReferenceSpaceManager.h"
#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"
#include "Mapper.h"

#include "Exception.h"
#include "FormulationUpdateOO2.h"

using namespace std;

FormulationUpdateOO2::
FormulationUpdateOO2(const GroupOfElement& domain,
                     const FunctionSpaceScalar& fs,
                     Complex a,
                     Complex b,
                     const std::map<Dof, Complex>& solution,
                     const std::map<Dof, Complex>& oldG){
  // Domain //
  goe = &domain;

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationUpdateOO2 needs a uniform mesh");

  // Save fspace //
  fspace = &fs;
  basis  = &fs.getBasis(eType);

  // a & b //
  this->a = a;
  this->b = b;

  // Gaussian Quadrature (Field - Field) //
  Quadrature gaussFF(eType, basis->getOrder(), 2);

  gCFF = new fullMatrix<double>(gaussFF.getPoints());
  gWFF = new fullVector<double>(gaussFF.getWeights());

  // Gaussian Quadrature (Grad - Grad) //
  Quadrature gaussGG(eType, basis->getOrder() - 1, 2);

  gCGG = new fullMatrix<double>(gaussGG.getPoints());
  gWGG = new fullVector<double>(gaussGG.getWeights());

  // Pre-evalution //
  basis->preEvaluateFunctions(*gCFF);
  basis->preEvaluateDerivatives(*gCGG);

  jacFF = new GroupOfJacobian(domain, *gCFF, "jacobian");
  jacGG = new GroupOfJacobian(domain, *gCGG, "invert");

  // DDM //
  this->solution = &solution;
  this->oldG     = &oldG;
}

FormulationUpdateOO2::~FormulationUpdateOO2(void){
  delete gCFF;
  delete gWFF;
  delete jacFF;

  delete gCGG;
  delete gWGG;
  delete jacGG;
}

Complex FormulationUpdateOO2::
weak(size_t dofI, size_t dofJ, size_t elementId) const{
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
    jacFF->getJacobian(elementId).getJacobianMatrix();

  // Loop over Integration Point //
  const size_t G = gWFF->size();

  for(size_t g = 0; g < G; g++){
    det   = allJac[g]->second;

    phiI = eFun(dofI, g);
    phiJ = eFun(dofJ, g);

    integral += phiI * phiJ * fabs(det) * (*gWFF)(g);
  }

  return Complex(integral, 0);
}

Complex FormulationUpdateOO2::rhs(size_t equationI, size_t elementId) const{
  // Init //
  size_t G;

  double det;
  const fullMatrix<double>* jac;

  double phi;
  fullVector<double> gradPhi(3);

  double pxyz[3];
  fullVector<double>  xyz(3);
  fullVector<Complex> gradValue;

  Complex oldGValue;
  Complex solutionValue;
  Complex sub;

  Complex integral = Complex(0, 0);

  // Get Element //
  const MElement& element = goe->get(elementId);

  // Get Basis Functions //
  const fullMatrix<double>& eFunFF =
    basis->getPreEvaluatedFunctions(element);

  // Get Grad Basis Functions //
  const fullMatrix<double>& eFunGG =
    basis->getPreEvaluatedDerivatives(element);

  // Get Jacobians //
  const vector<const pair<const fullMatrix<double>*, double>*>& allJacFF =
    jacFF->getJacobian(elementId).getJacobianMatrix();

  const vector<const pair<const fullMatrix<double>*, double>*>& allJacGG =
    jacGG->getJacobian(elementId).getInvertJacobianMatrix();

  // Loop over Integration Point (Field - Field) //
  G = gWFF->size();

  for(size_t g = 0; g < G; g++){
    // Compute phi
    det = allJacFF[g]->second;
    phi = eFunFF(equationI, g);

    // Get *physical* coordinate
    ReferenceSpaceManager::mapFromABCtoXYZ(element,
                                           (*gCFF)(g, 0),
                                           (*gCFF)(g, 1),
                                           (*gCFF)(g, 2),
                                           pxyz);
    xyz(0) = pxyz[0];
    xyz(1) = pxyz[1];
    xyz(2) = pxyz[2];

    // OldG & solution in *physical* coordinate
    oldGValue     = interpolate(element, xyz, *oldG);
    solutionValue = interpolate(element, xyz, *solution);
    sub           = Complex(2, 0) * a * solutionValue - oldGValue;

    // Integrate
    integral += sub * phi * fabs(det) * (*gWFF)(g);
  }

  // Loop over Integration Point (Grad - Grad) //
  G = gWGG->size();

  for(size_t g = 0; g < G; g++){
    // Compute gradPhi
    det = allJacGG[g]->second;
    jac = allJacGG[g]->first;

    Mapper::hCurl(eFunGG, equationI, g, *jac, gradPhi);

    // Get *physical* coordinate
    ReferenceSpaceManager::mapFromABCtoXYZ(element,
                                           (*gCGG)(g, 0),
                                           (*gCGG)(g, 1),
                                           (*gCGG)(g, 2),
                                           pxyz);
    xyz(0) = pxyz[0];
    xyz(1) = pxyz[1];
    xyz(2) = pxyz[2];

    // grad(solution) in *physical* coordinate
    gradValue = interpolateGrad(element, xyz, *solution);

    // Integrate
    integral +=
      Complex(-2, 0) * b *
      (gradValue(0) * Complex(gradPhi(0), 0) +
       gradValue(1) * Complex(gradPhi(1), 0) +
       gradValue(2) * Complex(gradPhi(2), 0)) * fabs(det) * (*gWGG)(g);
  }

  return integral;
}

Complex FormulationUpdateOO2::
interpolate(const MElement& element,
            const fullVector<double>& xyz,
            const std::map<Dof, Complex>& f) const{

  // Get Dofs associated to element //
  const vector<Dof>  dof = fspace->getKeys(element);
  const size_t      nDof = dof.size();

  // Get Values of these Dofs //
  map<Dof, Complex>::const_iterator end = f.end();
  map<Dof, Complex>::const_iterator it;
  vector<double> realCoef(nDof);
  vector<double> imagCoef(nDof);

  for(size_t i = 0; i < nDof; i++){
    it = f.find(dof[i]);
    if(it == end)
      throw Exception("Snif");

    realCoef[i] = it->second.real();
    imagCoef[i] = it->second.imag();
  }

  // Interpolate
  double real = fspace->interpolate(element, realCoef, xyz);
  double imag = fspace->interpolate(element, imagCoef, xyz);

  return Complex(real, imag);
}

fullVector<Complex> FormulationUpdateOO2::
interpolateGrad(const MElement& element,
                const fullVector<double>& xyz,
                const std::map<Dof, Complex>& f) const{

  // Get Dofs associated to element //
  const vector<Dof>  dof = fspace->getKeys(element);
  const size_t      nDof = dof.size();

  // Get Values of these Dofs //
  map<Dof, Complex>::const_iterator end = f.end();
  map<Dof, Complex>::const_iterator it;
  vector<double> realCoef(nDof);
  vector<double> imagCoef(nDof);

  for(size_t i = 0; i < nDof; i++){
    it = f.find(dof[i]);
    if(it == end)
      throw Exception("Snif");

    realCoef[i] = it->second.real();
    imagCoef[i] = it->second.imag();
  }

  // Interpolate
  fullVector<double> re = fspace->interpolateDerivative(element, realCoef, xyz);
  fullVector<double> im = fspace->interpolateDerivative(element, imagCoef, xyz);

  // Return //
  if(re.size() != 3)
    throw Exception("Snif");

  fullVector<Complex> ret(3);

  ret(0) = Complex(re(0), im(0));
  ret(1) = Complex(re(1), im(1));
  ret(2) = Complex(re(2), im(2));

  return ret;
}

const FunctionSpace& FormulationUpdateOO2::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationUpdateOO2::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationUpdateOO2::domain(void) const{
  return *goe;
}
