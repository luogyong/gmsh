#include "ReferenceSpaceManager.h"
#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationOO2.h"

using namespace std;

FormulationOO2::
FormulationOO2(const GroupOfElement& domain,
               const FunctionSpaceScalar& fs,
               Complex a,
               Complex b,
               const std::map<Dof, Complex>& ddmDof){

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = domain.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationOO2 needs a uniform mesh");

  // a & b //
  this->a = a;
  this->b = b;

  // Domain //
  goe = &domain;

  // Save FunctionSpace & Get Basis //
  fspace = &fs;
  basis  = &fs.getBasis(eType);

  // Gaussian Quadrature (Field - Field) //
  Quadrature gauss(eType, basis->getOrder(), 2);

  gC = new fullMatrix<double>(gauss.getPoints());
  gW = new fullVector<double>(gauss.getWeights());

  // Gaussian Quadrature (Grad - Grad) -- Do not need to be saved //
  Quadrature gaussGG(eType, basis->getOrder() - 1, 2);

  const fullMatrix<double>& gCGG = gaussGG.getPoints();
  const fullVector<double>& gWGG = gaussGG.getWeights();

  // Pre-evalution //
  basis->preEvaluateFunctions(*gC);
  basis->preEvaluateDerivatives(gCGG);

  jac = new GroupOfJacobian(domain, *gC, "jacobian");
  GroupOfJacobian jacGG(domain, gCGG, "invert");

  // Local Terms //
  localTermsUU = new TermFieldField(*jac, *basis, *gW);
  localTermsGG = new TermGradGrad(jacGG, *basis, gWGG);

  // DDM //
  this->ddmDof = &ddmDof;
}

FormulationOO2::~FormulationOO2(void){
  delete localTermsUU;
  delete localTermsGG;

  delete gC;
  delete gW;
  delete jac;
}

Complex FormulationOO2::weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return
    a * localTermsUU->getTerm(dofI, dofJ, elementId) -
    b * localTermsGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOO2::rhs(size_t equationI, size_t elementId) const{
  // Init //
  double phi;
  double det;

  double pxyz[3];
  fullVector<double> xyz(3);

  Complex ddmValue;
  Complex integral = Complex(0, 0);

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

    // Compute ddmValue in the *physical* coordinate
    ReferenceSpaceManager::mapFromABCtoXYZ(element,
                                           (*gC)(g, 0),
                                           (*gC)(g, 1),
                                           (*gC)(g, 2),
                                           pxyz);
    xyz(0) = pxyz[0];
    xyz(1) = pxyz[1];
    xyz(2) = pxyz[2];

    ddmValue = interpolate(element, xyz);

    // Integrate
    integral += ddmValue * phi * fabs(det) * (*gW)(g);
  }

  return integral;
}

Complex FormulationOO2::
interpolate(const MElement& element, const fullVector<double>& xyz) const{
  // Get Dofs associated to element //
  const vector<Dof>  dof = fspace->getKeys(element);
  const size_t      nDof = dof.size();

  // Get Values of these Dofs //
  map<Dof, Complex>::const_iterator end = ddmDof->end();
  map<Dof, Complex>::const_iterator it;
  vector<double> realCoef(nDof);
  vector<double> imagCoef(nDof);

  for(size_t i = 0; i < nDof; i++){
    it = ddmDof->find(dof[i]);
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

const FunctionSpace& FormulationOO2::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationOO2::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationOO2::domain(void) const{
  return *goe;
}
