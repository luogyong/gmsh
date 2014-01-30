#include "ReferenceSpaceManager.h"
#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationOO2.h"

using namespace std;

FormulationOO2::
FormulationOO2(GroupOfElement& goe,
               const FunctionSpaceScalar& fs,
               std::complex<double> a,
               std::complex<double> b,
               const std::map<Dof, std::complex<double> >& ddmDof){

  // Check GroupOfElement Stats: Uniform Mesh //
  const vector<size_t>& gType = goe.getTypeStats();
  const size_t nGType = gType.size();
  size_t eType = (size_t)(-1);

  for(size_t i = 0; i < nGType; i++)
    if((eType == (size_t)(-1)) && (gType[i] != 0))
      eType = i;
    else if((eType != (size_t)(-1)) && (gType[i] != 0))
      throw Exception("FormulationOO2 needs a uniform mesh");

  // a & b //
  this->a = a;
  this->b = b;

  // Domain //
  this->goe = &goe;

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

  jac = new GroupOfJacobian(goe, *gC, "jacobian");
  GroupOfJacobian jacGG(goe, gCGG, "invert");

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

complex<double> FormulationOO2::weak(size_t dofI, size_t dofJ,
                                     size_t elementId) const{
  return
    a * localTermsUU->getTerm(dofI, dofJ, elementId) -
    b * localTermsGG->getTerm(dofI, dofJ, elementId);
}

complex<double> FormulationOO2::rhs(size_t equationI, size_t elementId) const{
  // Init //
  double phi;
  double det;

  double pxyz[3];
  fullVector<double> xyz(3);

  complex<double> ddmValue;
  complex<double> integral = complex<double>(0, 0);

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

std::complex<double> FormulationOO2::
interpolate(const MElement& element, const fullVector<double>& xyz) const{
  // Get Dofs associated to element //
  const vector<Dof>  dof = fspace->getKeys(element);
  const size_t      nDof = dof.size();

  // Get Values of these Dofs //
  map<Dof, complex<double> >::const_iterator end = ddmDof->end();
  map<Dof, complex<double> >::const_iterator it;
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

  return complex<double>(real, imag);
}

bool FormulationOO2::isGeneral(void) const{
  return false;
}

complex<double>  FormulationOO2::weakB(size_t dofI, size_t dofJ,
                                        size_t elementId) const{
  return complex<double>(0, 0);
}

const FunctionSpace& FormulationOO2::fs(void) const{
  return *fspace;
}
