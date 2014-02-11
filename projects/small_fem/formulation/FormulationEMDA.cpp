#include "ReferenceSpaceManager.h"
#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationEMDA.h"

using namespace std;

FormulationEMDA::
FormulationEMDA(const GroupOfElement& goe,
                const FunctionSpaceScalar& fs,
                double k,
                double chi,
                const std::map<Dof, std::complex<double> >& ddmDof){

  // Check GroupOfElement Stats: Uniform Mesh //
  pair<bool, size_t> uniform = goe.isUniform();
  size_t               eType = uniform.second;

  if(!uniform.first)
    throw Exception("FormulationEMDA needs a uniform mesh");

  // Wavenumber & Chi //
  this->k   = k;
  this->chi = chi;

  // Domain //
  this->goe = &goe;

  // Save FunctionSpace & Get Basis //
  fspace = &fs;
  basis  = &fs.getBasis(eType);

  // Gaussian Quadrature //
  Quadrature gauss(eType, basis->getOrder(), 2);

  gC = new fullMatrix<double>(gauss.getPoints());
  gW = new fullVector<double>(gauss.getWeights());

  // Pre-evalution //
  basis->preEvaluateFunctions(*gC);
  jac = new GroupOfJacobian(goe, *gC, "jacobian");

  // Local Terms //
  localTerms = new TermFieldField(*jac, *basis, *gW);

  // DDM //
  this->ddmDof = &ddmDof;
}

FormulationEMDA::~FormulationEMDA(void){
  delete localTerms;

  delete gC;
  delete gW;
  delete jac;
}

complex<double> FormulationEMDA::weak(size_t dofI, size_t dofJ,
                                      size_t elementId) const{
  return
    complex<double>(chi, -k) * localTerms->getTerm(dofI, dofJ, elementId);
}

complex<double> FormulationEMDA::rhs(size_t equationI, size_t elementId) const{
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

std::complex<double> FormulationEMDA::
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

const FunctionSpace& FormulationEMDA::field(void) const{
  return *fspace;
}

const FunctionSpace& FormulationEMDA::test(void) const{
  return *fspace;
}

const GroupOfElement& FormulationEMDA::domain(void) const{
  return *goe;
}
