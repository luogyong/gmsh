#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "Exception.h"
#include "FormulationPoisson.h"

using namespace std;

// Poisson //
FormulationPoisson::FormulationPoisson(GroupOfElement& goe,
                                       FunctionSpaceScalar& fs,
                                       double (*f)(fullVector<double>& xyz)){

  // Check GroupOfElement Stats: Uniform Mesh //
  const vector<size_t>& gType = goe.getTypeStats();
  const size_t nGType = gType.size();
  size_t eType = (size_t)(-1);

  for(size_t i = 0; i < nGType; i++)
    if((eType == (size_t)(-1)) && (gType[i] != 0))
      eType = i;
    else if((eType != (size_t)(-1)) && (gType[i] != 0))
      throw Exception("FormulationPoisson needs a uniform mesh");

  // Save FunctionSpace (TO BE REMOVED) & Get Basis //
  const Basis& basis = fs.getBasis(eType);
  const size_t order = basis.getOrder();
  self               = false;
  fspace             = &fs;

  // Source Term //
  fSource = f;

  // Gaussian Quadrature //
  Quadrature gaussGradGrad(eType, order - 1, 2);
  Quadrature gaussFF(eType, order, 2);

  const fullMatrix<double>& gCL = gaussGradGrad.getPoints();
  const fullVector<double>& gWL = gaussGradGrad.getWeights();

  const fullMatrix<double>& gCR = gaussFF.getPoints();
  const fullVector<double>& gWR = gaussFF.getWeights();

  // Local Terms //
  basis.preEvaluateDerivatives(gCL);
  basis.preEvaluateFunctions(gCR);

  GroupOfJacobian jacL(goe, gCL, "invert");
  GroupOfJacobian jacR(goe, gCR, "jacobian");

  localTermsL = new TermGradGrad(jacL, basis, gWL);
  localTermsR = new TermProjectionField(jacR, basis, gWR, gCR, fSource);
}

FormulationPoisson::FormulationPoisson(GroupOfElement& goe,
                                       double (*f)(fullVector<double>& xyz),
                                       size_t order){
  // Can't have 0th order //
  if(order == 0)
    throw
      Exception("Can't have a Poisson formulation of order 0");

  // Function Space & Basis //
  self   = true;
  basis  = BasisGenerator::generate(goe.get(0).getType(),
                                    0, order, "hierarchical");

  fspace = new FunctionSpaceScalar(goe, *basis);

  // Source Term //
  fSource = f;

  // Gaussian Quadrature //
  Quadrature gaussGradGrad(goe.get(0).getType(), order - 1, 2);
  Quadrature gaussFF(goe.get(0).getType(), order, 2);

  const fullMatrix<double>& gCL = gaussGradGrad.getPoints();
  const fullVector<double>& gWL = gaussGradGrad.getWeights();

  const fullMatrix<double>& gCR = gaussFF.getPoints();
  const fullVector<double>& gWR = gaussFF.getWeights();

  // Local Terms //
  basis->preEvaluateDerivatives(gCL);
  basis->preEvaluateFunctions(gCR);

  GroupOfJacobian jacL(goe, gCL, "invert");
  GroupOfJacobian jacR(goe, gCR, "jacobian");

  localTermsL = new TermGradGrad(jacL, *basis, gWL);
  localTermsR = new TermProjectionField(jacR, *basis, gWR, gCR, fSource);
}

FormulationPoisson::~FormulationPoisson(void){
  if(self){
    delete basis;
    delete fspace;
  }

  delete localTermsL;
  delete localTermsR;
}

double FormulationPoisson::weak(size_t dofI, size_t dofJ,
                                size_t elementId) const{

  return localTermsL->getTerm(dofI, dofJ, elementId);
}

double FormulationPoisson::rhs(size_t equationI,
                               size_t elementId) const{

  return localTermsR->getTerm(0, equationI, elementId);
}

bool FormulationPoisson::isGeneral(void) const{
  return false;
}

double FormulationPoisson::weakB(size_t dofI,
                                 size_t dofJ,
                                 size_t elementId) const{
  return 0;
}

const FunctionSpace& FormulationPoisson::fs(void) const{
  return *fspace;
}
