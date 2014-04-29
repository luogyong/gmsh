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

  // Save FunctionSpace & Domain //
  fspace = &fs;
  goe    = &domain;

  // Get Basis //
  const Basis& basis = fs.getBasis(eType);

  // Gaussian Quadrature (Field - Field & Field - Projection) //
  Quadrature gaussFF(eType, basis.getOrder(), 2);
  const fullMatrix<double>& gCFF = gaussFF.getPoints();

  // Gaussian Quadrature (Grad - Grad) //
  Quadrature gaussGG(eType, basis.getOrder() - 1, 2);
  const fullMatrix<double>& gCGG = gaussGG.getPoints();

  // Pre-evalution //
  basis.preEvaluateFunctions(gCFF);
  basis.preEvaluateDerivatives(gCGG);

  // Jacobians //
  GroupOfJacobian jacFF(domain, gCFF, "jacobian");
  GroupOfJacobian jacGG(domain, gCGG, "invert");

  // Local Terms //
  localTermsFF = new TermFieldField<double>(jacFF, basis, gaussFF);
  localTermsGG = new TermGradGrad<double>(jacGG, basis, gaussGG);
  localTermsPr =
    new TermProjectionField<Complex>(jacFF, basis, gaussFF, fs, ddmDof);
}

FormulationOO2::~FormulationOO2(void){
  delete localTermsFF;
  delete localTermsGG;
  delete localTermsPr;
}

Complex FormulationOO2::weak(size_t dofI, size_t dofJ, size_t elementId) const{
  return
    a * localTermsFF->getTerm(dofI, dofJ, elementId) -
    b * localTermsGG->getTerm(dofI, dofJ, elementId);
}

Complex FormulationOO2::rhs(size_t equationI, size_t elementId) const{
  return localTermsPr->getTerm(0, equationI, elementId);
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

bool FormulationOO2::isBlock(void) const{
  return true;
}
