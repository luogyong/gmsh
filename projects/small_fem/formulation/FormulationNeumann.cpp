#include "BasisGenerator.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"

#include "FormulationNeumann.h"

using namespace std;

FormulationNeumann::FormulationNeumann(GroupOfElement& goe,
                                       const FunctionSpaceScalar& fs,
                                       double k){

  // Check GroupOfElement Stats: Uniform Mesh //
  const vector<size_t>& gType = goe.getTypeStats();
  const size_t nGType = gType.size();
  size_t eType = (size_t)(-1);

  for(size_t i = 0; i < nGType; i++)
    if((eType == (size_t)(-1)) && (gType[i] != 0))
      eType = i;
    else if((eType != (size_t)(-1)) && (gType[i] != 0))
      throw Exception("FormulationNeumann needs a uniform mesh");

  // Wavenumber //
  this->k = k;

  // Save FunctionSpace & Get Basis //
  const Basis& basis = fs.getBasis(eType);
  const size_t order = basis.getOrder();
  fspace             = &fs;

  // Gaussian Quadrature //
  Quadrature gaussFF(eType, order, 2);
  const fullMatrix<double>& gC = gaussFF.getPoints();
  const fullVector<double>& gW = gaussFF.getWeights();

  // Local Terms //
  basis.preEvaluateFunctions(gC);

  GroupOfJacobian jac(goe, gC, "jacobian");

  localTerms = new TermFieldField(jac, basis, gW);
}

FormulationNeumann::~FormulationNeumann(void){
  delete localTerms;
}

complex<double> FormulationNeumann::weak(size_t dofI, size_t dofJ,
                                         size_t elementId) const{
  return
    complex<double>(0, -1 * k * localTerms->getTerm(dofI, dofJ, elementId));
}

complex<double> FormulationNeumann::rhs(size_t equationI,
                                        size_t elementId) const{
  return complex<double>(0, 0);
}

bool FormulationNeumann::isGeneral(void) const{
  return false;
}

complex<double> FormulationNeumann::weakB(size_t dofI, size_t dofJ,
                                          size_t elementId) const{
  return complex<double>(0, 0);
}


const FunctionSpace& FormulationNeumann::fs(void) const{
  return *fspace;
}
