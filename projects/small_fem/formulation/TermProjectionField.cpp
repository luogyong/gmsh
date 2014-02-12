#include <vector>

#include "TermProjectionField.h"
#include "SmallFem.h"

using namespace std;

// Real Version //
// ------------ //
template<>
void TermProjectionField<double>::computeC(const Basis& basis,
                                           const fullVector<double>& gW,
                                           fullMatrix<double>**& cM){
  const size_t nG = gW.size();

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(nG, nFunction);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis.getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(size_t g = 0; g < nG; g++)
      for(size_t i = 0; i < nFunction; i++)
        (*cM[s])(g, i) = gW(g) * phi(i, g);
  }
}

template<>
double TermProjectionField<double>::
interpolate(const MElement& element,
            const fullVector<double>& xyz) const{

  // Get Dofs associated to element //
  const vector<Dof>  dof = fsScalar->getKeys(element);
  const size_t      nDof = dof.size();

  // Get Values of these Dofs //
  map<Dof, double>::const_iterator end = dofValue->end();
  map<Dof, double>::const_iterator it;
  vector<double> coef(nDof);

  for(size_t i = 0; i < nDof; i++){
    it = dofValue->find(dof[i]);
    if(it == end)
      throw Exception("TermProjectionField::interpolate unknown dof %s",
                      dof[i].toString().c_str());

    coef[i] = it->second;
  }

  // Interpolate & Return
  return fsScalar->interpolate(element, coef, xyz);
}

// Complex Version //
// --------------- //
template<>
void TermProjectionField<Complex>::computeC(const Basis& basis,
                                            const fullVector<double>& gW,
                                            fullMatrix<Complex>**& cM){
  const size_t nG = gW.size();

  // Alloc //
  cM = new fullMatrix<Complex>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<Complex>(nG, nFunction);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      basis.getPreEvaluatedFunctions(s);

    // Loop on Gauss Points
    for(size_t g = 0; g < nG; g++)
      for(size_t i = 0; i < nFunction; i++)
        (*cM[s])(g, i) = Complex(gW(g) * phi(i, g), 0);
  }
}

template<>
Complex TermProjectionField<Complex>::
interpolate(const MElement& element,
            const fullVector<double>& xyz) const{

  // Get Dofs associated to element //
  const vector<Dof>  dof = fsScalar->getKeys(element);
  const size_t      nDof = dof.size();

  // Get Values of these Dofs //
  map<Dof, Complex>::const_iterator end = dofValue->end();
  map<Dof, Complex>::const_iterator it;
  vector<double> realCoef(nDof);
  vector<double> imagCoef(nDof);

  for(size_t i = 0; i < nDof; i++){
    it = dofValue->find(dof[i]);
    if(it == end)
      throw Exception("TermProjectionField::interpolate unknown dof %s",
                      dof[i].toString().c_str());

    realCoef[i] = it->second.real();
    imagCoef[i] = it->second.imag();
  }

  // Interpolate
  double real = fsScalar->interpolate(element, realCoef, xyz);
  double imag = fsScalar->interpolate(element, imagCoef, xyz);

  return Complex(real, imag);
}
