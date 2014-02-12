#include <vector>
#include "TermProjectionGrad.h"
#include "SmallFem.h"

using namespace std;

// Real Version //
// ------------ //
template<>
void TermProjectionGrad<double>::computeC(const Basis& basis,
                                          const BFunction& getFunction,
                                          const fullVector<double>& gW,
                                          fullMatrix<double>**& cM){
  const size_t nG = gW.size();
  size_t k;

  // Alloc //
  cM = new fullMatrix<double>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<double>(3 * nG, nFunction);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      (basis.*getFunction)(s);

    // Loop on Gauss Points
    k = 0;

    for(size_t g = 0; g < nG; g++){
      for(size_t a = 0; a < 3; a++){
        for(size_t i = 0; i < nFunction; i++)
          (*cM[s])(k, i) = gW(g) * phi(i, k);

        k++;
      }
    }
  }
}

template<>
fullVector<double> TermProjectionGrad<double>::
interpolateGrad(const MElement& element,
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
      throw Exception("TermProjectionGrad::interpolateGrad unknown dof %s",
                      dof[i].toString().c_str());

    coef[i] = it->second;
  }

  // Interpolate & Return //
  return fsScalar->interpolateDerivative(element, coef, xyz);
}

// Complex Version //
// --------------- //
template<>
void TermProjectionGrad<Complex>::computeC(const Basis& basis,
                                           const BFunction& getFunction,
                                           const fullVector<double>& gW,
                                           fullMatrix<Complex>**& cM){
  const size_t nG = gW.size();
  size_t k;

  // Alloc //
  cM = new fullMatrix<Complex>*[nOrientation];

  for(size_t s = 0; s < nOrientation; s++)
    cM[s] = new fullMatrix<Complex>(3 * nG, nFunction);

  // Fill //
  for(size_t s = 0; s < nOrientation; s++){
    // Get functions for this Orientation
    const fullMatrix<double>& phi =
      (basis.*getFunction)(s);

    // Loop on Gauss Points
    k = 0;

    for(size_t g = 0; g < nG; g++){
      for(size_t a = 0; a < 3; a++){
        for(size_t i = 0; i < nFunction; i++)
          (*cM[s])(k, i) = Complex(gW(g) * phi(i, k), 0);

        k++;
      }
    }
  }
}

template<>
fullVector<Complex> TermProjectionGrad<Complex>::
interpolateGrad(const MElement& element,
                const fullVector<double>& xyz) const{

  // Get Dofs associated to element //
  const vector<Dof>  dof = fsScalar->getKeys(element);
  const size_t      nDof = dof.size();

  // Get Values of these Dofs //
  map<Dof, Complex>::const_iterator end = dofValue->end();
  map<Dof, Complex>::const_iterator it;
  vector<double> reCoef(nDof);
  vector<double> imCoef(nDof);

  for(size_t i = 0; i < nDof; i++){
    it = dofValue->find(dof[i]);
    if(it == end)
      throw Exception("TermProjectionGrad::interpolateGrad unknown dof %s",
                      dof[i].toString().c_str());

    reCoef[i] = it->second.real();
    imCoef[i] = it->second.imag();
  }

  // Interpolate
  fullVector<double> re = fsScalar->interpolateDerivative(element, reCoef, xyz);
  fullVector<double> im = fsScalar->interpolateDerivative(element, imCoef, xyz);

  // Return //
  fullVector<Complex> ret(3);

  ret(0) = Complex(re(0), im(0));
  ret(1) = Complex(re(1), im(1));
  ret(2) = Complex(re(2), im(2));

  return ret;
}
