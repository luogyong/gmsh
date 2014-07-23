#include <vector>
#include "TermProjectionGrad.h"
#include "SmallFem.h"

using namespace std;

// Real Version //
// ------------ //
template<>
fullVector<double> TermProjectionGrad<double>::
interpolateGrad(const MElement& element,
                const fullVector<double>& xyz) const{

  // Get Dofs associated to element //
  vector<Dof> dof;
  fsScalar->getKeys(element, dof);

  // Get Values of these Dofs //
  const size_t nDof = dof.size();

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
fullVector<Complex> TermProjectionGrad<Complex>::
interpolateGrad(const MElement& element,
                const fullVector<double>& xyz) const{

  // Get Dofs associated to element //
  vector<Dof> dof;
  fsScalar->getKeys(element, dof);

  // Get Values of these Dofs //
  const size_t nDof = dof.size();

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
