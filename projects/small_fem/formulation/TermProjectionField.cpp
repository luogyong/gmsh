#include <vector>

#include "TermProjectionField.h"
#include "SmallFem.h"

using namespace std;

// Real Version //
// ------------ //
template<>
double TermProjectionField<double>::
interpolate(const MElement& element,
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
Complex TermProjectionField<Complex>::
interpolate(const MElement& element,
            const fullVector<double>& xyz) const{

  // Get Dofs associated to element //
  vector<Dof> dof;
  fsScalar->getKeys(element, dof);

  // Get Values of these Dofs //
  const size_t nDof = dof.size();

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
