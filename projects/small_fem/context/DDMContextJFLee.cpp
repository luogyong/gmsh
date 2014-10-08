#include <cmath>
#include "DDMContextJFLee.h"

using namespace std;

DDMContextJFLee::DDMContextJFLee(const GroupOfElement& domain,
                                 const FunctionSpace& fSpace,
                                 const FunctionSpace& fPhi,
                                 const FunctionSpace& fRho,
                                 double k, double lc){
  // Math //
  const double pi = atan(1.0) * 4;

  // Check //
  if(fSpace.isScalar())
    throw Exception("DDMContextJFLee: "
                    "Primary function space must be vectorial");

  if(fPhi.isScalar())
    throw Exception("DDMContextJFLee: "
                    "Auxiliary function space Phi must be vectorial");

  if(!fRho.isScalar())
    throw Exception("DDMContextJFLee: "
                    "Auxiliary function space Rho must be scalar");

  // Data for JFLee //
  this->domain = &domain;
  this->fSpace = &fSpace;
  this->fPhi   = static_cast<const FunctionSpaceVector*>(&fPhi);
  this->fRho   = static_cast<const FunctionSpaceScalar*>(&fRho);
  this->k      = k;

  // Coef //
  double kMax  = pi / lc;
  double delta = sqrt((kMax * kMax - k * k) / (k * k));

  C1 = Complex(+1, 0) / (Complex(+1, 0) + Complex(0, delta));
  C2 = Complex(-1, 0) * C1;
}

DDMContextJFLee::~DDMContextJFLee(void){
}
