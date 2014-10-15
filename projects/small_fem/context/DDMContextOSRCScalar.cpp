#include "DDMContextOSRCScalar.h"

using namespace std;

DDMContextOSRCScalar::
DDMContextOSRCScalar(const GroupOfElement& domain,
                     const FunctionSpace& fSpace,
                     const vector<const FunctionSpaceScalar*>& phi,
                     double k, Complex keps, int NPade){
  // Check if scalar //
  if(!fSpace.isScalar())
    throw Exception("DDMContextOSRCScalar: need a scalar function space");

  // Data for OSRCScalar //
  this->domain = &domain;
  this->fSpace = &fSpace;
  this->phi    = &phi;
  this->NPade  = NPade;
  this->k      = k;
  this->keps   = keps;
}

DDMContextOSRCScalar::~DDMContextOSRCScalar(void){
}
