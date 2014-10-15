#include "DDMContextOSRCVector.h"

using namespace std;

DDMContextOSRCVector::
DDMContextOSRCVector(const GroupOfElement& domain,
                     const FunctionSpace& fSpace,
                     const vector<const FunctionSpaceVector*>& phi,
                     const vector<const FunctionSpaceScalar*>& rho,
                     const FunctionSpaceVector& r,
                     double k, Complex keps, int NPade){
  // Check if vector //
  if(!fSpace.isVector())
    throw Exception("DDMContextOSRCVector: need a vector function space");

  // Data for OSRCVector //
  this->domain = &domain;
  this->fSpace = &fSpace;
  this->phi    = &phi;
  this->rho    = &rho;
  this->r      = &r;
  this->NPade  = NPade;
  this->k      = k;
  this->keps   = keps;
}

DDMContextOSRCVector::~DDMContextOSRCVector(void){
}
