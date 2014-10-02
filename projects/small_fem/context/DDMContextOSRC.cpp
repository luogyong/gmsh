#include "DDMContextOSRC.h"

using namespace std;

DDMContextOSRC::DDMContextOSRC(const GroupOfElement& domain,
                               const FunctionSpace& fSpace,
                               const vector<const FunctionSpaceScalar*>& phi,
                               double k, Complex keps, int NPade){
  // Check if scalar //
  if(!fSpace.isScalar())
    throw Exception("DDMContextOSRC: need a scalar function space");

  // Data for OSRC //
  this->domain = &domain;
  this->fSpace = &fSpace;
  this->phi    = &phi;
  this->NPade  = NPade;
  this->k      = k;
  this->keps   = keps;
}

DDMContextOSRC::~DDMContextOSRC(void){
}
