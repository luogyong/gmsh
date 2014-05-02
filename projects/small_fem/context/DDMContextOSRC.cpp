#include "DDMContextOSRC.h"

using namespace std;

DDMContextOSRC::DDMContextOSRC(const GroupOfElement& domain,
                               const FunctionSpaceScalar& fSpace,
                               const vector<const FunctionSpaceScalar*>& phi,
                               double k, Complex keps, int NPade){
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
