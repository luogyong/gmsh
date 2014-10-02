#include "DDMContextEMDA.h"

using namespace std;

DDMContextEMDA::DDMContextEMDA(const GroupOfElement& domain,
                               const FunctionSpace& fSpace,
                               double k, double chi){
  // Data for EMDA //
  this->domain = &domain;
  this->fSpace = &fSpace;
  this->k      = k;
  this->chi    = chi;
}

DDMContextEMDA::~DDMContextEMDA(void){
}
