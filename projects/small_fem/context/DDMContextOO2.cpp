#include "DDMContextOO2.h"

using namespace std;

DDMContextOO2::DDMContextOO2(const GroupOfElement& domain,
                             const FunctionSpaceScalar& fSpace,
                             Complex a, Complex b){
  // Data for OO2 //
  this->domain = &domain;
  this->fSpace = &fSpace;
  this->a      = a;
  this->b      = b;
}

DDMContextOO2::~DDMContextOO2(void){
}
