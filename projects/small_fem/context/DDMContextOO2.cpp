#include "DDMContextOO2.h"

using namespace std;

DDMContextOO2::DDMContextOO2(const GroupOfElement& domain,
                             const FunctionSpace& fSpace,
                             Complex a, Complex b){
  // Check if scalar //
  if(!fSpace.isScalar())
    throw Exception("DDMContextOO2: need a scalar function space");

  // Data for OO2 //
  this->domain = &domain;
  this->fSpace = &fSpace;
  this->a      = a;
  this->b      = b;
}

DDMContextOO2::~DDMContextOO2(void){
}
