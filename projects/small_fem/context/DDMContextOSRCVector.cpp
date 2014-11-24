#include "DDMContextOSRCVector.h"

using namespace std;

DDMContextOSRCVector::
DDMContextOSRCVector(const GroupOfElement& domain,
                     vector<const GroupOfElement*>& dirichlet,
                     const FunctionSpace& fSpace,
                     const FunctionSpace& fSpaceG,
                     const vector<const FunctionSpaceVector*>& phi,
                     const vector<const FunctionSpaceScalar*>& rho,
                     const FunctionSpaceVector& r,
                     double k, Complex keps,
                     int NPade, double theta){
  // Check if vector //
  if(fSpace.isScalar())
    throw Exception("DDMContextOSRCVector: need a vector function space");

  // Data for OSRCVector //
  this->domain    = &domain;
  this->fSpace    = &fSpace;
  this->fSpaceG   = &fSpaceG;
  this->dirichlet = dirichlet;
  this->phi       = &phi;
  this->rho       = &rho;
  this->r         = &r;
  this->NPade     = NPade;
  this->theta     = theta;
  this->k         = k;
  this->keps      = keps;
}

DDMContextOSRCVector::~DDMContextOSRCVector(void){
}
