#include "DDMContext.h"

using namespace std;

const string DDMContext::typeNULL = string("null");
const string DDMContext::typeEMDA = string("emda");
const string DDMContext::typeOO2  = string("oo2");
const string DDMContext::typeOSRC = string("osrc");


DDMContext::DDMContext(void){
  // Init to zero //
  typeDDM    = typeNULL;
  system     = NULL;
  domain     = NULL;
  fSpace     = NULL;
  phi        = NULL;
  ddm        = NULL;
  k          = 0;
  EMDA_Chi   = 0;
  OO2_A      = Complex(0, 0);
  OO2_B      = Complex(0, 0);
  OSRC_keps  = Complex(0, 0);
  OSRC_NPade = 0;
}

DDMContext::~DDMContext(void){
}

void DDMContext::setToEMDA(const GroupOfElement& domain,
                           const FunctionSpaceScalar& fSpace,
                           double k, double chi){
  // DDM Type //
  typeDDM = typeEMDA;

  // Data for EMDA //
  this->domain   = &domain;
  this->fSpace   = &fSpace;
  this->k        = k;
  this->EMDA_Chi = chi;
}

void DDMContext::setToOO2(const GroupOfElement& domain,
                          const FunctionSpaceScalar& fSpace,
                          Complex a, Complex b){
  // DDM Type //
  typeDDM = typeOO2;

  // Data for OO2 //
  this->domain = &domain;
  this->fSpace = &fSpace;
  this->OO2_A  = a;
  this->OO2_B  = b;
}

void DDMContext::setToOSRC(const GroupOfElement& domain,
                           const FunctionSpaceScalar& fSpace,
                           const vector<const FunctionSpaceScalar*>& phi,
                           double k, Complex keps, int NPade){
  // DDM Type //
  typeDDM = typeOSRC;

  // Data for OSRC //
  this->domain     = &domain;
  this->fSpace     = &fSpace;
  this->phi        = &phi;
  this->k          = k;
  this->OSRC_keps  = keps;
  this->OSRC_NPade = NPade;
}
