#include "DDMContext.h"

using namespace std;

DDMContext::DDMContext(void){
  // Init to zero //
  system     = NULL;
  domain     = NULL;
  fSpace     = NULL;
  fSpaceG    = NULL;
  ddm        = NULL;

  dirichlet.clear();
}

DDMContext::~DDMContext(void){
}
