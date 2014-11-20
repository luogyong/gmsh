#ifndef _DDMCONTEXTEMDA_H_
#define _DDMCONTEXTEMDA_H_

#include "DDMContext.h"

/**
   @class DDMContextEMDA
   @brief Context for DDM EMDA

   This class is a DDMContext for EMDA.
   In addition to a DDMContext, this class alows the definition
   of the EMDA wavenumber and its chi value.
 */

class DDMContextEMDA: public DDMContext{
 private:
  double k;
  double chi;

 public:
  DDMContextEMDA(const GroupOfElement& domain,
                 std::vector<const GroupOfElement*>& dirichlet,
                 const FunctionSpace& fSpace,
                 double k, double chi);

  virtual ~DDMContextEMDA(void);

  double getWavenumber(void) const;
  double getChi(void) const;
};

/**
   @fn DDMContextEMDA::DDMContextEMDA
   @param domain The DDM border for the EMDA problem
   @param dirichlet The Dirichlet border for the EMDA problem
   @param fSpace The FunctionSpace for the EMDA problem
   @param k The wavenumber of the EMDA problem
   @param chi The chi of the EMDA problem

   Instanciates a new DDMContextEMDA with the given parameters
   **

   @fn DDMContextEMDA::~DDMContext
   Deletes this DDMContextEMDA
   **

   @fn DDMContextEMDA::getWavenumer
   @return Returns the wavenumber of this DDMContextEMDA
   **

   @fn DDMContextEMDA::getChi
   @return Returns the chi of this DDMContextEMDA
 */

//////////////////////
// Inline Functions //
//////////////////////

inline double DDMContextEMDA::getWavenumber(void) const{
  return k;
}

inline double DDMContextEMDA::getChi(void) const{
  return chi;
}

#endif
