#ifndef _DDMCONTEXTOSRCSCALAR_H_
#define _DDMCONTEXTOSRCSCALAR_H_

#include "DDMContext.h"

/**
   @class DDMContextOSRCScalar
   @brief Context for DDM scalar OSRC

   This class is a DDMContext for scalar OSRC.
   In addition to a DDMContext, this class alows the definition
   of the scalar OSRC wavenumber and its chi value.
 */

class DDMContextOSRCScalar: public DDMContext{
 private:
  const std::vector<const FunctionSpaceScalar*>* phi;

  int     NPade;
  double  k;
  Complex keps;

 public:
  DDMContextOSRCScalar(const GroupOfElement& domain,
                       const FunctionSpace& fSpace,
                       const std::vector<const FunctionSpaceScalar*>& phi,
                       double k, Complex keps, int NPade);

  virtual ~DDMContextOSRCScalar(void);

  int     getNPade(void) const;
  double  getWavenumber(void) const;
  Complex getComplexWavenumber(void) const;

  const std::vector<const FunctionSpaceScalar*>&
    getAuxFunctionSpace(void) const;
};

/**
   @fn DDMContextOSRCScalar::DDMContextOSRCScalar
   @param domain The border for the scalar OSRC problem
   @param fSpace The FunctionSpace for the scalar OSRC problem
   @param phi A vector with the auxiliary FunctionSpace for scalar OSRC problem
   @param k The wavenumber of the scalar OSRC problem
   @param keps The complexified wavenumber of the scalar OSRC problem
   @param NPade The number of Pade term tu use

   Instanciates a new DDMContextOSRCScalar with the given parameters
   **

   @fn DDMContextOSRCScalar::~DDMContext
   Deletes this DDMContextOSRCScalar
   **

   @fn DDMContextOSRCScalar::getNPade
   @return Returns the number of Pade term of this DDMContextOSRCScalar
   **

   @fn DDMContextOSRCScalar::getWavenumer
   @return Returns the wavenumber of this DDMContextOSRCScalar
   **

   @fn DDMContextOSRCScalar::getComplexWavenumer
   @return Returns the complexified wavenumber of this DDMContextOSRCScalar
   **

   @fn DDMContextOSRCScalar::getAuxFunctionSpace
   @return Returns the vector of auxiliary FunctionSpace
   of this DDMContextOSRCScalar
 */

//////////////////////
// Inline Functions //
//////////////////////

inline int DDMContextOSRCScalar::getNPade(void) const{
  return NPade;
}

inline double DDMContextOSRCScalar::getWavenumber(void) const{
  return k;
}

inline Complex DDMContextOSRCScalar::getComplexWavenumber(void) const{
  return keps;
}

inline
const std::vector<const FunctionSpaceScalar*>&
DDMContextOSRCScalar::getAuxFunctionSpace(void) const{
  return *phi;
}

#endif
