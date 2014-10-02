#ifndef _DDMCONTEXTOSRC_H_
#define _DDMCONTEXTOSRC_H_

#include "DDMContext.h"

/**
   @class DDMContextOSRC
   @brief Context for DDM OSRC

   This class is a DDMContext for OSRC.
   In addition to a DDMContext, this class alows the definition
   of the OSRC wavenumber and its chi value.
 */

class DDMContextOSRC: public DDMContext{
 private:
  const std::vector<const FunctionSpaceScalar*>* phi;

  int     NPade;
  double  k;
  Complex keps;

 public:
  DDMContextOSRC(const GroupOfElement& domain,
                 const FunctionSpace& fSpace,
                 const std::vector<const FunctionSpaceScalar*>& phi,
                 double k, Complex keps, int NPade);

  virtual ~DDMContextOSRC(void);

  int     getNPade(void) const;
  double  getWavenumber(void) const;
  Complex getComplexWavenumber(void) const;

  const std::vector<const FunctionSpaceScalar*>&
    getAuxFunctionSpace(void) const;
};

/**
   @fn DDMContextOSRC::DDMContextOSRC
   @param domain The border for the OSRC problem
   @param fSpace The FunctionSpace for the OSRC problem
   @param phi A vector with the auxiliary FunctionSpace for the OSRC problem
   @param k The wavenumber of the OSRC problem
   @param keps The complexified wavenumber of the OSRC problem
   @param NPade The number of Pade term tu use

   Instanciates a new DDMContextOSRC with the given parameters
   **

   @fn DDMContextOSRC::~DDMContext
   Deletes this DDMContextOSRC
   **

   @fn DDMContextOSRC::getNPade
   @return Returns the number of Pade term of this DDMContextOSRC
   **

   @fn DDMContextOSRC::getWavenumer
   @return Returns the wavenumber of this DDMContextOSRC
   **

   @fn DDMContextOSRC::getComplexWavenumer
   @return Returns the complexified wavenumber of this DDMContextOSRC
   **

   @fn DDMContextOSRC::getAuxFunctionSpace
   @return Returns the vector of auxiliary FunctionSpace of this DDMContextOSRC
 */

//////////////////////
// Inline Functions //
//////////////////////

inline int DDMContextOSRC::getNPade(void) const{
  return NPade;
}

inline double DDMContextOSRC::getWavenumber(void) const{
  return k;
}

inline Complex DDMContextOSRC::getComplexWavenumber(void) const{
  return keps;
}

inline
const std::vector<const FunctionSpaceScalar*>&
DDMContextOSRC::getAuxFunctionSpace(void) const{
  return *phi;
}

#endif
