#ifndef _DDMCONTEXTOSRCVECTOR_H_
#define _DDMCONTEXTOSRCVECTOR_H_

#include "DDMContext.h"

/**
   @class DDMContextOSRCVector
   @brief Context for DDM vectorial OSRC

   This class is a DDMContext for vectorial OSRC.
   In addition to a DDMContext, this class alows the definition
   of the vectorial OSRC wavenumber and its chi value.
 */

class DDMContextOSRCVector: public DDMContext{
 private:
  const std::vector<const FunctionSpaceVector*>* phi;
  const std::vector<const FunctionSpaceScalar*>* rho;
  const FunctionSpaceVector*                     r;

  int     NPade;
  double  k;
  Complex keps;

 public:
  DDMContextOSRCVector(const GroupOfElement& domain,
                       const FunctionSpace& fSpace,
                       const std::vector<const FunctionSpaceVector*>& phi,
                       const std::vector<const FunctionSpaceScalar*>& rho,
                       const FunctionSpaceVector& r,
                       double k, Complex keps, int NPade);

  virtual ~DDMContextOSRCVector(void);

  int     getNPade(void) const;
  double  getWavenumber(void) const;
  Complex getComplexWavenumber(void) const;

  const std::vector<const FunctionSpaceVector*>&
    getPhiFunctionSpace(void) const;

  const std::vector<const FunctionSpaceScalar*>&
    getRhoFunctionSpace(void) const;

  const FunctionSpaceVector& getRFunctionSpace(void) const;
};

/**
   @fn DDMContextOSRCVector::DDMContextOSRCVector
   @param domain The border for the vectorial OSRC problem
   @param fSpace The FunctionSpace for the vectorial OSRC problem
   @param phi A vector with the auxiliary FunctionSpace for vector OSRC problem
   @param k The wavenumber of the vectorial OSRC problem
   @param keps The complexified wavenumber of the vectorial OSRC problem
   @param NPade The number of Pade term tu use

   Instanciates a new DDMContextOSRCVector with the given parameters
   **

   @fn DDMContextOSRCVector::~DDMContext
   Deletes this DDMContextOSRCVector
   **

   @fn DDMContextOSRCVector::getNPade
   @return Returns the number of Pade term of this DDMContextOSRCVector
   **

   @fn DDMContextOSRCVector::getWavenumer
   @return Returns the wavenumber of this DDMContextOSRCVector
   **

   @fn DDMContextOSRCVector::getComplexWavenumer
   @return Returns the complexified wavenumber of this DDMContextOSRCVector
   **

   @fn DDMContextOSRCVector::getPhiFunctionSpace
   @return Returns the vector of auxiliary FunctionSpaceVector
   of this DDMContextOSRCVector
   **

   @fn DDMContextOSRCVector::getRhoFunctionSpace
   @return Returns the vector of auxiliary FunctionSpaceScalar
   of this DDMContextOSRCVector
   **

   @fn DDMContextOSRCVector::getRFunctionSpace
   @return Returns the auxiliary FunctionSpace of this DDMContextOSRCVector
 */

//////////////////////
// Inline Functions //
//////////////////////

inline int DDMContextOSRCVector::getNPade(void) const{
  return NPade;
}

inline double DDMContextOSRCVector::getWavenumber(void) const{
  return k;
}

inline Complex DDMContextOSRCVector::getComplexWavenumber(void) const{
  return keps;
}

inline const std::vector<const FunctionSpaceVector*>&
DDMContextOSRCVector::getPhiFunctionSpace(void) const{
  return *phi;
}

inline const std::vector<const FunctionSpaceScalar*>&
DDMContextOSRCVector::getRhoFunctionSpace(void) const{
  return *rho;
}

inline const FunctionSpaceVector&
DDMContextOSRCVector::getRFunctionSpace(void) const{
  return *r;
}

#endif
