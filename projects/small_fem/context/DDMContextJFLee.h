#ifndef _DDMCONTEXTJFLEE_H_
#define _DDMCONTEXTJFLEE_H_

#include "DDMContext.h"

#include "FunctionSpaceVector.h"
#include "FunctionSpaceScalar.h"

/**
   @class DDMContextJFLee
   @brief Context for DDM JFLee

   This class is a DDMContext for JFLee.
   In addition to a DDMContext, this class alows the definition
   of the JFLee coefficients.
 */

class DDMContextJFLee: public DDMContext{
 private:
  const FunctionSpaceVector* fPhi;
  const FunctionSpaceScalar* fRho;

  double  k;
  Complex C1;
  Complex C2;

 public:
  DDMContextJFLee(const GroupOfElement& domain,
                  std::vector<const GroupOfElement*>& dirichlet,
                  const FunctionSpace& fSpace,
                  const FunctionSpace& fPhi,
                  const FunctionSpace& fRho,
                  double k, double lc);

  virtual ~DDMContextJFLee(void);

  const FunctionSpaceVector& getPhi(void) const;
  const FunctionSpaceScalar& getRho(void) const;

  double  getK(void)  const;
  Complex getC1(void) const;
  Complex getC2(void) const;
};

/**
   @fn DDMContextJFLee::DDMContextJFLee
   @param domain The DDM border for the JFLee problem
   @param dirichlet The Dirichlet border for the JFLee problem
   @param fSpace The primary FunctionSpaceVector for the JFLee problem
   @param fPhi The auxiliary FunctionSpaceVector for the JFLee problem
   @param fRho The auxiliary FunctionSpaceScalar for the JFLee problem
   @param k The wavenumber of the JFLee problem
   @param lc The caracteristic mesh size of the JFLee problem

   Instanciates a new DDMContextJFLee with the given parameters
   **

   @fn DDMContextJFLee::~DDMContext
   Deletes this DDMContextJFLee
   **

   @fn DDMContextJFLee::getPhi
   @return Returns the vectorial auxiliary FunctionSpace of this DDMContextJFLee
   **

   @fn DDMContextJFLee::getRho
   @return Returns the scalar auxiliary FunctionSpace of this DDMContextJFLee
   **

   @fn DDMContextJFLee::getK
   @return Returns the wavenumber of this DDMContextJFLee
   **

   @fn DDMContextJFLee::getC1
   @return Returns the first coefficient of this DDMContextJFLee
   **

   @fn DDMContextJFLee::getC2
   @return Returns the second coefficient value of this DDMContextJFLee
 */

//////////////////////
// Inline Functions //
//////////////////////

inline const FunctionSpaceVector& DDMContextJFLee::getPhi(void) const{
  return *fPhi;
}

inline const FunctionSpaceScalar& DDMContextJFLee::getRho(void) const{
  return *fRho;
}

inline double DDMContextJFLee::getK(void) const{
  return k;
}

inline Complex DDMContextJFLee::getC1(void) const{
  return C1;
}

inline Complex DDMContextJFLee::getC2(void) const{
  return C2;
}

#endif
