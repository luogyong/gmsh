#ifndef _DDMCONTEXTOO2_H_
#define _DDMCONTEXTOO2_H_

#include "DDMContext.h"

/**
   @class DDMContextOO2
   @brief Context for DDM OO2

   This class is a DDMContext for OO2.
   In addition to a DDMContext, this class alows the definition
   of the OO2 A and B values.
 */

class DDMContextOO2: public DDMContext{
 private:
  Complex a;
  Complex b;

 public:
  DDMContextOO2(const GroupOfElement& domain,
                std::vector<const GroupOfElement*>& dirichlet,
                const FunctionSpace& fSpace,
                const FunctionSpace& fSpaceG,
                Complex a, Complex b);

  virtual ~DDMContextOO2(void);

  Complex getA(void) const;
  Complex getB(void) const;
};

/**
   @fn DDMContextOO2::DDMContextOO2
   @param domain The DDM border for the OO2 problem
   @param dirichlet The Dirichlet border for the OO2 problem
   @param fSpace The FunctionSpace for the OO2 field
   @param fSpaceG The FunctionSpace for the OO2 DDM
   @param a The A value of the OO2 problem
   @param b The B value of the OO2 problem

   Instanciates a new DDMContextOO2 with the given parameters
   **

   @fn DDMContextOO2::~DDMContext
   Deletes this DDMContextOO2
   **

   @fn DDMContextOO2::getA
   @return Returns the A value of this DDMContextOO2
   **

   @fn DDMContextOO2::getB
   @return Returns the B value of this DDMContextOO2
 */

//////////////////////
// Inline Functions //
//////////////////////

inline Complex DDMContextOO2::getA(void) const{
  return a;
}

inline Complex DDMContextOO2::getB(void) const{
  return b;
}

#endif
