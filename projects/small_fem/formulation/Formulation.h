#ifndef _FORMULATION_H_
#define _FORMULATION_H_

#include <string>
#include "FunctionSpace.h"

/**
   @interface Formulation
   @brief Base interface of a finite element formulation

   This is the base interface of a finite element formulation.

   A Formulation is able to compute a weak term of a finite element matrix
   and its right hand side.

   A Formulation is defined with two FunctionSpaces:
   @li One to insterpolate the unknown field
   @li One for the test functions

   Finaly a Formulation must be defined on a geomtrical domain.
 */

template<typename scalar>
class Formulation{
 public:
  virtual ~Formulation(void);

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId) const = 0;
  virtual scalar rhs(size_t equationI, size_t elementId)          const = 0;

  virtual const FunctionSpace&  field(void)  const = 0;
  virtual const FunctionSpace&  test(void)   const = 0;
  virtual const GroupOfElement& domain(void) const = 0;
};

/**
   @fn Formulation::~Formulation
   Deletes this Formulation
   **

   @fn Formulation::weak
   @param dofI The first index of the formulation term
   @param dofJ The second index of the formulation term
   @param elementId The element ID associated with the formulation term
   @return The value of the requested formulation term
   **

   @fn Formulation::rhs
   @param equationI The ith equation of the formulation
   @param elementId The element ID associated
   with the ith equation of the formulation
   @return The value of the ith equation right hand side
   **

   @fn Formulation::field
   @return Returns the FunctionSpace used to interpolate the unknown field
   **

   @fn Formulation::test
   @return Returns the FunctionSpace used for the test functions
   **

   @fn Formulation::domain
   @return Returns the domain of definition of this Formulation
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationInclusion.h"

#endif
