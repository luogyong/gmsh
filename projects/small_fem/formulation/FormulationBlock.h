#ifndef _FORMULATIONBLOCK_H_
#define _FORMULATIONBLOCK_H_

#include <string>
#include "FunctionSpace.h"
#include "Formulation.h"

/**
   @interface FormulationBlock
   @brief Base interface of a finite element formulation block

   This is the base interface of a finite element formulation block.

   A FormulationBlock is able to compute a weak term of a finite element matrix
   and its right hand side.

   A FormulationBlock is defined with two FunctionSpaces:
   @li One to insterpolate the unknown field
   @li One for the test functions

   Finaly a FormulationBlock must be defined on a geomtrical domain.
 */

template<typename scalar>
class FormulationBlock: public Formulation<scalar>{
 public:
  virtual ~FormulationBlock(void);

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId) const = 0;
  virtual scalar rhs(size_t equationI, size_t elementId)          const = 0;

  virtual const FunctionSpace&  field(void)  const = 0;
  virtual const FunctionSpace&  test(void)   const = 0;
  virtual const GroupOfElement& domain(void) const = 0;
};

/**
   @fn FormulationBlock::~FormulationBlock
   Deletes this FormulationBlock
   **

   @fn FormulationBlock::weak
   @param dofI The first index of this formulation block term
   @param dofJ The second index of this formulation block term
   @param elementId The element ID associated with this formulation block term
   @return The value of the requested formulation block term
   **

   @fn FormulationBlock::rhs
   @param equationI The ith equation of this formulation block
   @param elementId The element ID associated
   with the ith equation of this formulation block
   @return The value of the ith equation right hand side
   **

   @fn FormulationBlock::field
   @return Returns the FunctionSpace used to interpolate the unknown field
   **

   @fn FormulationBlock::test
   @return Returns the FunctionSpace used for the test functions
   **

   @fn FormulationBlock::domain
   @return Returns the domain of definition of this FormulationBlock
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationBlockInclusion.h"

#endif
