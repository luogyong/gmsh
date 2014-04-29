#ifndef _FORMULATION_H_
#define _FORMULATION_H_

/**
   @interface Formulation
   @brief Base interface of a finite element formulation

   Base interface of a finite element formulation.

   A Formulation can be eather a FormulationBlock or a FormulationCoupled.
 */

template<typename scalar>
class Formulation{
 public:
  virtual ~Formulation(void);

  virtual bool isBlock(void) const = 0;
};

/**
   @fn Formulation::~Formulation
   Deletes this Formulation
   **

   @fn Formulation::isBlock
   @return
   Returns true if this Formulation is a FormulationBlock
   and false if this Formulation is a FormulationCoupled.
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationInclusion.h"

#endif
