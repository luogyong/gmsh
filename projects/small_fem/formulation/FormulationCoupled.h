#ifndef _FORMULATIONCOUPLED_H_
#define _FORMULATIONCOUPLED_H_

#include <list>
#include "FormulationBlock.h"

/**
   @interface FormulationCoupled
   @brief Base interface for coupled Formulations%s

   This is the base interface for coupled Formulations%s.

   A coupled Formulation is a @em container for Formulation%s
   that means something only if taken as a group.

   For instance, Formulation%s that use coupled FunctionSpace%s
   shall be embeded is a FormulationCoupled.
   For those Formulations%s,
   it is recommanded to implement them with @em private constructors and
   as a friends of the particular FormulationCoupled.

   A FormulationCoupled can give access to a list of Formulation%s.
*/

template<typename scalar>
class FormulationCoupled: public Formulation<scalar>{
 public:
  virtual ~FormulationCoupled(void);

  virtual
    const std::list<const FormulationBlock<scalar>*>&
                                           getFormulationBlocks(void) const = 0;
};

/**
   @fn FormulationCoupled::~FormulationCoupled
   Deletes this FormulationCoupled
   **

   @fn FormulationCoupled::getFormulationBlocks
   @return Returns a list of the FormulationBlock%s of this FormulationCoupled
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationCoupledInclusion.h"

#endif
