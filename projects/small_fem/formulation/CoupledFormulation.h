#ifndef _COUPLEDFORMULATION_H_
#define _COUPLEDFORMULATION_H_

#include <list>
#include "Formulation.h"

/**
   @interface CoupledFormulation
   @brief Base interface for coupled Formulations%s

   This is the base interface for coupled Formulations%s.

   A coupled Formulation is a @em container for Formulation%s
   that means something only if taken as a group.

   For instance, Formulation%s that use coupled FunctionSpace%s
   shall be embeded is a CoupledFormulation.
   For those Formulations%s,
   it is recommanded to implement them with @em private constructors and
   as a friends of the particular CoupledFormulation.

   A CoupledFormulation can give access to a list of Formulation%s.
*/

template<typename scalar>
class CoupledFormulation{
 public:
  virtual ~CoupledFormulation(void);

  virtual
    const std::list<const Formulation<scalar>*>& getFormulations(void) const = 0;
};

/**
   @fn CoupledFormulation::~CoupledFormulation
   Deletes this CoupledFormulation
   **

   @fn CoupledFormulation::getFormulations
   @return Returns a list of the Formulation%s of this CoupledFormulation
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "CoupledFormulationInclusion.h"

#endif
