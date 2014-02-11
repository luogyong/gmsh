#ifndef _FORMULATIONLAGRANGEFIELD_H_
#define _FORMULATIONLAGRANGEFIELD_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "Formulation.h"

/**
   @class FormulationLagrangeField
   @brief Lagrange multipliers Formulation (Tested by Field)

   Lagrange multipliers Formulation (Tested by Field)
 */

class FormulationLagrangeField: public Formulation<Complex>{
 private:
  // Function Space & Domain //
  const FunctionSpaceScalar* fsF;
  const FunctionSpaceScalar* fsT;
  const GroupOfElement*      goe;

  // Local Terms //
  TermFieldField* localTerms;

 public:
  FormulationLagrangeField(const GroupOfElement& domain,
                           const FunctionSpaceScalar& field,
                           const FunctionSpaceScalar& test);

  virtual ~FormulationLagrangeField(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationLagrangeField::FormulationLagrangeField
   @param domain A GroupOfElement for the domain
   @param field FunctionSpace for the computed field
   @param test FunctionSpace for the test functions

   Instantiates a new FormulationLagrangeField
   **

   @fn FormulationLagrangeField::~FormulationLagrangeField
   Deletes this FormulationLagrangeField
*/

#endif
