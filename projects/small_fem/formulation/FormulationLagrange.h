#ifndef _FORMULATIONLAGRANGE_H_
#define _FORMULATIONLAGRANGE_H_

#include "SmallFem.h"
#include "CoupledFormulation.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"

/**
   @class FormulationLagrange
   @brief Lagrange multipliers CoupledFormulation

   Lagrange multipliers CoupledFormulation
 */

class FormulationLagrange: public CoupledFormulation<Complex>{
 private:
  // Local Terms //
  TermFieldField<double>*      local; // Lagrange . Field && Field . Lagrange
  TermProjectionField<double>* proj;  // Proj     . Lagrange

  // Formulations //
  std::list<const Formulation<Complex>*> fList;

 public:
  FormulationLagrange(const GroupOfElement& domain,
                      const FunctionSpaceScalar& field,
                      const FunctionSpaceScalar& lagrange,
                      double (*f)(fullVector<double>& xyz));

  virtual ~FormulationLagrange(void);

  virtual
    const std::list<const Formulation<Complex>*>& getFormulations(void) const;
};

/**
   @fn FormulationLagrange::FormulationLagrange
   @param domain A GroupOfElement for the domain
   @param field FunctionSpace for the computed field (non Lagrange variable)
   @param test FunctionSpace for the test functions (Lagrange variable)
   @param f The function to impose by Lagrange multipliers

   Instantiates a new FormulationLagrange
   **

   @fn FormulationLagrange::~FormulationLagrange
   Deletes this FormulationLagrange
*/

#endif
