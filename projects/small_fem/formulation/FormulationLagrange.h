#ifndef _FORMULATIONLAGRANGE_H_
#define _FORMULATIONLAGRANGE_H_

#include "SmallFem.h"
#include "FormulationCoupled.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"

/**
   @class FormulationLagrange
   @brief Lagrange multipliers FormulationCoupled

   Lagrange multipliers FormulationCoupled
 */

class FormulationLagrange: public FormulationCoupled<Complex>{
 private:
  // Local Terms //
  TermFieldField<double>*      local; // Lagrange . Field && Field . Lagrange
  TermProjectionField<double>* proj;  // Proj     . Lagrange

  // Formulations //
  std::list<const FormulationBlock<Complex>*> fList;

 public:
  FormulationLagrange(const GroupOfElement& domain,
                      const FunctionSpaceScalar& field,
                      const FunctionSpaceScalar& lagrange,
                      double (*f)(fullVector<double>& xyz));

  virtual ~FormulationLagrange(void);

  virtual
    const std::list<const FormulationBlock<Complex>*>&
                                               getFormulationBlocks(void) const;
  virtual bool isBlock(void) const;
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
