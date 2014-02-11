#ifndef _FORMULATIONFIELDLAGRANGE_H_
#define _FORMULATIONFIELDLAGRANGE_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "TermProjectionField.h"
#include "Formulation.h"

/**
   @class FormulationFieldLagrange
   @brief Lagrange multipliers Formulation (Tested by Lagrange)

   Lagrange multipliers Formulation (Tested by Lagrange)
 */

class FormulationFieldLagrange: public Formulation<Complex>{
 private:
  // Function Space & Domain //
  const FunctionSpaceScalar* fsF;
  const FunctionSpaceScalar* fsT;
  const GroupOfElement*      goe;

  // Local Terms //
  TermFieldField*      localTerms;
  TermProjectionField* projectionTerms;

 public:
  FormulationFieldLagrange(const GroupOfElement& domain,
                           const FunctionSpaceScalar& field,
                           const FunctionSpaceScalar& test,
                           double (*f)(fullVector<double>& xyz));

  virtual ~FormulationFieldLagrange(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationFieldLagrange::FormulationFieldLagrange
   @param domain A GroupOfElement for the domain
   @param field A FunctionSpace for the computed field
   @param test A FunctionSpace for the test functions
   @param f A function to impose by Lagrange multipliers

   Instantiates a new FormulationFieldLagrange
   **

   @fn FormulationFieldLagrange::~FormulationFieldLagrange
   Deletes this FormulationFieldLagrange
*/

#endif
