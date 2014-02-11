#ifndef _FORMULATIONFIELDLAGRANGE_H_
#define _FORMULATIONFIELDLAGRANGE_H_

#include <complex>
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "TermProjectionField.h"
#include "Formulation.h"

/**
   @class FormulationFieldLagrange
   @brief Lagrange multipliers Formulation (Tested by Lagrange)

   Lagrange multipliers Formulation (Tested by Lagrange)
 */

class FormulationFieldLagrange: public Formulation<std::complex<double> >{
 private:
  // Function Space & Domain //
  const FunctionSpaceScalar* fsF;
  const FunctionSpaceScalar* fsT;
  const GroupOfElement*      goe;

  // Local Terms //
  TermFieldField*      localTerms;
  TermProjectionField* projectionTerms;

 public:
  FormulationFieldLagrange(const GroupOfElement& goe,
                           const FunctionSpaceScalar& fsField,
                           const FunctionSpaceScalar& fsTest,
                           double (*f)(fullVector<double>& xyz));

  virtual ~FormulationFieldLagrange(void);

  virtual std::complex<double>
    weak(size_t dofI, size_t dofJ, size_t elementId) const;

  virtual std::complex<double>
    rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationFieldLagrange::FormulationFieldLagrange
   @param goe A GroupOfElement
   @param fsField FunctionSpace for the computed field
   @param fsTest FunctionSpace for the test functions
   @param f Function imposed by Lagrange multipliers

   Instantiates a new FormulationFieldLagrange
   **

   @fn FormulationFieldLagrange::~FormulationFieldLagrange
   Deletes this FormulationFieldLagrange
*/

#endif
