#ifndef _FORMULATIONLAGRANGEFIELD_H_
#define _FORMULATIONLAGRANGEFIELD_H_

#include <complex>
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "Formulation.h"

/**
   @class FormulationLagrangeField
   @brief Lagrange multipliers Formulation (Tested by Field)

   Lagrange multipliers Formulation (Tested by Field)
 */

class FormulationLagrangeField: public Formulation<std::complex<double> >{
 private:
  // Function Space & Domain //
  const FunctionSpaceScalar* fsF;
  const FunctionSpaceScalar* fsT;
  const GroupOfElement*      goe;

  // Local Terms //
  TermFieldField* localTerms;

 public:
  FormulationLagrangeField(const GroupOfElement& goe,
                           const FunctionSpaceScalar& fsField,
                           const FunctionSpaceScalar& fsTest);

  virtual ~FormulationLagrangeField(void);

  virtual bool isGeneral(void) const;

  virtual std::complex<double>
    weak(size_t dofI, size_t dofJ, size_t elementId)  const;

  virtual std::complex<double>
    weakB(size_t dofI, size_t dofJ, size_t elementId) const;

  virtual std::complex<double>
    rhs(size_t equationI, size_t elementId)           const;

  virtual const FunctionSpace&  fsField(void) const;
  virtual const FunctionSpace&  fsTest(void)  const;
  virtual const GroupOfElement& domain(void)  const;
};

/**
   @fn FormulationLagrangeField::FormulationLagrangeField
   @param goe A GroupOfElement
   @param fsField FunctionSpace for the computed field
   @param fsTest FunctionSpace for the test functions

   Instantiates a new FormulationLagrangeField
   **

   @fn FormulationLagrangeField::~FormulationLagrangeField
   Deletes this FormulationLagrangeField
*/

#endif
