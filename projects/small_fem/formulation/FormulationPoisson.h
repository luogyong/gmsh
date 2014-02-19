#ifndef _FORMULATIONPOISSON_H_
#define _FORMULATIONPOISSON_H_

#include "FunctionSpaceScalar.h"
#include "fullMatrix.h"

#include "TermGradGrad.h"
#include "TermProjectionField.h"

#include "Formulation.h"

/**
   @class FormulationPoisson
   @brief Formulation for the Poisson problem

   Formulation for the Poisson problem
 */

class FormulationPoisson: public Formulation<double>{
 private:
  // Function Space & Domain //
  const FunctionSpaceScalar* fspace;
  const GroupOfElement*      goe;

  // Local Terms //
  TermGradGrad<double>*        localTermsL;
  TermProjectionField<double>* localTermsR;

  // Source Term //
  double (*fSource)(fullVector<double>& xyz);

 public:
  FormulationPoisson(const GroupOfElement& domain,
                     const FunctionSpaceScalar& fs,
                     double (*f)(fullVector<double>& xyz));

  virtual ~FormulationPoisson(void);

  virtual double weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual double rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationPoisson::FormulationPoisson
   @param domain A GroupOfElement for the domain
   @param fs A FunctionSpace used for both test and unknown field
   @param f A function for the source term

   Instantiates a new FormulationPoisson
   **

   @fn FormulationPoisson::~FormulationPoisson
   Deletes this FormulationPoisson
*/

#endif
