#ifndef _FORMULATIONPROJECTIONSCALAR_H_
#define _FORMULATIONPROJECTIONSCALAR_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "fullMatrix.h"

#include "TermFieldField.h"
#include "TermProjectionField.h"

#include "Formulation.h"

/**
   @class FormulationProjectionScalar
   @brief Formulation for the L2 projection of a scalar function

   Scalar Formulation for the L2 projection problem
 */

template<typename scalar>
class FormulationProjectionScalar: public Formulation<scalar>{
 private:
  // Function Space & Basis & Domain //
  const FunctionSpaceScalar* fspace;
  const GroupOfElement*      goe;
  const Basis*               basis;

  // For real version (Local Terms) //
  TermFieldField*      localTerms1;
  TermProjectionField* localTerms2;

  // For complex version //
  Complex (*f)(fullVector<double>& xyz);
  fullMatrix<double>*   gC;
  fullVector<double>*   gW;
  GroupOfJacobian*      jac;

 public:
  FormulationProjectionScalar(const GroupOfElement& domain,
                              const FunctionSpaceScalar& fs,
                              scalar (*f)(fullVector<double>& xyz));

  virtual ~FormulationProjectionScalar(void);

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual scalar rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationProjectionScalar::FormulationProjectionScalar
   @param domain The domain of this Formulation
   @param fs A FunctionSpaceScalar for both unknown and test field
   @param f The function to project

   Instantiates a new FormulationProjectionScalar to project the given function
   **

   @fn FormulationProjectionScalar::~FormulationProjectionScalar
   Deletes the this FormulationProjectionScalar
*/

#endif
