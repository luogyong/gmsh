#ifndef _FORMULATIONPROJECTIONVECTOR_H_
#define _FORMULATIONPROJECTIONVECTOR_H_

#include <complex>
#include "FunctionSpaceVector.h"
#include "fullMatrix.h"

#include "TermGradGrad.h"
#include "TermProjectionGrad.h"

#include "Formulation.h"

/**
   @class FormulationProjectionVector
   @brief Formulation for the Projection of Vectorial Function problem

   Vectorial Formulation for the L2 projection problem
 */

template<typename scalar>
class FormulationProjectionVector: public Formulation<scalar>{
 private:
  // Function Space & Basis & Domain //
  const FunctionSpaceVector* fspace;
  const GroupOfElement*      goe;
  const Basis*               basis;

  // For real version (Local Terms) //
  TermGradGrad*       localTerms1;
  TermProjectionGrad* localTerms2;

  // For complex version //
  fullVector<std::complex<double> > (*f)(fullVector<double>& xyz);
  fullMatrix<double>*   gC;
  fullVector<double>*   gW;
  GroupOfJacobian*      jac;

 public:
  FormulationProjectionVector(const GroupOfElement& goe,
                              const FunctionSpaceVector& fs,
                              fullVector<scalar> (*f)(fullVector<double>& xyz));

  virtual ~FormulationProjectionVector(void);

  virtual bool isGeneral(void) const;

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId)  const;
  virtual scalar weakB(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual scalar rhs(size_t equationI, size_t elementId)           const;

  virtual const FunctionSpace&  fs(void)     const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationProjectionVector::FormulationProjectionVector
   @param f The function to project
   @param fs A FunctionSpaceEdge

   Instantiates a new FormulationProjectionVector to project
   the given function@n

   FormulationProjectionVector will use the given FunctionSpace
   for the projection
   **

   @fn FormulationProjectionVector::~FormulationProjectionVector
   Deletes the this FormulationProjectionVector
*/

#endif
