#ifndef _FORMULATIONPROJECTIONVECTOR_H_
#define _FORMULATIONPROJECTIONVECTOR_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "fullMatrix.h"

#include "TermGradGrad.h"
#include "TermProjectionGrad.h"

#include "Formulation.h"

/**
   @class FormulationProjectionVector
   @brief Formulation for the L2 projection of a vectorial function

   Vectorial Formulation for the L2 projection problem
 */

template<typename scalar>
class FormulationProjectionVector: public Formulation<scalar>{
 private:
  // Function Space & Basis & Domain //
  const FunctionSpaceVector* fspace;
  const GroupOfElement*      ddomain;

  // For real version (Local Terms) //
  TermGradGrad*               localTerms1;
  TermProjectionGrad<scalar>* localTerms2;

 public:
  FormulationProjectionVector(const GroupOfElement& domain,
                              const FunctionSpaceVector& fs,
                              fullVector<scalar> (*f)(fullVector<double>& xyz));

  virtual ~FormulationProjectionVector(void);

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual scalar rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationProjectionVector::FormulationProjectionVector
   @param domain The domain of this Formulation
   @param fs A FunctionSpaceVector  for both unknown and test field
   @param f The function to project

   Instantiates a new FormulationProjectionScalar to project the given function
   **

   @fn FormulationProjectionVector::~FormulationProjectionVector
   Deletes the this FormulationProjectionVector
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationProjectionVectorInclusion.h"

#endif
