#ifndef _FORMULATIONSTIFFNESS_H_
#define _FORMULATIONSTIFFNESS_H_

#include "FunctionSpace.h"
#include "FormulationBlock.h"
#include "Term.h"

/**
   @class FormulationStiffness
   @brief Formulation for Stiffness terms

   Formulation for Stiffness terms
 */

template<typename scalar>
class FormulationStiffness: public FormulationBlock<scalar>{
 private:
  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  Term<scalar>* localTerms;

 public:
  FormulationStiffness(const GroupOfElement& domain,
                       const FunctionSpace& field,
                       const FunctionSpace& test);

  FormulationStiffness(const GroupOfElement& domain,
                       const FunctionSpace& field,
                       const FunctionSpace& test,
                       void (*m)(fullVector<double>&, fullMatrix<scalar>&));

  virtual ~FormulationStiffness(void);

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual scalar rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationStiffness::FormulationStiffness
   @param domain A GroupOfElement
   @param field A FunctionSpace
   @param test A FunctionSpace

   Instantiates a new FormulationStiffness with given parametres:
   @li domain for the domain of this Formulation
   @li field for the function space used for the unknown field
   @li test for the function space used for the test functions
   **

   @fn FormulationStiffness::~FormulationStiffness
   Deletes this FormulationStiffness
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationStiffnessInclusion.h"

#endif
