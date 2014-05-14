#ifndef _FORMULATIONMASS_H_
#define _FORMULATIONMASS_H_

#include "FunctionSpace.h"
#include "FormulationBlock.h"
#include "Term.h"

/**
   @class FormulationMass
   @brief Formulation for Mass terms

   Formulation for Mass terms
 */

template<typename scalar>
class FormulationMass: public FormulationBlock<scalar>{
 private:
  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  Term<scalar>* localTerms;

 public:
  FormulationMass(const GroupOfElement& domain,
                  const FunctionSpace& field,
                  const FunctionSpace& test);

  FormulationMass(const GroupOfElement& domain,
                  const FunctionSpace& field,
                  const FunctionSpace& test,
                  void (*m)(fullVector<double>&, fullMatrix<scalar>&));

  FormulationMass(const GroupOfElement& domain,
                  const FunctionSpace& field,
                  const FunctionSpace& test,
                  scalar (*m)(fullVector<double>&));

  virtual ~FormulationMass(void);

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual scalar rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationMass::FormulationMass
   @param domain A GroupOfElement
   @param field A FunctionSpace
   @param test A FunctionSpace

   Instantiates a new FormulationMass with given parametres:
   @li domain for the domain of this Formulation
   @li field for the function space used for the unknown field
   @li test for the function space used for the test functions
   **

   @fn FormulationMass::~FormulationMass
   Deletes this FormulationMass
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationMassInclusion.h"

#endif
