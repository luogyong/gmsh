#ifndef _FORMULATIONLAGRANGETWO_H_
#define _FORMULATIONLAGRANGETWO_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "FormulationBlock.h"

#include "FormulationLagrange.h"

/**
   @class FormulationLagrangeTwo
   @brief Helping class for FormulationLagrange (Lagrange is unknown)

   Helping class for FormulationLagrange (Lagrange is unknown)

   FormulationLagrangeTwo is a friend of FormulationLagrange
 */

class FormulationLagrangeTwo: public FormulationBlock<Complex>{
 private:
  friend class FormulationLagrange;

 private:
  // Function Space & Domain //
  const FunctionSpaceScalar* ffield;
  const FunctionSpaceScalar* ttest;
  const GroupOfElement*      ddomain;

  // Local Terms //
  const TermFieldField<double>* localTerm;

 private:
  FormulationLagrangeTwo(void);
  FormulationLagrangeTwo(const GroupOfElement& domain,
                         const FunctionSpaceScalar& field,
                         const FunctionSpaceScalar& test,
                         const TermFieldField<double>& localTerm);

 public:
  virtual ~FormulationLagrangeTwo(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationLagrangeTwo::~FormulationLagrangeTwo
   Deletes this FormulationLagrangeTwo
*/

#endif
