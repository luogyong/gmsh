#ifndef _FORMULATIONLAGRANGEONE_H_
#define _FORMULATIONLAGRANGEONE_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "Formulation.h"

#include "FormulationLagrange.h"

/**
   @class FormulationLagrangeOne
   @brief Helping class for FormulationLagrange (field is unknown)

   Helping class for FormulationLagrange (field is unknown)

   FormulationLagrangeOne is a friend of FormulationLagrange
*/

class FormulationLagrangeOne: public Formulation<Complex>{
 private:
  friend class FormulationLagrange;

 private:
  // Function Space & Domain //
  const FunctionSpaceScalar* ffield;
  const FunctionSpaceScalar* ttest;
  const GroupOfElement*      ddomain;

  // Local Terms //
  const TermFieldField<double>*      localTerm;
  const TermProjectionField<double>* projectionTerm;

 private:
  FormulationLagrangeOne(void);
  FormulationLagrangeOne(const GroupOfElement& domain,
                         const FunctionSpaceScalar& field,
                         const FunctionSpaceScalar& test,
                         const TermFieldField<double>& localTerm,
                         const TermProjectionField<double>& projectionTerm);

 public:
  virtual ~FormulationLagrangeOne(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationLagrangeOne::~FormulationLagrangeOne
   Deletes this FormulationLagrangeOne
*/

#endif
