#ifndef _FORMULATIONJFLEETWO_H_
#define _FORMULATIONJFLEETWO_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationJFLee.h"

/**
   @class FormulationJFLeeTwo
   @brief Helping class for FormulationJFLee (<phi, e> term)

   Helping class for FormulationJFLee (<phi, e> term)

   FormulationJFLee is a friend of FormulationJFLeeTwo
 */

class FormulationJFLeeTwo: public FormulationBlock<Complex>{
 private:
  friend class FormulationJFLee;

 private:
  // Wavenumber //
  double k;

  // Function Space & Domain //
  const FunctionSpaceVector* ffield;
  const FunctionSpaceVector* ttest;
  const GroupOfElement*      ddomain;

  // Local Terms //
  const TermGradGrad<double>* local;

 private:
  FormulationJFLeeTwo(void);
  FormulationJFLeeTwo(const GroupOfElement& domain,
                      const FunctionSpaceVector& field,
                      const FunctionSpaceVector& test,
                      double k,
                      const TermGradGrad<double>& local);

 public:
  virtual ~FormulationJFLeeTwo(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationJFLeeTwo::~FormulationJFLeeTwo
   Deletes this FormulationJFLeeTwo
*/

#endif
