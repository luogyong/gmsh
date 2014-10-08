#ifndef _FORMULATIONJFLEEEIGHT_H_
#define _FORMULATIONJFLEEEIGHT_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "FunctionSpaceScalar.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationJFLee.h"

/**
   @class FormulationJFLeeEight
   @brief Helping class for FormulationJFLee (phi, grad(rho)> term)

   Helping class for FormulationJFLee (<phi, grad(rho)> term)

   FormulationJFLee is a friend of FormulationJFLeeEight
 */

class FormulationJFLeeEight: public FormulationBlock<Complex>{
 private:
  friend class FormulationJFLee;

 private:
  // Function Space & Domain //
  const FunctionSpaceVector* ffield;
  const FunctionSpaceScalar* ttest;
  const GroupOfElement*      ddomain;

  // Local Terms //
  const TermGradGrad<double>* local;

 private:
  FormulationJFLeeEight(void);
  FormulationJFLeeEight(const GroupOfElement& domain,
                        const FunctionSpaceVector& field,
                        const FunctionSpaceScalar& test,
                        const TermGradGrad<double>& local);

 public:
  virtual ~FormulationJFLeeEight(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationJFLeeEight::~FormulationJFLeeEight
   Deletes this FormulationJFLeeEight
*/

#endif
