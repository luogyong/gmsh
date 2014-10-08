#ifndef _FORMULATIONJFLEETHREE_H_
#define _FORMULATIONJFLEETHREE_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "FunctionSpaceScalar.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationJFLee.h"

/**
   @class FormulationJFLeeThree
   @brief Helping class for FormulationJFLee (<grad(rho), phi> term)

   Helping class for FormulationJFLee (<grad(rho), phi> term)

   FormulationJFLee is a friend of FormulationJFLeeThree
 */

class FormulationJFLeeThree: public FormulationBlock<Complex>{
 private:
  friend class FormulationJFLee;

 private:
  // JF Lee coef 2 //
  Complex C2;

  // Function Space & Domain //
  const FunctionSpaceScalar* ffield;
  const FunctionSpaceVector* ttest;
  const GroupOfElement*      ddomain;

  // Local Terms //
  const TermGradGrad<double>* local;

 private:
  FormulationJFLeeThree(void);
  FormulationJFLeeThree(const GroupOfElement& domain,
                        const FunctionSpaceScalar& field,
                        const FunctionSpaceVector& test,
                        Complex C2,
                        const TermGradGrad<double>& local);

 public:
  virtual ~FormulationJFLeeThree(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationJFLeeThree::~FormulationJFLeeThree
   Deletes this FormulationJFLeeThree
*/

#endif
