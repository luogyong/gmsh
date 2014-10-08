#ifndef _FORMULATIONJFLEESIX_H_
#define _FORMULATIONJFLEESIX_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermCurlCurl.h"

#include "FormulationBlock.h"
#include "FormulationJFLee.h"

/**
   @class FormulationJFLeeSix
   @brief Helping class for FormulationJFLee (<curl(e), curl(phi)> term)

   Helping class for FormulationJFLee (<curl(e), curl(phi)> term)

   FormulationJFLee is a friend of FormulationJFLeeSix
 */

class FormulationJFLeeSix: public FormulationBlock<Complex>{
 private:
  friend class FormulationJFLee;

 private:
  // JF Lee coef //
  Complex C1;

  // Function Space & Domain //
  const FunctionSpaceVector* ffield;
  const FunctionSpaceVector* ttest;
  const GroupOfElement*      ddomain;

  // Local Terms //
  const TermCurlCurl<double>* local;

 private:
  FormulationJFLeeSix(void);
  FormulationJFLeeSix(const GroupOfElement& domain,
                      const FunctionSpaceVector& field,
                      const FunctionSpaceVector& test,
                      Complex C1,
                      const TermCurlCurl<double>& local);

 public:
  virtual ~FormulationJFLeeSix(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationJFLeeSix::~FormulationJFLeeSix
   Deletes this FormulationJFLeeSix
*/

#endif
