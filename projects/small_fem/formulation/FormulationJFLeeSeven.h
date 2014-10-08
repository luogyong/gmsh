#ifndef _FORMULATIONJFLEESEVEN_H_
#define _FORMULATIONJFLEESEVEN_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"

#include "FormulationBlock.h"
#include "FormulationJFLee.h"

/**
   @class FormulationJFLeeSeven
   @brief Helping class for FormulationJFLee (rho, rho> term)

   Helping class for FormulationJFLee (<rho, rho> term)

   FormulationJFLee is a friend of FormulationJFLeeSeven
 */

class FormulationJFLeeSeven: public FormulationBlock<Complex>{
 private:
  friend class FormulationJFLee;

 private:
  // Function Space & Domain //
  const FunctionSpaceScalar* ffield;
  const GroupOfElement*      ddomain;

  // Local Terms //
  const TermFieldField<double>* local;

 private:
  FormulationJFLeeSeven(void);
  FormulationJFLeeSeven(const GroupOfElement& domain,
                        const FunctionSpaceScalar& field,
                        const TermFieldField<double>& local);

 public:
  virtual ~FormulationJFLeeSeven(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationJFLeeSeven::~FormulationJFLeeSeven
   Deletes this FormulationJFLeeSeven
*/

#endif
