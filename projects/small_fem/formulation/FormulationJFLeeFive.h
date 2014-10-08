#ifndef _FORMULATIONJFLEEFIVE_H_
#define _FORMULATIONJFLEEFIVE_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationJFLee.h"

/**
   @class FormulationJFLeeFive
   @brief Helping class for FormulationJFLee (<e, phi> term)

   Helping class for FormulationJFLee (<e, phi> term)

   FormulationJFLee is a friend of FormulationJFLeeFive
 */

class FormulationJFLeeFive: public FormulationBlock<Complex>{
 private:
  friend class FormulationJFLee;

 private:
  // Wavenumber //
  double kSquare;

  // Function Space & Domain //
  const FunctionSpaceVector* ffield;
  const FunctionSpaceVector* ttest;
  const GroupOfElement*      ddomain;

  // Local Terms //
  const TermGradGrad<double>* local;

 private:
  FormulationJFLeeFive(void);
  FormulationJFLeeFive(const GroupOfElement& domain,
                       const FunctionSpaceVector& field,
                       const FunctionSpaceVector& test,
                       double k,
                       const TermGradGrad<double>& local);

 public:
  virtual ~FormulationJFLeeFive(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationJFLeeFive::~FormulationJFLeeFive
   Deletes this FormulationJFLeeFive
*/

#endif
