#ifndef _FORMULATIONJFLEEFOUR_H_
#define _FORMULATIONJFLEEFOUR_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationJFLee.h"

/**
   @class FormulationJFLeeFour
   @brief Helping class for FormulationJFLee (<phi, phi> term)

   Helping class for FormulationJFLee (<phi, phi> term)

   FormulationJFLee is a friend of FormulationJFLeeFour
 */

class FormulationJFLeeFour: public FormulationBlock<Complex>{
 private:
  friend class FormulationJFLee;

 private:
  // Wavenumber //
  double kSquare;

  // Function Space & Domain //
  const FunctionSpaceVector* ffield;
  const GroupOfElement*      ddomain;

  // Local Terms //
  const TermGradGrad<double>* local;

 private:
  FormulationJFLeeFour(void);
  FormulationJFLeeFour(const GroupOfElement& domain,
                       const FunctionSpaceVector& field,
                       double k,
                       const TermGradGrad<double>& local);

 public:
  virtual ~FormulationJFLeeFour(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationJFLeeFour::~FormulationJFLeeFour
   Deletes this FormulationJFLeeFour
*/

#endif
