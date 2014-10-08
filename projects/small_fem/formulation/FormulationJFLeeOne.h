#ifndef _FORMULATIONJFLEEONE_H_
#define _FORMULATIONJFLEEONE_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermProjectionGrad.h"

#include "FormulationBlock.h"
#include "FormulationJFLee.h"

/**
   @class FormulationJFLeeOne
   @brief Helping class for FormulationJFLee (<g, e> term)

   Helping class for FormulationJFLee (<g, e> term)

   FormulationJFLee is a friend of FormulationJFLeeOne
 */

class FormulationJFLeeOne: public FormulationBlock<Complex>{
 private:
  friend class FormulationJFLee;

 private:
  // Wavenumber //
  double k;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermProjectionGrad<Complex>* localRHS;

 private:
  FormulationJFLeeOne(void);
  FormulationJFLeeOne(const GroupOfElement& domain,
                      const FunctionSpace& field,
                      double k,
                      const TermProjectionGrad<Complex>& localRHS);

 public:
  virtual ~FormulationJFLeeOne(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;

 private:
  void update(TermProjectionGrad<Complex>& localRHS);
};

/**
   @fn FormulationJFLeeOne::~FormulationJFLeeOne
   Deletes this FormulationJFLeeOne
*/

#endif
