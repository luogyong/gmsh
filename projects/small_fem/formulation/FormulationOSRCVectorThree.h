#ifndef _FORMULATIONOSRCVECTORTHREE_H_
#define _FORMULATIONOSRCVECTORTHREE_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermProjectionGrad.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorThree
   @brief Helping class for FormulationOSRCVector <r, r> - <g, r>

   Helping class for FormulationOSRCVector <r, r> - <g, r>

   FormulationOSRCVector is a friend of FormulationOSRCVectorThree
 */

class FormulationOSRCVectorThree: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Pade //
  Complex oneOverR0;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermGradGrad<double>*       localGG;
  const TermProjectionGrad<Complex>* localRHS;

 private:
  FormulationOSRCVectorThree(void);
  FormulationOSRCVectorThree(const GroupOfElement& domain,
                             const FunctionSpace& field,
                             const FunctionSpace& test,
                             Complex R0,
                             const TermGradGrad<double>& localGG,
                             const TermProjectionGrad<Complex>& localRHS);

 public:
  virtual ~FormulationOSRCVectorThree(void);

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
   @fn FormulationOSRCVectorThree::~FormulationOSRCVectorThree
   Deletes this FormulationOSRCVectorThree
*/

#endif
