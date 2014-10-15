#ifndef _FORMULATIONOSRCVECTOREIGHT_H_
#define _FORMULATIONOSRCVECTOREIGHT_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorEight
   @brief Helping class for FormulationOSRCVector <phi, grad(rho)>

   Helping class for FormulationOSRCVector <phi, grad(rho)>

   FormulationOSRCVector is a friend of FormulationOSRCVectorEight
 */

class FormulationOSRCVectorEight: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermGradGrad<double>* localGG;

 private:
  FormulationOSRCVectorEight(void);
  FormulationOSRCVectorEight(const GroupOfElement& domain,
                             const FunctionSpace& field,
                             const FunctionSpace& test,
                             const TermGradGrad<double>& localGG);

 public:
  virtual ~FormulationOSRCVectorEight(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCVectorEight::~FormulationOSRCVectorEight
   Deletes this FormulationOSRCVectorEight
*/

#endif
