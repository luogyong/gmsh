#ifndef _FORMULATIONOSRCVECTORSEVEN_H_
#define _FORMULATIONOSRCVECTORSEVEN_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorSeven
   @brief Helping class for FormulationOSRCVector <grad(rho), phi>

   Helping class for FormulationOSRCVector <grad(rho), phi>

   FormulationOSRCVector is a friend of FormulationOSRCVectorSeven
 */

class FormulationOSRCVectorSeven: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Pade //
  Complex minusBi;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermGradGrad<double>* localGG;

 private:
  FormulationOSRCVectorSeven(void);
  FormulationOSRCVectorSeven(const GroupOfElement& domain,
                             const FunctionSpace& field,
                             const FunctionSpace& test,
                             Complex Bi,
                             const TermGradGrad<double>& localGG);

 public:
  virtual ~FormulationOSRCVectorSeven(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCVectorSeven::~FormulationOSRCVectorSeven
   Deletes this FormulationOSRCVectorSeven
*/

#endif
