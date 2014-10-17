#ifndef _FORMULATIONOSRCVECTORFOUR_H_
#define _FORMULATIONOSRCVECTORFOUR_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorFour
   @brief Helping class for FormulationOSRCVector <phi, r>

   Helping class for FormulationOSRCVector <phi, r>

   FormulationOSRCVector is a friend of FormulationOSRCVectorFour
 */

class FormulationOSRCVectorFour: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Pade //
  Complex alpha;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermGradGrad<double>* localGG;

 private:
  FormulationOSRCVectorFour(void);
  FormulationOSRCVectorFour(const GroupOfElement& domain,
                            const FunctionSpace& field,
                            const FunctionSpace& test,
                            Complex R0,
                            Complex Ai,
                            Complex Bi,
                            const TermGradGrad<double>& localGG);

 public:
  virtual ~FormulationOSRCVectorFour(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCVectorFour::~FormulationOSRCVectorFour
   Deletes this FormulationOSRCVectorFour
*/

#endif
