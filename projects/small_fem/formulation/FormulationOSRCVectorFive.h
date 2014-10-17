#ifndef _FORMULATIONOSRCVECTORFIVE_H_
#define _FORMULATIONOSRCVECTORFIVE_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorFive
   @brief Helping class for FormulationOSRCVector <r, phi>

   Helping class for FormulationOSRCVector <r, phi>

   FormulationOSRCVector is a friend of FormulationOSRCVectorFive
 */

class FormulationOSRCVectorFive: public FormulationBlock<Complex>{
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
  FormulationOSRCVectorFive(void);
  FormulationOSRCVectorFive(const GroupOfElement& domain,
                            const FunctionSpace& field,
                            const FunctionSpace& test,
                            const TermGradGrad<double>& localGG);

 public:
  virtual ~FormulationOSRCVectorFive(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCVectorFive::~FormulationOSRCVectorFive
   Deletes this FormulationOSRCVectorFive
*/

#endif
