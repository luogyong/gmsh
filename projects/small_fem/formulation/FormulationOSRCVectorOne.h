#ifndef _FORMULATIONOSRCVECTORONE_H_
#define _FORMULATIONOSRCVECTORONE_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorOne
   @brief Helping class for FormulationOSRCVector <r, e>

   Helping class for FormulationOSRCVector <r, e>

   FormulationOSRCVector is a friend of FormulationOSRCVectorOne
 */

class FormulationOSRCVectorOne: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Wavenumber //
  Complex jK;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermGradGrad<double>* localGG;

 private:
  FormulationOSRCVectorOne(void);
  FormulationOSRCVectorOne(const GroupOfElement& domain,
                           const FunctionSpace& field,
                           const FunctionSpace& test,
                           double  k,
                           const TermGradGrad<double>& localGG);

 public:
  virtual ~FormulationOSRCVectorOne(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCVectorOne::~FormulationOSRCVectorOne
   Deletes this FormulationOSRCVectorOne
*/

#endif
