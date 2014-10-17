#ifndef _FORMULATIONOSRCVECTOREIGHT_H_
#define _FORMULATIONOSRCVECTOREIGHT_H_

#include "SmallFem.h"
#include "FunctionSpaceVector.h"
#include "TermFieldField.h"

#include "FormulationBlock.h"
#include "FormulationOSRCVector.h"

/**
   @class FormulationOSRCVectorEight
   @brief Helping class for FormulationOSRCVector <rho, rho>

   Helping class for FormulationOSRCVector <rho, rho>

   FormulationOSRCVector is a friend of FormulationOSRCVectorEight
 */

class FormulationOSRCVectorEight: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCVector;

 private:
  // Wavenumber //
  Complex kEpsSquare;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  ttest;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermFieldField<double>* localFF;

 private:
  FormulationOSRCVectorEight(void);
  FormulationOSRCVectorEight(const GroupOfElement& domain,
                             const FunctionSpace& field,
                             const FunctionSpace& test,
                             Complex kEps,
                             const TermFieldField<double>& localFF);

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
