#ifndef _FORMULATIONOSRCSCALARTWO_H_
#define _FORMULATIONOSRCSCALARTWO_H

#include "SmallFem.h"
#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationOSRCScalar.h"

/**
   @class FormulationOSRCScalarTwo
   @brief Helping class for FormulationOSRCScalar (coupled with auxiliary)

   Helping class for FormulationOSRCScalar
   (auxiliary is unknown and tested by field)

   FormulationOSRCScalar is a friend of FormulationOSRCScalarTwo
 */

class FormulationOSRCScalarTwo: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCScalar;

 private:
  // Wavenumber (normal and complexified) //
  double  k;
  Complex keps;

  // Pade Aj //
  Complex Aj;

  // Function Space (field and aux) & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  faux;
  const GroupOfElement* ddomain;

  // Local Term //
  const TermGradGrad<double>* localTerm;

 private:
  FormulationOSRCScalarTwo(void);
  FormulationOSRCScalarTwo(const GroupOfElement& domain,
                           const FunctionSpace& auxiliary,
                           const FunctionSpace& field,
                           double  k,
                           Complex keps,
                           Complex Aj,
                           const TermGradGrad<double>& localTerm);

 public:
  virtual ~FormulationOSRCScalarTwo(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCScalarTwo::~FormulationOSRCScalarTwo
   Deletes this FormulationOSRCScalarTwo
*/

#endif
