#ifndef _FORMULATIONOSRCSCALARFOUR_H_
#define _FORMULATIONOSRCSCALARFOUR_H_

#include "SmallFem.h"
#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"

#include "FormulationBlock.h"
#include "FormulationOSRCScalar.h"

/**
   @class FormulationOSRCScalarFour
   @brief Helping class for FormulationOSRCScalar
   (coupled with field)

   Helping class for FormulationOSRCScalar
   (field is unknown and tested by auxiliary)

   FormulationOSRCScalar is a friend of FormulationOSRCScalarFour
 */

class FormulationOSRCScalarFour: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCScalar;

 private:
  // FunctionSpace (field and aux) & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  faux;
  const GroupOfElement* ddomain;

  // Local Term //
  const TermFieldField<double>* localTerm;

 private:
  FormulationOSRCScalarFour(void);
  FormulationOSRCScalarFour(const GroupOfElement& domain,
                            const FunctionSpace& field,
                            const FunctionSpace& auxiliary,
                            const TermFieldField<double>& localTerm);

 public:
  virtual ~FormulationOSRCScalarFour(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCScalarFour::~FormulationOSRCScalarFour
   Deletes this FormulationOSRCScalarFour
*/

#endif
