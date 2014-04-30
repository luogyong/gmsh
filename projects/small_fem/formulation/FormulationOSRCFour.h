#ifndef _FORMULATIONOSRCFOUR_H_
#define _FORMULATIONOSRCFOUR_H_

#include "SmallFem.h"
#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"

#include "FormulationBlock.h"
#include "FormulationOSRC.h"

/**
   @class FormulationOSRCFour
   @brief Helping class for FormulationOSRC (coupled with field as unknowns)

   Helping class for FormulationOSRC (field is unknown and tested by auxiliary)

   FormulationOSRC is a friend of FormulationOSRCFour
 */

class FormulationOSRCFour: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRC;

 private:
  // FunctionSpace (field and aux) & Domain //
  const FunctionSpace*  ffield;
  const FunctionSpace*  faux;
  const GroupOfElement* ddomain;

  // Local Term //
  const TermFieldField<double>* localTerm;

 private:
  FormulationOSRCFour(void);
  FormulationOSRCFour(const GroupOfElement& domain,
                      const FunctionSpace& field,
                      const FunctionSpace& auxiliary,
                      const TermFieldField<double>& localTerm);

 public:
  virtual ~FormulationOSRCFour(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCFour::~FormulationOSRCFour
   Deletes this FormulationOSRCFour
*/

#endif
