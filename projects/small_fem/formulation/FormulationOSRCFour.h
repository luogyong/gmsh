#ifndef _FORMULATIONOSRCFOUR_H_
#define _FORMULATIONOSRCFOUR_H_

#include <map>

#include "SmallFem.h"
#include "Formulation.h"
#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"

/**
   @class FormulationOSRCFour
   @brief FormulationOSRCFour

   FormulationOSRCFour
 */

class FormulationOSRCFour: public Formulation<Complex>{
 private:
  // FunctionSpace (field and test) & Domain //
  const FunctionSpaceScalar* ffField;
  const FunctionSpaceScalar* ffTest;
  const GroupOfElement*      ddomain;

  // Local Term //
  TermFieldField<double>* local;

 public:
  FormulationOSRCFour(const GroupOfElement& domain,
                      const FunctionSpaceScalar& fField,
                      const FunctionSpaceScalar& fTest);

  virtual ~FormulationOSRCFour(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

#endif
