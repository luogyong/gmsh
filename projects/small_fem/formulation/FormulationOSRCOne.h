#ifndef _FORMULATIONOSRCONE_H_
#define _FORMULATIONOSRCONE_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"

#include "FormulationBlock.h"
#include "FormulationOSRC.h"

/**
   @class FormulationOSRCOne
   @brief Helping class for FormulationOSRC (uncoupled with field as unknowns)

   Helping class for FormulationOSRC (field is unknown and tested by itself)

   FormulationOSRC is a friend of FormulationOSRCOne
 */

class FormulationOSRCOne: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRC;

 private:
  // Wavenumber //
  double k;

  // Pade C0 //
  Complex C0;

  // Function Space & Domain //
  const FunctionSpace*  ffield;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermFieldField<double>*       localLHS;
  const TermProjectionField<Complex>* localRHS;

 private:
  FormulationOSRCOne(void);
  FormulationOSRCOne(const GroupOfElement& domain,
                     const FunctionSpace& field,
                     double k,
                     int NPade,
                     const TermFieldField<double>& localLHS,
                     const TermProjectionField<Complex>& localRHS);

 public:
  virtual ~FormulationOSRCOne(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;

 private:
  void update(TermProjectionField<Complex>& localRHS);
};

/**
   @fn FormulationOSRCOne::~FormulationOSRCOne
   Deletes this FormulationOSRCOne
*/

#endif
