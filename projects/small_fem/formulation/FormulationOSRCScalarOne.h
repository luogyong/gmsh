#ifndef _FORMULATIONOSRCSCALARONE_H_
#define _FORMULATIONOSRCSCALARONE_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"

#include "FormulationBlock.h"
#include "FormulationOSRCScalar.h"

/**
   @class FormulationOSRCScalarOne
   @brief Helping class for FormulationOSRCScalar (uncoupled with field)

   Helping class for FormulationOSRCScalar
   (field is unknown and tested by itself)

   FormulationOSRCScalar is a friend of FormulationOSRCScalarOne
 */

class FormulationOSRCScalarOne: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCScalar;

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
  FormulationOSRCScalarOne(void);
  FormulationOSRCScalarOne(const GroupOfElement& domain,
                           const FunctionSpace& field,
                           double k,
                           Complex C0,
                           const TermFieldField<double>& localLHS,
                           const TermProjectionField<Complex>& localRHS);

 public:
  virtual ~FormulationOSRCScalarOne(void);

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
   @fn FormulationOSRCScalarOne::~FormulationOSRCScalarOne
   Deletes this FormulationOSRCScalarOne
*/

#endif
