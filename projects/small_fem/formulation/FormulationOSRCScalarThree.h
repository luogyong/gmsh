#ifndef _FORMULATIONOSRCSCALARTHREE_H_
#define _FORMULATIONOSRCSCALARTHREE_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"

#include "FormulationBlock.h"
#include "FormulationOSRCScalar.h"

/**
   @class FormulationOSRCScalarThree
   @brief Helping class for FormulationOSRCScalar (uncoupled with auxiliary)

   Helping class for FormulationOSRCScalar
   (auxiliary is unknown and tested by itself)

   FormulationOSRCScalar is a friend of FormulationOSRCScalarThree
*/

class FormulationOSRCScalarThree: public FormulationBlock<Complex>{
 private:
  friend class FormulationOSRCScalar;

 private:
  // Complexified wavenumber //
  Complex keps;

  // Pade Bj //
  Complex Bj;

  // FunctionSpace & Domain //
  const FunctionSpace*  faux;
  const GroupOfElement* ddomain;

  // Local Terms //
  const TermFieldField<double>* localFF;
  const TermGradGrad<double>*   localGG;

 private:
  FormulationOSRCScalarThree(void);
  FormulationOSRCScalarThree(const GroupOfElement& domain,
                             const FunctionSpace& auxiliary,
                             Complex keps,
                             int NPade,
                             int jPade,
                             const TermFieldField<double>& localFF,
                             const TermGradGrad<double>& localGG);

 public:
  virtual ~FormulationOSRCScalarThree(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationOSRCScalarThree::~FormulationOSRCScalarThree
   Deletes this FormulationOSRCScalarThree
*/

#endif
