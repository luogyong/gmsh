#ifndef _FORMULATIONOSRCTHREE_H_
#define _FORMULATIONOSRCTHREE_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"

#include "Formulation.h"
#include "FormulationOSRC.h"

/**
   @class FormulationOSRCThree
   @brief Helping class for FormulationOSRC (uncoupled with auxiliary as unknowns)

   Helping class for FormulationOSRC (auxiliary is unknown and tested by itself)

   FormulationOSRC is a friend of FormulationOSRCThree
*/

class FormulationOSRCThree: public Formulation<Complex>{
 private:
  friend class FormulationOSRC;

 private:
  // Complexified wavenumber //
  Complex keps;

  // Pade Bj //
  Complex Bj;

  // FunctionSpace & Domain //
  const FunctionSpaceScalar* faux;
  const GroupOfElement*      ddomain;

  // Local Terms //
  const TermFieldField<double>* localFF;
  const TermGradGrad<double>*   localGG;

 private:
  FormulationOSRCThree(void);
  FormulationOSRCThree(const GroupOfElement& domain,
                       const FunctionSpaceScalar& auxiliary,
                       Complex keps,
                       int NPade,
                       int jPade,
                       const TermFieldField<double>& localFF,
                       const TermGradGrad<double>& localGG);

 public:
  virtual ~FormulationOSRCThree(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationOSRCThree::~FormulationOSRCThree
   Deletes this FormulationOSRCThree
*/

#endif
