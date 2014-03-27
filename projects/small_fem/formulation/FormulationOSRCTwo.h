#ifndef _FORMULATIONOSRCTWO_H_
#define _FORMULATIONOSRCTWO_H

#include "SmallFem.h"
#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "TermGradGrad.h"

#include "Formulation.h"
#include "FormulationOSRC.h"

/**
   @class FormulationOSRCTwo
   @brief Helping class for FormulationOSRC (coupled with auxiliary as unknowns)

   Helping class for FormulationOSRC (auxiliary is unknown and tested by field)

   FormulationOSRC is a friend of FormulationOSRCTwo
 */

class FormulationOSRCTwo: public Formulation<Complex>{
 private:
  friend class FormulationOSRC;

 private:
  // Wavenumber (normal and complexified) //
  double  k;
  Complex keps;

  // Pade Aj //
  Complex Aj;

  // Function Space (field and aux) & Domain //
  const FunctionSpaceScalar* ffield;
  const FunctionSpaceScalar* faux;
  const GroupOfElement*      ddomain;

  // Local Term //
  const TermGradGrad<double>* localTerm;

 private:
  FormulationOSRCTwo(void);
  FormulationOSRCTwo(const GroupOfElement& domain,
                     const FunctionSpaceScalar& auxiliary,
                     const FunctionSpaceScalar& field,
                     double  k,
                     Complex keps,
                     int NPade,
                     int jPade,
                     const TermGradGrad<double>& localTerm);

 public:
  virtual ~FormulationOSRCTwo(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationOSRCTwo::~FormulationOSRCTwo
   Deletes this FormulationOSRCTwo
*/

#endif
