#ifndef _FORMULATIONOSRCTHREE_H_
#define _FORMULATIONOSRCTHREE_H_

#include <map>

#include "SmallFem.h"
#include "Formulation.h"
#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"

/**
   @class FormulationOSRCThree
   @brief FormulationOSRCThree

   FormulationOSRCThree
*/

class FormulationOSRCThree: public Formulation<Complex>{
 private:
  // Complexified wavenumber //
  Complex keps;

  // FunctionSpace & Domain //
  const FunctionSpaceScalar* ffspace;
  const GroupOfElement*      ddomain;

  // Pade B1 //
  Complex B1;

  // Local Terms //
  TermFieldField<double>* localFF;
  TermGradGrad<double>*   localGG;

 public:
  FormulationOSRCThree(const GroupOfElement& domain,
                       const FunctionSpaceScalar& fspace,
                       Complex keps);

  virtual ~FormulationOSRCThree(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

 private:
  static double pade_aj(int j, int N);
  static double pade_bj(int j, int N);

  static Complex padeBj(int j, int N, double theta);
};

#endif
