#ifndef _FORMULATIONOSRCTWO_H_
#define _FORMULATIONOSRCTWO_H

#include <map>

#include "SmallFem.h"
#include "Formulation.h"
#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "TermGradGrad.h"

/**
   @class FormulationOSRCTwo
   @brief FormulationOSRCTwo

   FormulationOSRCTwo
 */

class FormulationOSRCTwo: public Formulation<Complex>{
 private:
  // Wavenumber (normal and complexified) //
  double  k;
  Complex keps;

  // Function Space (field and test) & Domain //
  const FunctionSpaceScalar* ffField;
  const FunctionSpaceScalar* ffTest;
  const GroupOfElement*      ddomain;

  // Pade A1 //
  Complex A1;

  // Local Term //
  TermGradGrad<double>* local;

 public:
  FormulationOSRCTwo(const GroupOfElement& domain,
                     const FunctionSpaceScalar& fField,
                     const FunctionSpaceScalar& fTest,
                     double  k,
                     Complex keps);

  virtual ~FormulationOSRCTwo(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

 private:
  static double pade_aj(int j, int N);
  static double pade_bj(int j, int N);

  static Complex padeAj(int j, int N, double theta);
};

#endif
