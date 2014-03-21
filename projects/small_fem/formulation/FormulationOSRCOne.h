#ifndef _FORMULATIONOSRCONE_H_
#define _FORMULATIONOSRCONE_H_

#include <map>

#include "SmallFem.h"
#include "Formulation.h"
#include "GroupOfElement.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"

/**
   @class FormulationOSRCOne
   @brief FormulationOSRCOne

   FormulationOSRCOne
 */

class FormulationOSRCOne: public Formulation<Complex>{
 private:
  // Wavenumber //
  double k;

  // Function Space & Domain //
  const FunctionSpaceScalar* fspace;
  const GroupOfElement*      ddomain;

  // Pade C0 //
  Complex C0;

  // Local Terms //
  TermFieldField<double>*       localLHS;
  TermProjectionField<Complex>* localRHS;

 public:
  FormulationOSRCOne(const GroupOfElement& domain,
                     const FunctionSpaceScalar& fs,
                     double k,
                     const std::map<Dof, Complex>& ddmDof);

  virtual ~FormulationOSRCOne(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

 private:
  static double pade_aj(int j, int N);
  static double pade_bj(int j, int N);

  static Complex padeC0(int N, double theta);
};

#endif
