#ifndef _FORMULATIONUPDATEOSRC_H_
#define _FORMULATIONUPDATEOSRC_H_

#include <map>

#include "SmallFem.h"
#include "Formulation.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"

/**
   @class FormulationUpdateOSRC
   @brief Update Formulation for FormulationOSRC

   Update Formulation for FormulationOSRC
*/

class FormulationUpdateOSRC: public Formulation<Complex>{
 private:
  // Wavenumber //
  double k;

  // Function Space & Domain //
  const FunctionSpaceScalar* ffspace;
  const GroupOfElement*      ddomain;

  // Pade //
  Complex C0;
  Complex A1;
  Complex B1;

  // Local Terms //
  TermFieldField<double>*       lGout;
  TermProjectionField<Complex>* lGin;
  TermProjectionField<Complex>* lC0;
  TermProjectionField<Complex>* lAB;

 public:
  FormulationUpdateOSRC(const GroupOfElement& domain,
                        const FunctionSpaceScalar& fspace,
                        double k,
                        const std::map<Dof, Complex>& solU,
                        const std::map<Dof, Complex>& solPhi,
                        const std::map<Dof, Complex>& oldG);

  virtual ~FormulationUpdateOSRC(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationUpdateOSRC::FormulationUpdateOSRC
   @todo TODO
   **

   @fn FormulationUpdateOSRC::~FormulationUpdateOSRC
   Deletes this FormulationUpdateOSRC
*/

#endif
