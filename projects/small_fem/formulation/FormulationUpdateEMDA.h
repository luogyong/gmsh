#ifndef _FORMULATIONUPDATEEMDA_H_
#define _FORMULATIONUPDATEEMDA_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"
#include "Formulation.h"

/**
   @class FormulationUpdateEMDA
   @brief Update Formulation for FormulationEMDA

   Update Formulation for FormulationEMDA
 */

class FormulationUpdateEMDA: public Formulation<Complex>{
 private:
  // Wavenumber & Chi //
  double k;
  double chi;

  // Function Space & Domain //
  const FunctionSpaceScalar* fspace;
  const GroupOfElement*      goe;

  // Local Terms //
  TermFieldField*               lGout;
  TermProjectionField<Complex>* lGin;
  TermProjectionField<Complex>* lU;

 public:
  FormulationUpdateEMDA(const GroupOfElement& domain,
                        const FunctionSpaceScalar& fs,
                        double k,
                        double chi,
                        const std::map<Dof, Complex>& sol,
                        const std::map<Dof, Complex>& oldG);

  virtual ~FormulationUpdateEMDA(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationUpdateEMDA::FormulationUpdateEMDA
   @todo TODO
   **

   @fn FormulationUpdateEMDA::~FormulationUpdateEMDA
   Deletes this FormulationUpdateEMDA
*/

#endif
