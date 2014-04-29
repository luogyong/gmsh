#ifndef _FORMULATIONUPDATEOSRC_H_
#define _FORMULATIONUPDATEOSRC_H_

#include <map>

#include "SmallFem.h"
#include "FormulationBlock.h"
#include "FunctionSpaceScalar.h"
#include "TermProjectionField.h"
#include "TermFieldField.h"

/**
   @class FormulationUpdateOSRC
   @brief Update Formulation for FormulationOSRC

   Update Formulation for FormulationOSRC
*/

class FormulationUpdateOSRC: public FormulationBlock<Complex>{
 private:
  // Wavenumber //
  double k;

  // Function Space & Domain //
  const FunctionSpaceScalar* ffspace;
  const GroupOfElement*      ddomain;

  // Pade //
  Complex C0;
  std::vector<Complex> A;
  std::vector<Complex> B;

  // Local Terms //
  TermFieldField<double>*       lGout;
  TermProjectionField<Complex>* lGin;
  TermProjectionField<Complex>* lC0;
  TermProjectionField<Complex>* lAB;

 public:
  FormulationUpdateOSRC(const GroupOfElement& domain,
                        const FunctionSpaceScalar& fspace,
                        double k,
                        int NPade,
                        const std::map<Dof, Complex>& solU,
                        const std::vector<std::map<Dof, Complex> >& solPhi,
                        const std::map<Dof, Complex>& oldG);

  virtual ~FormulationUpdateOSRC(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationUpdateOSRC::FormulationUpdateOSRC
   @todo TODO
   **

   @fn FormulationUpdateOSRC::~FormulationUpdateOSRC
   Deletes this FormulationUpdateOSRC
*/

#endif
