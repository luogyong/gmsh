#ifndef _FORMULATIONSILVERMULLER_H_
#define _FORMULATIONSILVERMULLER_H_

#include "SmallFem.h"
#include "Term.h"
#include "FunctionSpace.h"
#include "FormulationBlock.h"

/**
   @class FormulationSilverMuller
   @brief Silver-Muller radiation formulation

   Weak formulation for the Silver-Muller radiation condition
 */

class FormulationSilverMuller: public FormulationBlock<Complex>{
 private:
  // Wavenumber //
  double k;

  // Function Space & Domain //
  const FunctionSpace*  fspace;
  const GroupOfElement* goe;

  // Local Terms //
  Term<double>* localTerms;

 public:
  FormulationSilverMuller(const GroupOfElement& domain,
                          const FunctionSpace& fs,
                          double k);

  virtual ~FormulationSilverMuller(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationSilverMuller::FormulationSilverMuller
   @param domain A GroupOfElement for the domain
   @param fs A FunctionSpace for both the unknown and test field
   @param k A real number

   Instantiates a new FormulationSilverMuller with wavenumber k
   **

   @fn FormulationSilverMuller::~FormulationSilverMuller
   Deletes this FormulationSilverMuller
*/

#endif
