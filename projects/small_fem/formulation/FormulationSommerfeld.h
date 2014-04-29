#ifndef _FORMULATIONSOMMERFELD_H_
#define _FORMULATIONSOMMERFELD_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "FormulationBlock.h"

/**
   @class FormulationSommerfeld
   @brief Sommerfeld radiation formulation

   Weak formulation for the Sommerfeld radiation condition
 */

class FormulationSommerfeld: public FormulationBlock<Complex>{
 private:
  // Wavenumber //
  double k;

  // Function Space & Domain //
  const FunctionSpaceScalar* fspace;
  const GroupOfElement*      goe;

  // Local Terms //
  TermFieldField<double>* localTerms;

 public:
  FormulationSommerfeld(const GroupOfElement& domain,
                        const FunctionSpaceScalar& fs,
                        double k);

  virtual ~FormulationSommerfeld(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationSommerfeld::FormulationSommerfeld
   @param domain A GroupOfElement for the domain
   @param fs A FunctionSpace for both the unknown and test field
   @param k A real number

   Instantiates a new FormulationSommerfeld with wavenumber k
   **

   @fn FormulationSommerfeld::~FormulationSommerfeld
   Deletes this FormulationSommerfeld
*/

#endif
