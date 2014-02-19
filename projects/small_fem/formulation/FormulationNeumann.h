#ifndef _FORMULATIONNEUMANN_H_
#define _FORMULATIONNEUMANN_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "Formulation.h"

/**
   @class FormulationNeumann
   @brief Neumann Formulation

   Neumann Formulation
 */

class FormulationNeumann: public Formulation<Complex>{
 private:
  // Wavenumber //
  double k;

  // Function Space & Domain //
  const FunctionSpaceScalar* fspace;
  const GroupOfElement*      goe;

  // Local Terms //
  TermFieldField<double>* localTerms;

 public:
  FormulationNeumann(const GroupOfElement& domain,
                     const FunctionSpaceScalar& fs,
                     double k);

  virtual ~FormulationNeumann(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationNeumann::FormulationNeumann
   @param domain A GroupOfElement for the domain
   @param fs A FunctionSpace for both the unknown and test field
   @param k A real number

   Instantiates a new FormulationNeumann with wavenumber k
   **

   @fn FormulationNeumann::~FormulationNeumann
   Deletes this FormulationNeumann
*/

#endif
