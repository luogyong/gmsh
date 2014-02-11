#ifndef _FORMULATIONNEUMANN_H_
#define _FORMULATIONNEUMANN_H_

#include <complex>
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "Formulation.h"

/**
   @class FormulationNeumann
   @brief Neumann Formulation

   Neumann Formulation
 */

class FormulationNeumann: public Formulation<std::complex<double> >{
 private:
  // Wavenumber //
  double k;

  // Function Space & Domain //
  const FunctionSpaceScalar* fspace;
  const GroupOfElement*      goe;

  // Local Terms //
  TermFieldField* localTerms;

 public:
  FormulationNeumann(const GroupOfElement& goe,
                     const FunctionSpaceScalar& fs,
                     double k);

  virtual ~FormulationNeumann(void);

  virtual std::complex<double>
    weak(size_t dofI, size_t dofJ, size_t elementId) const;

  virtual std::complex<double>
    rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationNeumann::FormulationNeumann
   @param goe A GroupOfElement
   @param k A real number
   @param order A natural number

   Instantiates a new FormulationNeumann of the given
   order and wavenumber (k)@n

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationNeumann::~FormulationNeumann
   Deletes this FormulationNeumann
*/

#endif
