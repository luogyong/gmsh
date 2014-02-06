#ifndef _FORMULATIONSTEADYWAVESCALAR_H_
#define _FORMULATIONSTEADYWAVESCALAR_H_

#include "FunctionSpaceScalar.h"

#include "TermGradGrad.h"
#include "TermFieldField.h"

#include "Formulation.h"

/**
   @class FormulationSteadyWaveScalar
   @brief Scalar Formulation for the Steady Wave problem

   Scalar Formulation for the steady wave problem
 */

template<typename scalar>
class FormulationSteadyWaveScalar: public Formulation<scalar>{
 private:
  // Wavenumber squared //
  double kSquare;

  // Function Space & Domain //
  const FunctionSpaceScalar* fspace;
  const GroupOfElement*      goe;

  // Local Terms //
  TermGradGrad*   localTerms1;
  TermFieldField* localTerms2;

 public:
  FormulationSteadyWaveScalar(const GroupOfElement& goe,
                              const FunctionSpaceScalar& fs,
                              double k);

  virtual ~FormulationSteadyWaveScalar(void);

  virtual bool isGeneral(void) const;

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId)  const;
  virtual scalar weakB(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual scalar rhs(size_t equationI, size_t elementId)           const;

  virtual const FunctionSpace&  fsField(void) const;
  virtual const FunctionSpace&  fsTest(void)  const;
  virtual const GroupOfElement& domain(void)  const;
};

/**
   @fn FormulationSteadyWaveScalar::FormulationSteadyWaveScalar
   @param goe A GroupOfElement
   @param k A real number
   @param order A natural number

   Instantiates a new FormulationSteadyWaveScalar of the given
   order and wavenumber (k)@n

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationSteadyWaveScalar::~FormulationSteadyWaveScalar
   Deletes this FormulationSteadyWaveScalar
*/

#endif
