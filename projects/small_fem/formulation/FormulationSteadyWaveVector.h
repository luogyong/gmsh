#ifndef _FORMULATIONSTEADYWAVEVECTOR_H_
#define _FORMULATIONSTEADYWAVEVECTOR_H_

#include "FunctionSpaceVector.h"

#include "TermCurlCurl.h"
#include "TermGradGrad.h"

#include "Formulation.h"

/**
   @class FormulationSteadyWaveVector
   @brief Vectorial Formulation for the Steady Wave problem

   Vectorial Formulation for the steady wave problem
 */

class FormulationSteadyWaveVector: public Formulation<double>{
 private:
  // Wavenumber Squared //
  double kSquare;

  // Function Space & Domain //
  const FunctionSpaceVector* fspace;
  const GroupOfElement*      goe;

  // Local Terms //
  TermCurlCurl* localTerms1;
  TermGradGrad* localTerms2;

 public:
  FormulationSteadyWaveVector(const GroupOfElement& goe,
                              const FunctionSpaceVector& fs,
                              double k);

  virtual ~FormulationSteadyWaveVector(void);

  virtual bool isGeneral(void) const;

  virtual double weak(size_t dofI, size_t dofJ, size_t elementId)  const;
  virtual double weakB(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual double rhs(size_t equationI, size_t elementId)           const;

  virtual const FunctionSpace&  fs(void)     const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationSteadyWaveVector::FormulationSteadyWaveVector
   @param goe A GroupOfElement
   @param k A real number
   @param order A natural number

   Instantiates a new FormulationSteadyWaveVector of the given
   order and wavenumber (k)@n

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationSteadyWaveVector::~FormulationSteadyWaveVector
   Deletes this FormulationSteadyWaveVector
*/

#endif
