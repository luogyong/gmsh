#ifndef _FORMULATIONPML_H_
#define _FORMULATIONPML_H_

#include "SmallFem.h"
#include "Formulation.h"
#include "FormulationSteadyWave.h"

/**
   @class FormulationPML
   @brief Formulation for the PML problem

   Formulation for the PML problem
 */

class FormulationPML: public Formulation<Complex>{
 private:
  // Wavenumber //
  double k;

  // Formulations for steady wave //
  FormulationSteadyWave<Complex>* wave;

 public:
  FormulationPML(const GroupOfElement& domain,
                 const FunctionSpace& fs,
                 double k);

  virtual ~FormulationPML(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;
};

/**
   @fn FormulationPML::FormulationPML
   @param domain A GroupOfElement
   @param fs A FunctionSpace  for both unknown and test field
   @param k a scalar

   Instantiates a new FormulationPML with given parametres:
   @li domain for the domain of this Formulation
   @li fs for the function space used for the unknown field
       and the test functions
   @li k for wavenumber
   **

   @fn FormulationPML::~FormulationPML
   Deletes this FormulationPML
*/

#endif
