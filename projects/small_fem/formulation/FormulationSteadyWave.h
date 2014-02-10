#ifndef _FORMULATIONSTEADYWAVE_H_
#define _FORMULATIONSTEADYWAVE_H_

#include "Formulation.h"
#include "FormulationMass.h"
#include "FormulationStiffness.h"

/**
   @class FormulationSteadyWave
   @brief Formulation for the steady wave problem

   Formulation for the steady wave problem
 */

template<typename scalar>
class FormulationSteadyWave: public Formulation<scalar>{
 private:
  // Wavenumber Squared //
  scalar kSquare;

  // Formulations for stiffness and mass terms //
  FormulationStiffness<scalar>* stiff;
  FormulationMass<scalar>*      mass;

 public:
  FormulationSteadyWave(const GroupOfElement& domain,
                        const FunctionSpace& fs,
                        scalar k);

  virtual ~FormulationSteadyWave(void);

  virtual bool isGeneral(void) const;

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId)  const;
  virtual scalar weakB(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual scalar rhs(size_t equationI, size_t elementId)           const;

  virtual const FunctionSpace&  fsField(void) const;
  virtual const FunctionSpace&  fsTest(void)  const;
  virtual const GroupOfElement& domain(void)  const;
};

/**
   @fn FormulationSteadyWave::FormulationSteadyWave
   @param domain A GroupOfElement
   @param fs A FunctionSpace
   @param k a scalar

   Instantiates a new FormulationSteadyWave with given parametres:
   @li domain for the domain of this Formulation
   @li fs for the function space used for the unknown field
       and the test functions
   @li k for wavenumber
   **

   @fn FormulationSteadyWave::~FormulationSteadyWave
   Deletes this FormulationSteadyWave
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationSteadyWaveInclusion.h"

#endif
