#ifndef _FORMULATIONSTEADYWAVEVECTORSLOW_H_
#define _FORMULATIONSTEADYWAVEVECTORSLOW_H_

#include <vector>

#include "FunctionSpaceVector.h"
#include "GroupOfJacobian.h"
#include "Formulation.h"

/**
   @class FormulationSteadyWaveVectorSlow
   @brief Vectorial Formulation for the Steady Wave problem (Slow version)

   Same as FormulationSteadyWaveVector, but this version don't use
   the fast integration algorithm
 */

class FormulationSteadyWaveVectorSlow: public Formulation<double>{
 private:
  // Wavenumber Squared //
  double kSquare;

  // Gaussian Quadrature Data (Term One) //
  int G1;
  fullMatrix<double>* gC1;
  fullVector<double>* gW1;

  // Gaussian Quadrature Data (Term Two) //
  int G2;
  fullMatrix<double>* gC2;
  fullVector<double>* gW2;

  // Domain //
  const GroupOfElement* goe;

  // Jacobians //
  GroupOfJacobian* jac1;
  GroupOfJacobian* jac2;

  // Function Space & Basis //
  const FunctionSpaceVector* fspace;
  const Basis*               basis;

 public:
  FormulationSteadyWaveVectorSlow(const GroupOfElement& goe,
                                  const FunctionSpaceVector& fs,
                                  double k);

  virtual ~FormulationSteadyWaveVectorSlow(void);

  virtual bool isGeneral(void) const;

  virtual double weak(size_t dofI, size_t dofJ, size_t elementId)  const;
  virtual double weakB(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual double rhs(size_t equationI, size_t elementId)           const;

  virtual const FunctionSpace&  fsField(void) const;
  virtual const FunctionSpace&  fsTest(void)  const;
  virtual const GroupOfElement& domain(void)  const;
};

/**
   @fn FormulationSteadyWaveVectorSlow::FormulationSteadyWaveVectorSlow
   @param goe A GroupOfElement
   @param k A real number
   @param order A natural number

   Instantiates a new FormulationSteadyWaveVectorSlow of the given
   order and wavenumber (k)@n

   The given GroupOfElement will be used as the geomtrical domain
   **

   @fn FormulationSteadyWaveVectorSlow::~FormulationSteadyWaveVectorSlow
   Deletes this FormulationSteadyWaveVectorSlow
*/

#endif
