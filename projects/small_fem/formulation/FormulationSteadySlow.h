#ifndef _FORMULATIONSTEADYSLOW_H_
#define _FORMULATIONSTEADYSLOW_H_

#include <vector>

#include "FunctionSpaceVector.h"
#include "GroupOfJacobian.h"
#include "FormulationBlock.h"

/**
   @class FormulationSteadySlow
   @brief Vectorial Formulation for the steady wave problem (slow version)

   Slow version of the vectorial Formulation for the steady wave problem.
   This version don't use the fast integration algorithm.
 */

template<typename scalar>
class FormulationSteadySlow: public FormulationBlock<scalar>{
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
  const GroupOfElement* ddomain;

  // Jacobians //
  GroupOfJacobian* jac1;
  GroupOfJacobian* jac2;

  // Function Space & Basis //
  const FunctionSpace* fspace;
  const Basis*         basis;

 public:
  FormulationSteadySlow(const GroupOfElement& domain,
                        const FunctionSpace& fs,
                        double k);

  virtual ~FormulationSteadySlow(void);

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual scalar rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

  virtual bool isBlock(void) const;
};

/**
   @fn FormulationSteadySlow::FormulationSteadySlow
   @param domain A GroupOfElement
   @param fs A FunctionSpaceVector for both unknown and test field
   @param k a scalar

   Instantiates a new FormulationSteadySlow with given parametres:
   @li domain for the domain of this Formulation
   @li fs for the function space used for the unknown field
       and the test functions
   @li k for wavenumber
   **

   @fn FormulationSteadySlow::~FormulationSteadySlow
   Deletes this FormulationSteadySlow
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationSteadySlowInclusion.h"

#endif
