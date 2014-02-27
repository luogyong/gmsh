#ifndef _FORMULATIONPML_H_
#define _FORMULATIONPML_H_

#include "TermFieldField.h"
#include "TermGradGrad.h"
#include "Formulation.h"
#include "SmallFem.h"

/**
   @class FormulationPML
   @brief Formulation for the PML problem

   Formulation for the PML problem
 */

class FormulationPML: public Formulation<Complex>{
 private:
  // Function Space & Domain //
  const FunctionSpace*  ffs;
  const GroupOfElement* ddomain;

  // Wavenumber squared //
  double kSquare;

  // Local Terms //
  TermGradGrad<Complex>*   stif;
  TermFieldField<Complex>* mass;

 public:
  FormulationPML(const GroupOfElement& domain,
                 const FunctionSpace& fs,
                 double k,
                 void    (*fS)(fullVector<double>& xyz, fullMatrix<Complex>& T),
                 Complex (*fM)(fullVector<double>& xyz));

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
   @param k A scalar
   @param fS A tensorial function
   @param fM A scalar function

   Instantiates a new FormulationPML with given parametres:
   @li domain for the domain of this Formulation
   @li fs for the function space used for the unknown field
       and the test functions
   @li k for wavenumber
   @li fS for the PML stiffness
   @li fM for the PML mass
   **

   @fn FormulationPML::~FormulationPML
   Deletes this FormulationPML
*/

#endif
