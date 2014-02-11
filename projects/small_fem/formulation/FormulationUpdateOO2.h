#ifndef _FORMULATIONUPDATEOO2_H_
#define _FORMULATIONUPDATEOO2_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "TermGradGrad.h"
#include "Formulation.h"

/**
   @class FormulationUpdateOO2
   @brief Update Formulation for FormulationOO2

   Update Formulation for FormulationOO2
 */

class FormulationUpdateOO2: public Formulation<Complex>{
 private:
  // a & b //
  Complex a;
  Complex b;

  // Function Space & Basis //
  const FunctionSpaceScalar* fspace;
  const Basis*                basis;

  // Domain //
  const GroupOfElement* goe;

  // Quadrature (Field - Field) //
  fullMatrix<double>* gCFF;
  fullVector<double>* gWFF;
  GroupOfJacobian*    jacFF;

  // Quadrature (Grad - Grad) //
  fullMatrix<double>* gCGG;
  fullVector<double>* gWGG;
  GroupOfJacobian*    jacGG;

  // DDM //
  const std::map<Dof, Complex>* solution;
  const std::map<Dof, Complex>* oldG;

 public:
  FormulationUpdateOO2(const GroupOfElement& domain,
                       const FunctionSpaceScalar& fs,
                       Complex a,
                       Complex b,
                       const std::map<Dof, Complex>& solution,
                       const std::map<Dof, Complex>& oldG);

  virtual ~FormulationUpdateOO2(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

 private:
  Complex interpolate(const MElement& element,
                      const fullVector<double>& xyz,
                      const std::map<Dof, Complex>& f) const;

  fullVector<Complex> interpolateGrad(const MElement& element,
                                      const fullVector<double>& xyz,
                                      const std::map<Dof, Complex>& f) const;
};

/**
   @fn FormulationUpdateOO2::FormulationUpdateOO2
   @todo TODO
   **

   @fn FormulationUpdateOO2::~FormulationUpdateOO2
   Deletes this FormulationUpdateOO2
*/

#endif
