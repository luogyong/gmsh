#ifndef _FORMULATIONUPDATEEMDA_H_
#define _FORMULATIONUPDATEEMDA_H_

#include "SmallFem.h"
#include "FunctionSpaceScalar.h"
#include "TermFieldField.h"
#include "Formulation.h"

/**
   @class FormulationUpdateEMDA
   @brief Update Formulation for FormulationEMDA

   Update Formulation for FormulationEMDA
 */

class FormulationUpdateEMDA: public Formulation<Complex>{
 private:
  // Wavenumber & Chi //
  double k;
  double chi;

  // Function Space & Basis //
  const FunctionSpaceScalar* fspace;
  const Basis*                basis;

  // Domain //
  const GroupOfElement* goe;

  // Quadrature //
  fullMatrix<double>* gC;
  fullVector<double>* gW;
  GroupOfJacobian*    jac;

  // DDM //
  const std::map<Dof, Complex>* solution;
  const std::map<Dof, Complex>* oldG;

 public:
  FormulationUpdateEMDA(const GroupOfElement& domain,
                        const FunctionSpaceScalar& fs,
                        double k,
                        double chi,
                        const std::map<Dof, Complex>& solution,
                        const std::map<Dof, Complex>& oldG);

  virtual ~FormulationUpdateEMDA(void);

  virtual Complex weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual Complex rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

 private:
  Complex
    interpolate(const MElement& element,
                const fullVector<double>& xyz,
                const std::map<Dof, Complex>& f) const;
};

/**
   @fn FormulationUpdateEMDA::FormulationUpdateEMDA
   @todo TODO
   **

   @fn FormulationUpdateEMDA::~FormulationUpdateEMDA
   Deletes this FormulationUpdateEMDA
*/

#endif
