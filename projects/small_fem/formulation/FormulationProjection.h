#ifndef _FORMULATIONPROJECTION_H_
#define _FORMULATIONPROJECTION_H_

#include "FunctionSpace.h"
#include "Formulation.h"
#include "fullMatrix.h"
#include "Term.h"

/**
   @class FormulationProjection
   @brief Formulation for the L2 projection of a function

    Formulation for the L2 projection problem
 */

template<typename scalar>
class FormulationProjection: public Formulation<scalar>{
 private:
  // Function Space & Basis & Domain //
  const FunctionSpace*  fspace;
  const GroupOfElement* ddomain;

  // Local Terms //
  Term<double>* localTerms1;
  Term<scalar>* localTerms2;

 public:
  FormulationProjection(const GroupOfElement& domain,
                        const FunctionSpace& fs,
                        scalar (*f)(fullVector<double>& xyz));

  FormulationProjection(const GroupOfElement& domain,
                        const FunctionSpace& fs,
                        fullVector<scalar> (*g)(fullVector<double>& xyz));

  virtual ~FormulationProjection(void);

  virtual scalar weak(size_t dofI, size_t dofJ, size_t elementId) const;
  virtual scalar rhs(size_t equationI, size_t elementId)          const;

  virtual const FunctionSpace&  field(void)  const;
  virtual const FunctionSpace&  test(void)   const;
  virtual const GroupOfElement& domain(void) const;

 private:
  const Basis& initCommon(const GroupOfElement& domain,
                          const FunctionSpace& fs,
                          fullMatrix<double>& gC,
                          fullVector<double>& gW);
};

/**
   @fn FormulationProjection::FormulationProjection(const GroupOfElement&,const FunctionSpace&,scalar(*f)(fullVector<double>& xyz))
   @param domain The domain of this Formulation
   @param fs A FunctionSpace for both unknown and test field
   @param f The @em scalar function to project

   Instantiates a new FormulationProjection to project the given scalar function
   **

   @fn FormulationProjection::FormulationProjection(const GroupOfElement&,const FunctionSpace&,fullVector<scalar>(*g)(fullVector<double>& xyz))
   @param domain The domain of this Formulation
   @param fs A FunctionSpace for both unknown and test field
   @param g The @em vector function to project

   Instantiates a new FormulationProjection to project the given vector function
   **

   @fn FormulationProjection::~FormulationProjection
   Deletes the this FormulationProjection
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FormulationProjectionInclusion.h"

#endif
