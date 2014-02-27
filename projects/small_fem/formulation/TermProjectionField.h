#ifndef _TERMPROJECTIONFIELD_H_
#define _TERMPROJECTIONFIELD_H_

#include <map>

#include "FunctionSpaceScalar.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermProjectionField
   @brief Term of a projection onto a Field

   Term of a projection onto a Field.
   The function to project may be defined by:
   @li A scalar function
   @li A map of (Dof, value) and a FunctionSpaceScalar
 */

template<typename scalar>
class TermProjectionField: public Term<scalar>{
 private:
  // Type def //
  typedef scalar
    (TermProjectionField<scalar>::*Eval)(const MElement& element,
                                         const fullVector<double>& xyz) const;
 private:
  // Function for fContainer //
  scalar (*f)(fullVector<double>& xyz);

  // Data for interpolation //
  const FunctionSpaceScalar*   fsScalar;
  const std::map<Dof, scalar>* dofValue;

 public:
  TermProjectionField(const GroupOfJacobian& goj,
                      const Basis& basis,
                      const Quadrature& quadrature,
                      scalar (*f)(fullVector<double>& xyz));

  TermProjectionField(const GroupOfJacobian& goj,
                      const Basis& basis,
                      const Quadrature& quadrature,
                      const FunctionSpaceScalar& fs,
                      const std::map<Dof, scalar>& dof);

  virtual ~TermProjectionField(void);

 private:
  void init(const GroupOfJacobian& goj,
            const Basis& basis,
            const Quadrature& quadrature,
            const Eval& evaluator);

  void computeC(const Basis& basis,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                const fullMatrix<double>& gC,
                const Eval& evaluator,
                fullMatrix<scalar>**& bM);

  scalar  fContainer(const MElement& elem, const fullVector<double>& xyz) const;
  scalar interpolate(const MElement& elem, const fullVector<double>& xyz) const;
};

/**
   @fn TermProjectionField<scalar>::TermProjectionField(const GroupOfJacobian& goj,const Basis& basis,const Quadrature& quadrature,scalar (*f)(fullVector<double>& xyz))
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param quadrature A Quadrature rule
   @param f A scalar function

   Instanciates a new Projection-Field Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The given Quadrature is used to compute the Term
   @li The Basis function must be pre-evaluated at the given integration points

   The projected function is f (0form)

   @todo Evaluate Basis in Term ?????
   **

   @fn TermProjectionField<scalar>::TermProjectionField(const GroupOfJacobian& goj,const Basis& basis,const Quadrature& quadrature,const FunctionSpaceScalar& fs,const std::map<Dof, scalar>& dof)
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param quadrature A Quadrature rule
   @param fs A FunctionSpaceScalar
   @param dof A map of (Dof, value)

   Instanciates a new Projection-Field Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The given Quadrature is used to compute the Term
   @li The Basis function must be pre-evaluated at the given integration points

   The projected function is defined by the given FunctionSpace
   and the map (Dof, value)

   @todo Evaluate Basis in Term ?????
   **

   @fn TermProjectionField<scalar>::~TermProjectionField
   Deletes this TermProjectionField
*/

/////////////////////
// Inline Function //
/////////////////////
template<typename scalar>
inline scalar TermProjectionField<scalar>::
fContainer(const MElement& element,
           const fullVector<double>& xyz) const{
  return f(const_cast<fullVector<double>&>(xyz));
}

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "TermProjectionFieldInclusion.h"

#endif
