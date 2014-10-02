#ifndef _TERMPROJECTIONGRAD_H_
#define _TERMPROJECTIONGRAD_H_

#include <map>

#include "FunctionSpaceScalar.h"
#include "FunctionSpaceVector.h"
#include "GroupOfJacobian.h"
#include "Quadrature.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermProjectionGrad
   @brief Term of a projection onto a Grad

   Term of a projection onto a Grad.
   The function to project may be defined by:
   @li A vectorial function
   @li A map of (Dof, value) and a gradient of a FunctionSpaceScalar
   @li A map of (Dof, value) and a FunctionSpaceVector
 */

template<typename scalar>
class TermProjectionGrad: public Term<scalar>{
 private:
  // Type defs //
  typedef const fullMatrix<double>& (Basis::*BFunction)(size_t s)const;
  typedef fullVector<scalar>
    (TermProjectionGrad<scalar>::*Eval)(const MElement& element,
                                        const fullVector<double>& xyz) const;
 private:
  // Function for fContainer //
  fullVector<scalar> (*f)(fullVector<double>& xyz);

  // Data for interpolation //
  Eval                        evaluator;
  const FunctionSpaceScalar*   fsScalar;
  const FunctionSpaceVector*   fsVector;
  const std::map<Dof, scalar>* dofValue;

 public:
  TermProjectionGrad(const GroupOfJacobian& goj,
                     const Basis& basis,
                     const Quadrature& quadrature,
                     fullVector<scalar> (*f)(fullVector<double>& xyz));

  TermProjectionGrad(const GroupOfJacobian& goj,
                     const Basis& basis,
                     const Quadrature& quadrature,
                     const FunctionSpace& fs,
                     const std::map<Dof, scalar>& dof);

  virtual ~TermProjectionGrad(void);

 private:
  void init(const GroupOfJacobian& goj,
            const Basis& basis,
            const Quadrature& quadrature,
            const Eval& evaluator);

  void computeC(const Basis& basis,
                const BFunction& getFunction,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                const fullMatrix<double>& gC,
                const Eval& evaluator,
                fullMatrix<scalar>**& bM);

  fullVector<scalar> fContainer(const MElement& element,
                                const fullVector<double>& xyz) const;

  fullVector<scalar> interpolateGrad(const MElement& element,
                                     const fullVector<double>& xyz) const;

  fullVector<scalar> interpolate(const MElement& element,
                                 const fullVector<double>& xyz) const;
};

/**
   @fn TermProjectionGrad<scalar>::TermProjectionGrad(const GroupOfJacobian& goj,const Basis& basis,const Quadrature& quadrature,fullVector<scalar> (*f)(fullVector<double>& xyz))
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param quadrature A Quadrature rule
   @param f A vectorial function

   Instanciates a new Projection-Grad Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The given Quadrature is used to compute the Term
   @li The Basis function must be pre-evaluated at the given integration points

   The projected function is f (1form)

   @todo Evaluate Basis in Term ?????
   **

   @fn TermProjectionGrad<scalar>::TermProjectionGrad(const GroupOfJacobian& goj,const Basis& basis,const Quadrature& quadrature,const FunctionSpaceScalar& fs,const std::map<Dof, scalar>& dof)
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param quadrature A Quadrature rule
   @param fs A FunctionSpaceScalar
   @param dof A map of (Dof, value)

   Instanciates a new Projection-Grad Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The given Quadrature is used to compute the Term
   @li The Basis function must be pre-evaluated at the given integration points

   The projected function is defined by the given FunctionSpaceScalar
   and the map (Dof, value).

   Since it a projection onto a Grad Space,
   the Grad of the given function space is used

   @todo Evaluate Basis in Term ?????
   **

   @fn TermProjectionGrad<scalar>::~TermProjectionGrad
   Deletes this TermProjectionGrad
*/

/////////////////////
// Inline Function //
/////////////////////
template<typename scalar>
inline fullVector<scalar> TermProjectionGrad<scalar>::
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

#include "TermProjectionGradInclusion.h"

#endif
