#ifndef _TERMGRADGRAD_H_
#define _TERMGRADGRAD_H_

#include "GroupOfJacobian.h"
#include "Quadrature.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermGradGrad
   @brief Term of the Grad Grad type

   Term of the Grad Grad type
 */

template<typename scalar>
class TermGradGrad: public Term<scalar>{
 private:
  // Type def //
  typedef const fullMatrix<double>& (Basis::*bFunction)(size_t s)const;

 public:
  TermGradGrad(const GroupOfJacobian& goj,
               const Basis& basis,
               const Quadrature& quadrature);

  virtual ~TermGradGrad(void);

 private:
  // Matrices
  void computeC(const Basis& basis,
                const bFunction& getFunction,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                size_t nG,
                fullMatrix<scalar>**& bM);
};

/**
   @fn TermGradGrad::TermGradGrad(const GroupOfJacobian&,const Basis&,const Quadrature&)
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param quadrature A Quadrature rule

   Instanciates a new Grad-Grad Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The given Quadrature is used to compute the Term
   @li The Basis function must be pre-evaluated at the integration points

   @todo Evaluate Basis in Term ?????
   **

   @fn TermGradGrad::~TermGradGrad
   Deletes this TermGradGrad
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "TermGradGradInclusion.h"

#endif
