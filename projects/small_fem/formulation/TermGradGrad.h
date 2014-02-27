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
  typedef const fullMatrix<double>& (Basis::*BFunction)(size_t s) const;

 private:
  std::vector<fullMatrix<scalar> > TJac; //Pre-evaluated tensor * jacobian

 public:
  TermGradGrad(const GroupOfJacobian& goj,
               const Basis& basis,
               const Quadrature& quadrature);

  TermGradGrad(const GroupOfJacobian& goj,
               const Basis& basis,
               const Quadrature& quadrature,
               void  (*f)(fullVector<double>& xyz, fullMatrix<scalar>& T));

  virtual ~TermGradGrad(void);

 private:
  // Init
  void init(const GroupOfJacobian& goj,
            const Basis& basis,
            const Quadrature& quadrature);

  // Matrices
  void computeC(const Basis& basis,
                const BFunction& getFunction,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                size_t nG,
                fullMatrix<scalar>**& bM);

  // Tensor Pre Evaluator
  void preEvalDummy(const GroupOfJacobian& goj,
                    const Quadrature& quadrature);

  void preEvalT(const GroupOfJacobian& goj,
                const Quadrature& quadrature,
                void  (*f)(fullVector<double>& xyz, fullMatrix<scalar>& T));
};

/**
   @fn TermGradGrad::TermGradGrad(const GroupOfJacobian& goj,const Basis& basis,const Quadrature& quadrature)
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

   @fn TermGradGrad::TermGradGrad(const GroupOfJacobian& goj,const Basis& basis,const Quadrature& quadrature,void (*f)(fullVector<double>& xyz, fullMatrix<scalar>& T))
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param quadrature A Quadrature rule
   @param f A multiplicative tensorial function

   Instanciates a new Grad-Grad Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The given Quadrature is used to compute the Term
   @li The Basis function must be pre-evaluated at the integration points
   @li The given function is a multiplicative tensor for the whole Term

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
