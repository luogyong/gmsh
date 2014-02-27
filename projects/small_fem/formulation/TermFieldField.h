#ifndef _TERMFIELDFIELD_H_
#define _TERMFIELDFIELD_H_

#include "GroupOfJacobian.h"
#include "Quadrature.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermFieldField
   @brief A Term of the Field Field type

   A Term of the Field Field type
 */

template<typename scalar>
class TermFieldField: public Term<scalar>{
 private:
  std::vector<scalar> alpha; //Pre-evaluated function

 public:
  TermFieldField(const GroupOfJacobian& goj,
                 const Basis& basis,
                 const Quadrature& quadrature);

  TermFieldField(const  GroupOfJacobian& goj,
                 const  Basis& basis,
                 const  Quadrature& quadrature,
                 scalar (*f)(fullVector<double>& xyz));

  virtual ~TermFieldField(void);

 private:
  // Init
  void init(const GroupOfJacobian& goj,
            const Basis& basis,
            const Quadrature& quadrature);

  // Matrices
  void computeC(const Basis& basis,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                size_t nG,
                fullMatrix<scalar>**& bM);

  // Tensor Pre Evaluator
  void preEvalDummy(const GroupOfJacobian& goj,
                    const Quadrature& quadrature);

  void preEvalF(const GroupOfJacobian& goj,
                const Quadrature& quadrature,
                scalar (*f)(fullVector<double>& xyz));
};

/**
   @fn TermFieldField::TermFieldField(const GroupOfJacobian& goj,const Basis& basis,const Quadrature& quadrature)
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param quadrature A Quadrature rule

   Instanciates a new Field-Field Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The given Quadrature is used to compute the Term
   @li The Basis function must be pre-evaluated at the integration points

   @todo Evaluate Basis in Term ?????
   **

   @fn TermFieldField::TermFieldField(const GroupOfJacobian& goj,const Basis& basis,const  Quadrature& quadrature,scalar (*f)(fullVector<double>& xyz))
   @param goj A GroupOfJacobian
   @param basis A Basis
   @param quadrature A Quadrature rule
   @param f A multiplicative scalar function

   Instanciates a new Field-Field Term:
   @li The geomtry and the Jacobians are given by the GroupOfJacobian
   @li The Basis functions to use are given by the Basis
   @li The given Quadrature is used to compute the Term
   @li The Basis function must be pre-evaluated at the integration points
   @li The given function is a multiplicative scalar for the whole Term

   @todo Evaluate Basis in Term ?????
   **

   @fn TermFieldField::~TermFieldField
   Deletes this TermFieldField
*/

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "TermFieldFieldInclusion.h"

#endif
