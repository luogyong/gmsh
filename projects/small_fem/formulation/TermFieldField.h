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
 public:
  TermFieldField(const GroupOfJacobian& goj,
                 const Basis& basis,
                 const Quadrature& quadrature);

  virtual ~TermFieldField(void);

 private:
  void computeC(const Basis& basis,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                size_t nG,
                fullMatrix<scalar>**& bM);
};

/**
   @fn TermFieldField::TermFieldField
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
