#ifndef _TERMPROJECTIONFIELD_H_
#define _TERMPROJECTIONFIELD_H_

#include "GroupOfJacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermProjectionField
   @brief Term of the Field (in physical space) Field (in reference space) type

   Term of the Field (in physical space) Field (in reference space) type
 */

template<typename scalar>
class TermProjectionField: public Term<scalar>{
 public:
  TermProjectionField(const GroupOfJacobian& goj,
                      const Basis& basis,
                      const fullVector<double>& integrationWeights,
                      const fullMatrix<double>& integrationPoints,
                      scalar (*f)(fullVector<double>& xyz));

  virtual ~TermProjectionField(void);

 private:
  void computeC(const Basis& basis,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                const Basis& basis,
                const fullMatrix<double>& gC,
                scalar (*f)(fullVector<double>& xyz),
                fullMatrix<scalar>**& bM);
};

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "TermProjectionFieldInclusion.h"

#endif
