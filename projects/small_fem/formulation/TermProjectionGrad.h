#ifndef _TERMPROJECTIONGRAD_H_
#define _TERMPROJECTIONGRAD_H_

#include <map>

#include "FunctionSpaceScalar.h"
#include "GroupOfJacobian.h"
#include "Basis.h"
#include "Term.h"

/**
   @class TermProjectionGrad
   @brief Term of the Field (in physical space) Grad (in reference space) type

   Term of the Field (in physical space) Grad (in reference space) type
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
  const FunctionSpaceScalar*   fsScalar;
  const std::map<Dof, scalar>* dofValue;

 public:
  TermProjectionGrad(const GroupOfJacobian& goj,
                     const Basis& basis,
                     const fullVector<double>& integrationWeights,
                     const fullMatrix<double>& integrationPoints,
                     fullVector<scalar> (*f)(fullVector<double>& xyz));

  TermProjectionGrad(const GroupOfJacobian& goj,
                     const Basis& basis,
                     const fullVector<double>& integrationWeights,
                     const fullMatrix<double>& integrationPoints,
                     const FunctionSpaceScalar& fs,
                     const std::map<Dof, scalar>& dof);

  virtual ~TermProjectionGrad(void);

 private:
  void init(const GroupOfJacobian& goj,
            const Basis& basis,
            const fullVector<double>& integrationWeights,
            const fullMatrix<double>& integrationPoints,
            const Eval& evaluator);

  void computeC(const Basis& basis,
                const BFunction& getFunction,
                const fullVector<double>& gW,
                fullMatrix<scalar>**& cM);

  void computeB(const GroupOfJacobian& goj,
                const Basis& basis,
                const fullMatrix<double>& gC,
                const Eval& evaluator,
                fullMatrix<scalar>**& bM);

  fullVector<scalar> fContainer(const MElement& element,
                                const fullVector<double>& xyz) const;

  fullVector<scalar> interpolateGrad(const MElement& element,
                                     const fullVector<double>& xyz) const;
};

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
