#ifndef _INTERPOLATOR_H_
#define _INTERPOLATOR_H_

#include <map>

#include "GroupOfElement.h"
#include "FunctionSpace.h"
#include "fullMatrix.h"

/**
   @class Interpolator
   @brief Interpolating method for FEM Solutions

   This class allows the interpolation of a map of (Dof, value)
   associated to a FunctionSpace

   @todo
   NEED TO BE MERGED WITH FEMSOLUTION: NEED HO PVIEW WITH ARBITRARY BASIS
 */

template <typename scalar>
class Interpolator{
 public:
   Interpolator(void);
  ~Interpolator(void);

  static void interpolate(const GroupOfElement& goe,
                          const FunctionSpace& fs,
                          const std::map<Dof, scalar>& coef,
                          const fullMatrix<double>& point,
                          fullMatrix<scalar>& values);

  static void interpolate(const GroupOfElement& goe,
                          const GroupOfElement& point,
                          const FunctionSpace& fs,
                          const std::map<Dof, scalar>& coef,
                          std::map<const MVertex*,
                                   std::vector<scalar> >& values);
 private:
  static void interpolate(const MElement& element,
                          const FunctionSpace& fs,
                          const std::vector<scalar>& coef,
                          const fullVector<double>& xyz,
                          fullVector<scalar>& value);
};


/**
   @fn Interpolator::Interpolator
   Instanciate a new Interpolator
   (this is not required, since Interpolator has only one class method)
   **

   @fn Interpolator::~Interpolator
   Deletes this Interpolator
   (this is not required, since Interpolator has only one class method)
   **

   @fn Interpolator::interpolate(const GroupOfElement&,const FunctionSpace&,const std::map<Dof, scalar>&,const fullMatrix<double>&,fullMatrix<scalar>&)
   @param goe A GroupOfElement on which the solution will be interpolated
   @param fs A FunctionSpace
   @param coef A map of (Dof, value) to be interpolated with the FunvtionSpace
   @param point A set of point coordinates (3D):
   one point per row and one coordinate per column
   @param values A matrix

   Interpolate the given parameters on the given set of point.

   The interpolated values are stored in values:
   @li Each row is point
   @li Each column is a coordinate of the solution
   (1 for scalar problems and 3 for vectorial ones)
 */

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "InterpolatorInclusion.h"

#endif
