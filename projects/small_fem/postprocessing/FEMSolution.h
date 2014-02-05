#ifndef _FEMSOLUTION_H_
#define _FEMSOLUTION_H_

#include <string>
#include <map>

#include "FunctionSpace.h"
#include "BasisLagrange.h"
#include "PViewDataGModel.h"

/**
   @class FEMSolution
   @brief Solution of a finite element problem

   This class represents a finite element problem solution.

   A FEMSolution is defined thanks to:
   @li A set of elements (MElement)
   @li An interpolation matrix (for Lagrange)
   @li A set of interpolation coefficient (for Lagrange) for each element

   A FEMSolution can write a post-processing map in the
   <a href="http://www.geuz.org/gmsh">gmsh</a> .msh file format.

   The same instance of FEMSolution may handle multiple sets of coefficients.
   Each set is associated to an integer called a 'step'.
   Each step is also associated to real value called 'time'.
   If two sets have the same step, the last one override to first one.
   Differents steps may have the same time.
 */

template<typename scalar>
class FEMSolution{
 private:
  PViewDataGModel* pView;

 public:
   FEMSolution(void);
  ~FEMSolution(void);

  void clear(void);
  void addCoefficients(size_t step,
                       double time,
                       const GroupOfElement& goe,
                       const FunctionSpace& fs,
                       const std::map<Dof, scalar>& coef);

  void write(std::string fileName) const;

 private:
  void toLagrange(const MElement& element,
                  const std::vector<BasisLagrange*>& lagrange,
                  const std::vector<scalar>& fsCoef,
                  const FunctionSpace& fs,
                  std::vector<scalar>& lCoef);

  void toPView(GModel& model, std::map<int, std::vector<scalar> >& data,
               size_t step, double time, int partition, int nComp);
};


/**
   @fn FEMSolution::FEMSolution
   Instanciates a new FEMSolution which is empty
   (no elements, interpolation matrix and coefficients)
   **

   @fn FEMSolution::~FEMSolution
   Deletes this FEMSolution
   **

   @fn FEMSolution::clear
   This FEMSolution is now empty
   **

   @fn FEMSolution::addCoefficients(size_t,double,const FunctionSpace&,const DofManager&,const fullVector<scalar>&)
   @param step An integer value
   @param time A real value
   @param fs A FunctionSpace
   @param dofM A DofManager
   @param coef A set of coefficient, whose indexes are related to the DofManager

   Adds the given set of coefficient to this FEMSolution
   with the given step and time.

   The set of elements is extracted from the FunctionSpace::getSupport().
   A Lagrange basis (BasisLagrange) is constructed for this set of elements.
   Its interpolation matrices will be used in this FEMSolution.
   The function space and the DofManager are used to compute a set
   of coefficient in the Lagrange basis for each elements.
   **

   @fn FEMSolution::write
   @param fileName A file name (without extension)

   Writes this FEMSolution in
   <a href="http://www.geuz.org/gmsh">gmsh</a> .msh file format
   into the given file
 */

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "FEMSolutionInclusion.h"

#endif
