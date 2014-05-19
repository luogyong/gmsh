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
  bool             isModulusPhase;
  PViewDataGModel* pView;

 public:
   FEMSolution(void);
  ~FEMSolution(void);

  void setModulusPhase(void);
  void setRealImaginary(void);

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
   Instanciates a new FEMSolution which is empty, and calls setRealImaginary
   **

   @fn FEMSolution::~FEMSolution
   Deletes this FEMSolution
   **

   @fn FEMSolution::setModulusPhase
   This FEMSolution will use complex number in a modulus - phase fashion
   (no impact on real FEMSolution)
   **

   @fn FEMSolution::setRealImaginary
   This FEMSolution will use complex number in a real - imaginary part fashion
   (no impact on real FEMSolution)
   **

   @fn FEMSolution::clear
   This FEMSolution is now empty
   **

   @fn FEMSolution::addCoefficients
   @param step An integer value
   @param time A real value
   @param goe A GroupOfElement
   @param fs A FunctionSpace
   @param coef A map associating some (or all) the Dof%s
   of the given FunctionSpace to a value

   Computes the FEM solution on the elements of the given GroupOfElement
   with the given FunctionSpace at map (Dof, value)

   This solution is added to this FEMSolution with the given step and time.
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
