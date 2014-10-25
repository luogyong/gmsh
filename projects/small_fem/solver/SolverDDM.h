#ifndef _SOLVERDDM_H_
#define _SOLVERDDM_H_

#include <string>
#include <vector>
#include <map>

#include "petscmat.h"
#include "petscvec.h"

#include "FunctionSpace.h"
#include "Formulation.h"
#include "DDMContext.h"
#include "System.h"

#include "SmallFem.h"

class SolverDDM{
 private:
  // MPI //
  int nProcs;
  int myProc;

  // DDM //
  const Formulation<Complex>* wave;
  const Formulation<Complex>* sommerfeld;
  Formulation<Complex>*       ddm;
  Formulation<Complex>*       upDdm;
  DDMContext*                 context;

  // Volume & Update System //
  System<Complex>* volume;
  System<Complex>* update;
  bool             once;

  // Dirichlet Domain & FunctionSpace //
  const GroupOfElement* dirichlet;
  const FunctionSpace*  fs;

  // DDM Dof values //
  std::map<Dof, Complex>* ddmG;
  std::vector<int>        neighbour;
  int                     neighbourOne;
  int                     neighbourTwo;

  // MPI out //
  std::vector<int>     outEntityOne;
  std::vector<int>     outTypeOne;
  std::vector<Complex> outValueOne;

  std::vector<int>     outEntityTwo;
  std::vector<int>     outTypeTwo;
  std::vector<Complex> outValueTwo;

  // MPI in //
  std::vector<int>     inEntityOne;
  std::vector<int>     inTypeOne;
  std::vector<Complex> inValueOne;

  std::vector<int>     inEntityTwo;
  std::vector<int>     inTypeTwo;
  std::vector<Complex> inValueTwo;

  // PETSc //
  Mat A;
  Vec x;
  Vec b;

  std::map<Dof, Complex>* rhs;

 public:
  SolverDDM(const Formulation<Complex>& wave,
            const Formulation<Complex>& sommerfeld,
            const GroupOfElement& dirichlet,
            DDMContext& context,
            Formulation<Complex>& ddm,
            Formulation<Complex>& update,
            std::map<Dof, Complex>& rhs);

  ~SolverDDM(void);

  void solve(int nStep);
  void constructIterationMatrix(std::string name, std::string filename);

  void getSolution(std::map<Dof, Complex>& ddm);

 private:
  void buildNeighbourhood(const std::map<Dof, Complex>& local);
  std::pair<size_t, size_t> splitSize(void);

  void serialize(std::map<Dof, Complex>& data, int neighbour,
                 std::vector<int>& entity, std::vector<int>& type,
                 std::vector<Complex>& value);
  void unserialize(std::map<Dof, Complex>& data,
                   std::vector<int>& entity, std::vector<int>& type,
                   std::vector<Complex>& value);
  template<typename T>
    void exchange(int target, std::vector<T>& out, std::vector<T>& in,
                  MPI_Request* send, MPI_Request* recv);

  void exchange(std::map<Dof, Complex>& data);

 private:
  static void setVecFromDof(Vec& v, std::map<Dof, Complex>& dof);
  static void setDofFromVec(Vec& v, std::map<Dof, Complex>& dof);

  static Complex             fZeroScal(fullVector<double>& xyz);
  static fullVector<Complex> fZeroVect(fullVector<double>& xyz);

  static PetscErrorCode matMult(Mat A, Vec x, Vec y);
};

#endif
