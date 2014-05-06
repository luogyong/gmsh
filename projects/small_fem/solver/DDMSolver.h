#ifndef _DDMSOLVER_H_
#define _DDMSOLVER_H_

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

class DDMSolver{
 private:
  // MPI //
  int numProcs;
  int myId;

  // DDM //
  const Formulation<Complex>* wave;
  const Formulation<Complex>* sommerfeld;
  DDMContext* context;
  Formulation<Complex>* ddm;
  Formulation<Complex>* upDdm;

  // Dirichlet Domain & FunctionSpace //
  const GroupOfElement* dirichlet;
  const FunctionSpace*  fs;

  // DDM Dof values //
  std::map<Dof, Complex>* ddmG;

  // MPI out //
  std::vector<int>     outEntity;
  std::vector<int>     outType;
  std::vector<Complex> outValue;

  // MPI in //
  std::vector<int>     inEntity;
  std::vector<int>     inType;
  std::vector<Complex> inValue;

  // PETSc //
  Mat A;
  Vec x;
  Vec b;

  std::map<Dof, Complex>* rhs;

  typedef struct{
    int myId;
    DDMContext* DDMctx;
    const GroupOfElement* dirichlet;
    const FunctionSpace* fs;
    const Formulation<Complex>* wave;
    const Formulation<Complex>* sommerfeld;
    Formulation<Complex>* ddm;
    Formulation<Complex>* upDdm;

    std::vector<Complex>* outValue;
    std::vector<Complex>* inValue;

  } FullContext;

  FullContext fullCtx;

 public:
  DDMSolver(const Formulation<Complex>& wave,
            const Formulation<Complex>& sommerfeld,
            const GroupOfElement& dirichlet,
            DDMContext& context,
            Formulation<Complex>& ddm,
            Formulation<Complex>& update,
            std::map<Dof, Complex>& rhs);

  ~DDMSolver(void);

  void solve(int nStep);

  void getSolution(std::map<Dof, Complex>& ddm);

 private:
  static void serialize(const std::map<Dof, Complex>& data,
                        std::vector<Complex>& value);

  static void serialize(const Vec& x,
                        std::vector<Complex>& value);

  static void unserialize(std::map<Dof, Complex>& data,
                          const std::vector<Complex>& value);

  static void unserialize(Vec& x,
                          const std::vector<Complex>& value);

  static void exchange(int myId,
                       std::vector<Complex>& outValue,
                       std::vector<Complex>& inValue);

  void displaySolution(const System<Complex>& system,
                       const FunctionSpace& fs,
                       const GroupOfElement& goe,
                       int id, int step, std::string name);

  static PetscErrorCode matMult(Mat A, Vec x, Vec y);

  static void setVecFromDof(Vec& v, std::map<Dof, Complex>& dof);
  static void setDofFromVec(Vec& v, std::map<Dof, Complex>& dof);

  static Complex fZero(fullVector<double>& xyz);

  static void see(Vec& v);
  static void see(const std::map<Dof, Complex>& dof, int id, size_t step, std::string name);
};

#endif
