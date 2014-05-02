#ifndef _DDMSOLVER_H_
#define _DDMSOLVER_H_

#include <string>
#include <vector>
#include <map>

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

  // FunctionSpace & Domain //
  const FunctionSpace*  fs;
  const GroupOfElement* ddmBorder;
  const GroupOfElement* source;

  // fSource //
  Complex (*fSource)(fullVector<double>& xyz);

  // DDM Dof values //
  std::map<Dof, Complex>* ddmG;

  // Systems //
  System<Complex>* volume;
  System<Complex>* update;

  // MPI out //
  std::vector<int>     outEntity;
  std::vector<int>     outType;
  std::vector<Complex> outValue;

  // MPI in //
  std::vector<int>     inEntity;
  std::vector<int>     inType;
  std::vector<Complex> inValue;

 public:
  DDMSolver(const Formulation<Complex>& wave,
            const Formulation<Complex>& sommerfeld,
            DDMContext& context,
            Formulation<Complex>& ddm,
            Formulation<Complex>& update,
            std::map<Dof, Complex>& ddmG,
            const GroupOfElement& source,
            Complex (*fSource)(fullVector<double>& xyz));

  ~DDMSolver(void);

  void solve(int nStep);

 private:
  void serialize(const std::map<Dof, Complex>& data,
                 std::vector<int>& entity,
                 std::vector<int>& type,
                 std::vector<Complex>& value);

  void unserialize(std::map<Dof, Complex>& data,
                   const std::vector<int>& entity,
                   const std::vector<int>& type,
                   const std::vector<Complex>& value);

  void exchange(int myId,
                std::vector<int>& outEntity,
                std::vector<int>& outType,
                std::vector<Complex>& outValue,
                std::vector<int>& inEntity,
                std::vector<int>& inType,
                std::vector<Complex>& inValue);

  void displaySolution(const System<Complex>& system,
                       const FunctionSpace& fs,
                       const GroupOfElement& goe,
                       int id, int step, std::string name);
};

#endif
