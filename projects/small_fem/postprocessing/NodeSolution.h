#ifndef _NODESOLUTION_H_
#define _NODESOLUTION_H_

#include <string>
#include <complex>
#include <map>

#include "Mesh.h"
#include "PViewDataGModel.h"

/**
   @class NodeSolution
   @brief TODO
 */

template<typename scalar>
class NodeSolution{
 private:
  PViewDataGModel* pView;

 public:
   NodeSolution(void);
  ~NodeSolution(void);

  void addNodeValue(size_t step,
                    double time,
                    const Mesh& mesh,
                    std::map<const MVertex*, scalar>& data);

  void addNodeValue(size_t step,
                    double time,
                    const Mesh& mesh,
                    std::map<const MVertex*, std::vector<scalar> >& data);

  void write(std::string fileName) const;
};

//////////////////////////////////////
// Templates Implementations:       //
// Inclusion compilation model      //
//                                  //
// Damn you gcc: we want 'export' ! //
//////////////////////////////////////

#include "NodeSolutionInclusion.h"

#endif
