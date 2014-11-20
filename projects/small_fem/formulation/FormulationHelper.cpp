#include <set>
#include "Exception.h"
#include "FormulationHelper.h"

using namespace std;

FormulationHelper::FormulationHelper(void){
}

FormulationHelper::~FormulationHelper(void){
}

void FormulationHelper::initDofMap(const FunctionSpace& fs,
                                   const GroupOfElement& goe,
                                   map<Dof, Complex>& data){
  // Get Keys from goe //
  set<Dof> dSet;
  fs.getKeys(goe, dSet);

  // Insert Dofs in data //
  set<Dof>::iterator it  = dSet.begin();
  set<Dof>::iterator end = dSet.end();

  for(; it != end; it++)
    data.insert(pair<Dof, Complex>(*it, 0));
}

void FormulationHelper::initDofMap(const vector<const FunctionSpace*>& fs,
                                   const GroupOfElement& goe,
                                   vector<map<Dof, Complex> >& data){

  const size_t size = getSize(fs, data);
  for(size_t i = 0; i < size; i++)
    initDofMap(*fs[i], goe, data[i]);
}

void FormulationHelper::initDofMap(const vector<const FunctionSpaceScalar*>& fs,
                                   const GroupOfElement& goe,
                                   vector<map<Dof, Complex> >& data){

  const size_t size = getSize(fs, data);
  for(size_t i = 0; i < size; i++)
    initDofMap(*fs[i], goe, data[i]);
}

void FormulationHelper::initDofMap(const vector<const FunctionSpaceVector*>& fs,
                                   const GroupOfElement& goe,
                                   vector<map<Dof, Complex> >& data){

  const size_t size = getSize(fs, data);
  for(size_t i = 0; i < size; i++)
    initDofMap(*fs[i], goe, data[i]);
}
