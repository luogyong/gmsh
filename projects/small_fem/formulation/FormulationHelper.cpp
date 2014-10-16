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
  set<Dof> dSet;
  fs.getKeys(goe, dSet);

  set<Dof>::iterator it  = dSet.begin();
  set<Dof>::iterator end = dSet.end();

  for(; it != end; it++)
    data.insert(pair<Dof, Complex>(*it, 0));
}

void FormulationHelper::initDofMap(const vector<const FunctionSpace*>& fs,
                                   const GroupOfElement& goe,
                                   vector<map<Dof, Complex> >& data){
  const size_t size = data.size();

  set<Dof> dSet;
  set<Dof>::iterator end;
  set<Dof>::iterator it;

  if(size != fs.size())
    throw Exception("FormulationHelper::initDofMap: %s",
                    "vector must have the same size");

  for(size_t i = 0; i < size; i++){
    dSet.clear();
    fs[i]->getKeys(goe, dSet);

    end = dSet.end();

    for(it = dSet.begin(); it != end; it++)
      data[i].insert(pair<Dof, Complex>(*it, 0));
  }
}

void FormulationHelper::initDofMap(const vector<const FunctionSpaceScalar*>& fs,
                                   const GroupOfElement& goe,
                                   vector<map<Dof, Complex> >& data){
  const size_t size = data.size();

  set<Dof> dSet;
  set<Dof>::iterator end;
  set<Dof>::iterator it;

  if(size != fs.size())
    throw Exception("FormulationHelper::initDofMap: %s",
                    "vector must have the same size");

  for(size_t i = 0; i < size; i++){
    dSet.clear();
    fs[i]->getKeys(goe, dSet);

    end = dSet.end();

    for(it = dSet.begin(); it != end; it++)
      data[i].insert(pair<Dof, Complex>(*it, 0));
  }
}

void FormulationHelper::initDofMap(const vector<const FunctionSpaceVector*>& fs,
                                   const GroupOfElement& goe,
                                   vector<map<Dof, Complex> >& data){
  const size_t size = data.size();

  set<Dof> dSet;
  set<Dof>::iterator end;
  set<Dof>::iterator it;

  if(size != fs.size())
    throw Exception("FormulationHelper::initDofMap: %s",
                    "vector must have the same size");

  for(size_t i = 0; i < size; i++){
    dSet.clear();
    fs[i]->getKeys(goe, dSet);

    end = dSet.end();

    for(it = dSet.begin(); it != end; it++)
      data[i].insert(pair<Dof, Complex>(*it, 0));
  }
}
