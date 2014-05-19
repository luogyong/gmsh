#include "NodeSolution.h"
#include "Exception.h"
#include "SmallFem.h"

using namespace std;

template<>
void NodeSolution<Complex>::
addNodeValue(size_t step,
             double time,
             const Mesh& mesh,
             std::map<const MVertex*, Complex>& data){

  // GModel //
  GModel& model = mesh.getModel();

  // Map with (Vertex Id, Node Value) //
  map<int, vector<double> > gmshDataReal;
  map<int, vector<double> > gmshDataImag;

  // Scalar Field //
  const size_t nComp = 1;

  // Populate gmshData //
  map<const MVertex*, Complex>::iterator it  = data.begin();
  map<const MVertex*, Complex>::iterator end = data.end();

  vector<double> tmpReal(nComp);
  vector<double> tmpImag(nComp);

  for(; it != end; it++){
    tmpReal[0] = it->second.real();
    tmpImag[0] = it->second.imag();

    gmshDataReal.insert
      (pair<int, vector<double> >(it->first->getNum(), tmpReal));
    gmshDataImag.insert
      (pair<int, vector<double> >(it->first->getNum(), tmpImag));
  }

  // Add map to PView //
  pView->addData(&model, gmshDataReal, 2 * step + 0, time, 0, nComp);
  pView->addData(&model, gmshDataImag, 2 * step + 1, time, 0, nComp);
}

template<>
void NodeSolution<Complex>::
addNodeValue(size_t step,
             double time,
             const Mesh& mesh,
             std::map<const MVertex*, std::vector<Complex> >& data){

  // GModel //
  GModel& model = mesh.getModel();

  // Map with (Vertex Id, Node Value) //
  map<int, vector<double> > gmshDataReal;
  map<int, vector<double> > gmshDataImag;

  // Vector Field //
  const size_t nComp = 3;

  // Populate gmshData //
  map<const MVertex*, vector<Complex> >::iterator it  = data.begin();
  map<const MVertex*, vector<Complex> >::iterator end = data.end();

  vector<double> tmpReal(nComp);
  vector<double> tmpImag(nComp);

  for(; it != end; it++){
    tmpReal[0] = std::abs(it->second.at(0));//.real();
    tmpReal[1] = std::abs(it->second.at(1));//.real();
    tmpReal[2] = std::abs(it->second.at(2));//.real();

    tmpImag[0] = std::arg(it->second.at(0));//.imag();
    tmpImag[1] = std::arg(it->second.at(1));//.imag();
    tmpImag[2] = std::arg(it->second.at(2));//.imag();

    gmshDataReal.insert
      (pair<int, vector<double> >(it->first->getNum(), tmpReal));
    gmshDataImag.insert
      (pair<int, vector<double> >(it->first->getNum(), tmpImag));
  }

  // Add map to PView //
  pView->addData(&model, gmshDataReal, 2 * step + 0, time, 0, nComp);
  pView->addData(&model, gmshDataImag, 2 * step + 1, time, 0, nComp);
}

template<>
void NodeSolution<double>::
addNodeValue(size_t step,
             double time,
             const Mesh& mesh,
             std::map<const MVertex*, double>& data){

  // GModel //
  GModel& model = mesh.getModel();

  // Map with (Vertex Id, Node Value) //
  map<int, vector<double> > gmshData;

  // Scalar Field //
  const size_t nComp = 1;

  // Populate gmshData //
  map<const MVertex*, double>::iterator it  = data.begin();
  map<const MVertex*, double>::iterator end = data.end();

  vector<double> tmp(nComp);

  for(; it != end; it++){
    tmp[0] = it->second;

    gmshData.insert(pair<int, vector<double> >(it->first->getNum(), tmp));
  }

  // Add map to PView //
  pView->addData(&model, gmshData,step, time, 0, nComp);
}

template<>
void NodeSolution<double>::
addNodeValue(size_t step,
             double time,
             const Mesh& mesh,
             std::map<const MVertex*, std::vector<double> >& data){

  // GModel //
  GModel& model = mesh.getModel();

  // Map with (Vertex Id, Node Value) //
  map<int, vector<double> > gmshData;

  // Vector Field //
  const size_t nComp = 3;

  // Populate gmshData //
  map<const MVertex*, vector<double> >::iterator it  = data.begin();
  map<const MVertex*, vector<double> >::iterator end = data.end();

  vector<double> tmp(nComp);

  for(; it != end; it++){
    tmp[0] = it->second.at(0);
    tmp[1] = it->second.at(1);
    tmp[2] = it->second.at(2);

    gmshData.insert(pair<int, vector<double> >(it->first->getNum(), tmp));
  }

  // Add map to PView //
  pView->addData(&model, gmshData, step, time, 0, nComp);
}
