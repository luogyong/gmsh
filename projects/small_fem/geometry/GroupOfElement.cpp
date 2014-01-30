#include <algorithm>
#include <sstream>
#include <list>

#include "GroupOfElement.h"

using namespace std;

const size_t GroupOfElement::nGeoType = 9;

GroupOfElement::GroupOfElement(std::list<const MElement*>& elementList,
                               const Mesh& mesh){
  // Size //
  const size_t nElement = elementList.size();

  // Get Elements //
  this->mesh    = &mesh;
  this->element.assign(elementList.begin(), elementList.end());

  // Elements Type //
  typeStat.resize(nGeoType, 0);

  for(size_t i = 0; i < nElement; i++)
    typeStat[element[i]->getType()]++;

  // Sort Elements //
  sort(element.begin(), element.end(), sortPredicate);

  // Get Orientation Stats //
  // Get some Data
  const size_t nOrient  =
    ReferenceSpaceManager::getNOrientation(element[0]->getType());

  // Init
  orientationStat.resize(nOrient);

  for(size_t i = 0; i < nOrient; i++)
    orientationStat[i] = 0;

  // Compute
  for(size_t i = 0; i < nElement; i++)
    orientationStat[ReferenceSpaceManager::getOrientation(*element[i])]++;
}

GroupOfElement::~GroupOfElement(void){
}

void GroupOfElement::
getAllVertex(std::set<const MVertex*, VertexComparator>& vertex) const{
  const size_t nElement = element.size();

  // Loop On element
  for(size_t i = 0; i < nElement; i++){
    // Get Vertex
    const size_t nVertex = element[i]->getNumVertices();

    for(size_t j = 0; j < nVertex; j++)
      vertex.insert(element[i]->getVertex(j));
  }
}

void GroupOfElement::getAllVertexCoordinate(fullMatrix<double>& coord) const{
  // Get Vertex
  set<const MVertex*, VertexComparator> vertex;
  getAllVertex(vertex);

  // Get Coordinates
  set<const MVertex*, VertexComparator>::iterator it  = vertex.begin();
  set<const MVertex*, VertexComparator>::iterator end = vertex.end();

  coord.resize(vertex.size(), 3);

  for(size_t i = 0; it != end; it++, i++){
    coord(i, 0) = (*it)->x();
    coord(i, 1) = (*it)->y();
    coord(i, 2) = (*it)->z();
  }
}

string GroupOfElement::toString(void) const{
  stringstream stream;

  stream << "***********************************************"
         << endl
         << "* Group Of Element"
         << endl
         << "***********************************************"
         << endl
         << "*                                             *"
         << endl
         << "* This group contains the following Elements: *"
         << endl;

  for(size_t i = 0; i < element.size(); i++)
    stream << "*   -- Element #"
           << mesh->getGlobalId(*element[i])
           << endl;

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;

  stream << "*                                             *"
         << endl
         << "* This group has the following Orientations:  *"
         << endl;

  for(size_t i = 0; i < orientationStat.size(); i++)
    stream << "*   -- Elements with Orientation " << i << " - "
           << orientationStat[i] << endl;

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;

  return stream.str();
}
