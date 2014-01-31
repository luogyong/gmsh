#include <algorithm>
#include <sstream>
#include <list>

#include "GroupOfElement.h"

using namespace std;

const size_t GroupOfElement::nGeoType = 9;

GroupOfElement::GroupOfElement(std::list<const MElement*>& element,
                               const Mesh& mesh){
  // Size //
  const size_t nElement = element.size();

  // Get Elements //
  this->mesh = &mesh;
  this->element.assign(element.begin(), element.end());

  // Elements Type //
  typeStat.resize(nGeoType, 0);

  for(size_t i = 0; i < nElement; i++)
    typeStat[this->element[i]->getType()]++;

  // Sort Elements //
  sort(this->element.begin(), this->element.end(), sortPredicate);

  // Get Orientation Stats //
  const MElement* elementI;

  orientationStat.resize(nGeoType);

  for(size_t i = 0; i < nGeoType; i++)
    if(typeStat[i])
      orientationStat[i].resize(ReferenceSpaceManager::getNOrientation(i), 0);
    else
      orientationStat[i].clear();

  for(size_t i = 0; i < nElement; i++){
    elementI = this->element[i];
    orientationStat[elementI->getType()]
                   [ReferenceSpaceManager::getOrientation(*elementI)]++;
  }
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

std::pair<bool, size_t> GroupOfElement::isUniform(void) const{
  size_t type  = (size_t)(-1);
  bool uniform = true;

  for(size_t i = 0; i < nGeoType && uniform; i++){
    if((type == (size_t)(-1)) && (typeStat[i] != 0))
      type = i;

    else if((type != (size_t)(-1)) && (typeStat[i] != 0))
      uniform = false;
  }

  return pair<bool, size_t>(uniform, type);
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

  return stream.str();
}
