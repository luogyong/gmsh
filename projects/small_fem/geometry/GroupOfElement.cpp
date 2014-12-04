#include <algorithm>
#include <sstream>
#include <list>

#include "GroupOfElement.h"

using namespace std;

const size_t GroupOfElement::nGeoType = 9;
      size_t GroupOfElement::nxtId    = 0;

GroupOfElement::GroupOfElement(std::list<const MElement*>& element,
                               const Mesh& mesh){
  // Get Elements //
  this->mesh = &mesh;
  this->element.assign(element.begin(), element.end());

  // Id //
  this->id = this->nxtId;
  this->nxtId++;

  // Init //
  init();
}

GroupOfElement::GroupOfElement(const Mesh& mesh){
  // Get Elements //
  this->mesh = &mesh;
  this->element.clear();

  // Id //
  this->id = this->nxtId;
  this->nxtId++;

  // Init //
  init();
}

void GroupOfElement::init(void){
  // Size //
  const size_t nElement = element.size();

  // Clear //
  orientationStat.clear();
  typeStat.clear();

  // Elements Type //
  typeStat.resize(nGeoType, 0);

  for(size_t i = 0; i < nElement; i++)
    typeStat[element[i]->getType()]++;

  // Sort Elements //
  sort(element.begin(), element.end(), sortPredicate);

  // Get Orientation Stats //
  orientationStat.resize(nGeoType);

  for(size_t i = 0; i < nGeoType; i++)
    if(typeStat[i])
      orientationStat[i].resize(ReferenceSpaceManager::getNOrientation(i), 0);
    else
      orientationStat[i].clear();

  for(size_t i = 0; i < nElement; i++)
    orientationStat[element[i]->getType()]
                   [ReferenceSpaceManager::getOrientation(*element[i])]++;
}

GroupOfElement::~GroupOfElement(void){
}

void GroupOfElement::add(const GroupOfElement& other){
  // Save this GoE Elements and clear //
  vector<const MElement*> old(element);
  element.clear();

  // Resize and populate //
  const size_t oldSize = old.size();
  const size_t newSize = old.size() + other.element.size();
  element.resize(newSize);

  for(size_t i = 0; i < oldSize; i++)
    element[i] = old[i];

  for(size_t i = oldSize; i < newSize; i++)
    element[i] = other.element[i - oldSize];

  // Init //
  init();
}

bool GroupOfElement::isMember(const MElement& element) const{
  const size_t size = this->element.size();

  for(size_t i = 0; i < size; i++)
    if(element.getNum() == this->element[i]->getNum())
      return true;

  return false;
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
