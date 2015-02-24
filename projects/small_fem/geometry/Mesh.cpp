#include <vector>
#include <sstream>

#include "GroupOfElement.h"
#include "GeoExtractor.h"
#include "Exception.h"
#include "Mesh.h"

using namespace std;

Mesh::Mesh(const std::string fileName){
  // Clear Maps //
  physical.clear();
  element.clear();
  vertex.clear();
  edge.clear();
  face.clear();

  // New Mode //
  model = new GModel("SmallFEM");

  // Read Mesh //
  if(!model->readMSH(fileName))
    throw Exception("Mesh: cannot open file: %s", fileName.c_str());

  // Get Entity //
  vector<GEntity*> entity;
  model->getEntities(entity);

  // Extract Element, Nodes, Edges and Faces //
  GeoExtractor::elementExtract(entity,  element, physical);
  GeoExtractor::vertexExtract(element, vertex);
  GeoExtractor::edgeExtract(element, edge);
  GeoExtractor::faceExtract(element, face);

  // Number Geometry //
  number();
}

Mesh::~Mesh(void){
  // Delete Model //
  delete model;
}

size_t Mesh::getGlobalId(const MElement& element) const{
  MapElement::const_iterator it = this->element.find(&element);

  if(it == this->element.end())
    throw Exception("Mesh::getGlobalId(): element not found");

  return it->second;
}

size_t Mesh::getGlobalId(const MVertex& vertex) const{
  MapVertex::const_iterator it = this->vertex.find(&vertex);

  if(it == this->vertex.end())
    throw Exception("Mesh::getGlobalId(): vertex not found");

  return it->second;
}

size_t Mesh::getGlobalId(const MEdge& edge) const{
  // Get Edge Vertices //
  vector<int> vertex(2);
  vertex[0] = edge.getVertex(0)->getNum();
  vertex[1] = edge.getVertex(1)->getNum();

  // Look for Edge //
  MapEntity::const_iterator it = this->edge.find(vertex);

  if(it == this->edge.end())
    throw Exception("Mesh::getGlobalId(): edge not found");

  return it->second;
}

size_t Mesh::getGlobalId(const MFace& face) const{
  // Get Face Vertices //
  const int  nVertex = face.getNumVertices();
  vector<int> vertex(nVertex);

  for(int i = 0; i < nVertex; i++)
    vertex[i] = face.getVertex(i)->getNum();

  // Look for Face //
  MapEntity::const_iterator it = this->face.find(vertex);

  if(it == this->face.end())
    throw Exception("Mesh::getGlobalId(): face not found");

  return it->second;
}

void Mesh::number(void){
  // Get Iterators //
  const MapElement::iterator endEl = element.end();
  const MapVertex::iterator  endV  = vertex.end();
  const MapEntity::iterator  endEd = edge.end();
  const MapEntity::iterator  endF  = face.end();

  MapElement::iterator itEl = element.begin();
  MapVertex::iterator  itV  = vertex.begin();
  MapEntity::iterator  itEd = edge.begin();
  MapEntity::iterator  itF  = face.begin();

  // Number Vertices //
  size_t nextId = 0;

  for(; itV != endV; itV++){
    itV->second = nextId;
    nextId++;
  }

  // Number Edges //
  for(; itEd != endEd; itEd++){
    itEd->second = nextId;
    nextId++;
  }

  // Number Faces //
  for(; itF != endF; itF++){
    itF->second = nextId;
    nextId++;
  }

  // Number Elements //
  for(; itEl != endEl; itEl++){
    itEl->second = nextId;
    nextId++;
  }
}

GroupOfElement Mesh::getFromPhysical(int physicalId) const{
  pair<multimap<int, const MElement*>::const_iterator,
       multimap<int, const MElement*>::const_iterator>
    p = physical.equal_range(physicalId);

  list<const MElement*> lst;

  for(; p.first != p.second; p.first++)
    lst.push_back(p.first->second);

  return GroupOfElement(lst, *this);
}

GroupOfElement Mesh::getFromPhysical(int physicalId, int partitionId) const{
  pair<multimap<int, const MElement*>::const_iterator,
       multimap<int, const MElement*>::const_iterator>
    p = physical.equal_range(physicalId);

  list<const MElement*> lst;

  for(; p.first != p.second; p.first++)
    if(p.first->second->getPartition() == partitionId)
      lst.push_back(p.first->second);

  return GroupOfElement(lst, *this);
}

string Mesh::toString(void) const{
  // Iterators //
  const MapElement::const_iterator endEl = element.end();
  const MapVertex::const_iterator  endV  = vertex.end();
  const MapEntity::const_iterator  endEd = edge.end();
  const MapEntity::const_iterator  endF  = face.end();

  MapElement::const_iterator itEl = element.begin();
  MapVertex::const_iterator  itV  = vertex.begin();
  MapEntity::const_iterator  itEd = edge.begin();
  MapEntity::const_iterator  itF  = face.begin();

  // Stream //
  stringstream stream;

  // Header //
  stream << "***********************************************"
         << endl
         << "*                     Mesh                    *"
         << endl
         << "***********************************************"
         << endl;


  // Vertices //
  stream << "*                                             *"
         << endl
         << "* This mesh contains the following Vertex:    *"
         << endl;

  for(; itV != endV; itV++)
    stream << "*   -- Vertex "
           << showpos
           << getGlobalId(*itV->first)
           << ": {"
           << itV->first->x()
           << ", "
           << itV->first->y()
           << ", "
           << itV->first->z()
           << "}"
           << endl;

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;


  // Edges //
  stream << "*                                             *"
         << endl
         << "* This mesh contains the following Edges:     *"
         << endl;

  for(; itEd != endEd; itEd++)
    stream << "*   -- Edge "
           << itEd->second << ": "
           << "[" << itEd->first[0] << ", " << itEd->first[1] << "]"
           << endl;

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;

  // Faces //
  stream << "*                                             *"
         << endl
         << "* This mesh contains the following Faces:     *"
         << endl;

  for(; itF != endF; itF++){
    stream << "*   -- Face "
           << itF->second << ": "
           << "["
           << itF->first[0] << ", "
           << itF->first[1] << ", "
           << itF->first[2];

    if(itF->first.size() == 4)
      stream << ", " << itF->first[3];

    stream << "]" << endl;
  }

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;


  // Elements //
  stream << "*                                             *"
         << endl
         << "* This mesh contains the following Elements:  *"
         << endl;

  for(; itEl != endEl; itEl++){
    int nVertex =
      const_cast<MElement*>(itEl->first)->getNumPrimaryVertices();

    int nVertexMinus = nVertex - 1;

    stream << "*   -- Element "
           << getGlobalId(*itEl->first)
           << ": [";

    for(int i = 0; i < nVertexMinus; i++)
      stream << getGlobalId(*const_cast<MElement*>(itEl->first)->getVertex(i))
             << ", ";

    stream <<
      getGlobalId(*const_cast<MElement*>(itEl->first)->getVertex(nVertexMinus))
           << "]"
           << endl;
  }

  stream << "*                                             *"
         << endl
         << "***********************************************"
         << endl;

  // Retrun //
  return stream.str();
}
