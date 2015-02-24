#include "GeoExtractor.h"

using namespace std;

GeoExtractor::GeoExtractor(void){
}

GeoExtractor::~GeoExtractor(void){
}

void GeoExtractor::elementExtract(const vector<GEntity*>& entity,
                                  MapElement& element,
                                  multimap<int, const MElement*>& physical){
  // Iterator on Entities
  const size_t nEntity = entity.size();

  for(size_t i = 0; i < nEntity; i++){
    // Get Mesh Elements
    const size_t nElement = entity[i]->getNumMeshElements();
    vector<MElement*> myElement(nElement);

    for(size_t j = 0; j < nElement; j++)
      myElement[j] = entity[i]->getMeshElement(j);

    // Get Physical
    const vector<int> myPhysical = entity[i]->getPhysicalEntities();
    const size_t nPhysical       = myPhysical.size();

    // Insert Element
    for(size_t j = 0; j < nElement; j++){
      bool isInsert = element.insert(myElement[j]);

      // If Insertion is a success: insert Physical
      if(isInsert)
        for(size_t k = 0; k < nPhysical; k++)
          physical.insert(pair<int, const MElement*>(myPhysical[k],
                                                     myElement[j]));
    }
  }
}

void GeoExtractor::vertexExtract(const MapElement& element, MapVertex& vertex){
  // Iterate on Elements
  MapElement::const_iterator end = element.end();
  MapElement::const_iterator  it = element.begin();

  for(; it != end; it++){
    // Get Current Element
    MElement* myElement = const_cast<MElement*>(it->first);

    // Iterate on Vertices
    const size_t N = myElement->getNumVertices();
    for(size_t i = 0; i < N; i++){
      // Take Current Vertex
      MVertex* myVertex = myElement->getVertex(i);

      // Insert
      vertex.insert(myVertex);
    }
  }
}

void GeoExtractor::edgeExtract(const MapElement& element, MapEntity& edge){
  // Iterate on Elements
  MapElement::const_iterator end = element.end();
  MapElement::const_iterator  it = element.begin();

  for(; it != end; it++){
    // Get Current Element
    MElement* myElement = const_cast<MElement*>(it->first);

    // Iterate on Edges
    const size_t N = myElement->getNumEdges();

    for(size_t i = 0; i < N; i++){
      // Take Current Edge
      const MEdge myEdge = myElement->getEdge(i);

      // Get Vertex Global ID
      vector<int> myVertex(2);
      myVertex[0] = myEdge.getVertex(0)->getNum();
      myVertex[1] = myEdge.getVertex(1)->getNum();

      // Try to Insert
      edge.insert(myVertex);
    }
  }
}

void GeoExtractor::faceExtract(const MapElement& element, MapEntity& face){
  // Iterate on Elements
  MapElement::const_iterator end = element.end();
  MapElement::const_iterator  it = element.begin();

  // Iterate on Elements
  for(; it != end; it++){
    // Get Current Element
    MElement* myElement = const_cast<MElement*>(it->first);

    // Iterate on Faces
    const size_t N = myElement->getNumFaces();

    for(size_t i = 0; i < N; i++){
      // Take Current Face
      const MFace myFace = myElement->getFace(i);

      // Get Vertex Global ID
      const int    nVertex = myFace.getNumVertices();
      vector<int> myVertex(nVertex);

      for(int i = 0; i < nVertex; i++)
        myVertex[i] = myFace.getVertex(i)->getNum();

      // Try to Insert
      face.insert(myVertex);
    }
  }
}
