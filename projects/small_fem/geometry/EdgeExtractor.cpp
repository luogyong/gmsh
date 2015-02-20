#include "GeoExtractor.h"

using namespace std;

set<const MEdge*, EdgeComparator>*
GeoExtractor::extractEdge(const map<const MElement*,
                                    size_t,
                                    ElementComparator>& element){
  // Init //
  set<const MEdge*, EdgeComparator>*
    edge = new set<const MEdge*, EdgeComparator>;

  // Get Edges //
  const map<const MElement*, size_t, ElementComparator>::const_iterator
    endE = element.end();

  map<const MElement*, size_t, ElementComparator>::const_iterator
    itE = element.begin();

  // Iterate on Elements
  for(; itE != endE; itE++){
    // Get Current Element
    MElement* myElement = const_cast<MElement*>(itE->first);

    // Iterate on Edges
    const size_t N = myElement->getNumEdges();

    for(size_t i = 0; i < N; i++){
      // Take Current Edge
      const MEdge myEdge = myElement->getEdge(i);

      // Make a copy (on heap)
      MEdge* edgeCopy = copy(myEdge);

      // Try to Insert
      pair<set<const MEdge*, EdgeComparator>::iterator, bool>
        insert = edge->insert(edgeCopy);

      // If Insertion is not a success,
      // Delete edgeCopy
      if(!insert.second)
        delete edgeCopy;
    }
  }

  // Return //
  return edge;
}
