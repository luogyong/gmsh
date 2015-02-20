#include "GeoExtractor.h"

using namespace std;

set<const MFace*, FaceComparator>*
GeoExtractor::extractFace(const map<const MElement*,
                                    size_t,
                                    ElementComparator>& element){
  // Init //
  set<const MFace*, FaceComparator>*
    face = new set<const MFace*, FaceComparator>;

  // Get Faces //
  const map<const MElement*, size_t, ElementComparator>::const_iterator
    endE = element.end();

  map<const MElement*, size_t, ElementComparator>::const_iterator
    itE = element.begin();

  // Iterate on Elements
  for(; itE != endE; itE++){
    // Get Current Element
    MElement* myElement = const_cast<MElement*>(itE->first);

    // Iterate on Faces
    const size_t N = myElement->getNumFaces();

    for(size_t i = 0; i < N; i++){
      // Take Current Face
      const MFace myFace = myElement->getFace(i);

      // Make a copy (on heap)
      MFace* faceCopy = copy(myFace);

      // Try to Insert
      pair<set<const MFace*, FaceComparator>::iterator, bool>
        insert = face->insert(faceCopy);

      // If Insertion is not a success,
      // Delete faceCopy
      if(!insert.second)
        delete faceCopy;
    }
  }

  // Return //
  return face;
}
