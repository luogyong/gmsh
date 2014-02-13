#ifndef _GEOEXTRACTOR_H_
#define _GEOEXTRACTOR_H_

#include <vector>
#include <map>

#include "Comparators.h"
#include "GEntity.h"
#include "MElement.h"
#include "MVertex.h"
#include "MEdge.h"
#include "MFace.h"

/**
   @class GeoExtractor
   @brief Extraction of geometrical entities

   This class allows the extraction of geomtrical entities:
   @li Element (MElement) extraction from an entity (GEntity)
   @li Vertex (MVertex) extraction from a collection of elements (MElement)
   @li Edge (MEdge) extraction from a collection of elements (MElement)
   @li Face (MFace) extraction from a collection of elements (MElement)

   This class got only class methods, so it is not requiered to instanciate it.
*/

class GeoExtractor{
 public:
   GeoExtractor(void);
  ~GeoExtractor(void);

  static std::pair<
    std::map<const MElement*, size_t, ElementComparator>*,
    std::multimap<int, const MElement*>*
    >
    extractElement(const std::vector<GEntity*>& entity);

  static std::map<const MVertex*, size_t, VertexComparator>*
    extractVertex(const std::map<const MElement*,
                                 size_t,
                                 ElementComparator>& element);

  static std::map<const MEdge*, size_t, EdgeComparator>*
    extractEdge(const std::map<const MElement*,
                               size_t,
                               ElementComparator>& element);

  static std::map<const MFace*, size_t, FaceComparator>*
    extractFace(const std::map<const MElement*,
                               size_t,
                               ElementComparator>& element);

 private:
  static MEdge* copy(const MEdge& edge);
  static MFace* copy(const MFace& face);
};


/**
   @fn GeoExtractor::GeoExtractor
   Instantiates a new GeoExtractor
   (unneeded since GeoExtractor got only class methods)
   **

   @fn GeoExtractor::~GeoExtractor
   Deletes this GeoExtractor
   **

   @fn GeoExtractor::extractElement
   @param entity A vector of GEntity
   @return Returns an std::pair with:
   @li The first field containing a map with the MElement%s in
   the given entities (the mapped values are set to zero)
   @li The second field containing a multimap with the MElement%s
   of the first field, and the physicals of the MElement%s

   @see
   See <a href="http://www.geuz.org/gmsh">gmsh</a> documentation for physcials
   **

   @fn GeoExtractor::extractVertex
   @param element A map with MElement%s
   @return Returns a map with the MVertices in the given MElement%s
   (the mapped values are set to zero)
   **

   @fn GeoExtractor::extractEdge
   @param element A map with MElement%s
   @return Returns a map with the MEdge%s in the given MElement%s
   (the mapped values are set to zero)
   **

   @fn GeoExtractor::extractFace
   @param element A map with MElement%s
   @return Returns a map with the MFace%s in the given MElement%s
   (the mapped values are set to zero)
 */


#endif
