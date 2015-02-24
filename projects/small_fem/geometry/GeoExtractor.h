#ifndef _GEOEXTRACTOR_H_
#define _GEOEXTRACTOR_H_

#include "MapElement.h"
#include "MapVertex.h"
#include "MapEntity.h"
#include "GEntity.h"

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

  static void elementExtract(const std::vector<GEntity*>& entity,
                             MapElement& element,
                             std::multimap<int, const MElement*>& physical);

  static void  vertexExtract(const MapElement& element, MapVertex& vertex);
  static void    edgeExtract(const MapElement& element, MapEntity& edge);
  static void    faceExtract(const MapElement& element, MapEntity& face);
};


/**
   @fn GeoExtractor::GeoExtractor
   Instantiates a new GeoExtractor
   (unneeded since GeoExtractor got only class methods)
   **

   @fn GeoExtractor::~GeoExtractor
   Deletes this GeoExtractor
   **

   @fn GeoExtractor::elementExtract
   @param entity A vector of GEntity
   @param element A MapElement
   @param physical A map linking a MElement to its physical

   Populates the MapElement with the MElement%s in the given GEntity%s
   (the mappend value is defaulted to 0).

   Populates the physical map with the MElement%s and
   the physical tags of the given GEntity%s.

   @see
   See <a href="http://www.geuz.org/gmsh">gmsh</a> documentation for physcials
   **

   @fn GeoExtractor::vertexExtract
   @param element A MapElement
   @param vertex A MapVertex

   Populate the MapVertex with the MVertex of the MElement
   in the given MapElement (the mapped values are set to zero)
   **

   @fn GeoExtractor::edgeExtract
   @param element A MapElement
   @param edge A MapEntity

   Populate the MapEntity with the MEdge of the MElement
   in the given MapElement (the mapped values are set to zero)
   **

   @fn GeoExtractor::faceExtract
   @param element A MapElement
   @param face A MapEntity

   Populate the MapEntity with the MFace of the MElement
   in the given MapElement (the mapped values are set to zero)
 */


#endif
