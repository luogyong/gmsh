#ifndef _GROUPOFELEMENT_H_
#define _GROUPOFELEMENT_H_

#include <string>
#include <vector>
#include <list>

#include "Mesh.h"
#include "MElement.h"
#include "ReferenceSpaceManager.h"

/**
   @class GroupOfElement
   @brief A Group of MElement%s

   This class is collection of discrete elements (MElement%s).

   The elements are sorted by geomtrical types and by orientations
*/

class Mesh;

class GroupOfElement{
 private:
  static const size_t nGeoType;

 private:
  const Mesh* mesh;

  std::vector<const MElement*>      element;
  std::vector<std::vector<size_t> > orientationStat;
  std::vector<size_t>               typeStat;

 public:
   GroupOfElement(std::list<const MElement*>& element, const Mesh& mesh);
   GroupOfElement(const Mesh& mesh);
  ~GroupOfElement(void);

  void add(const GroupOfElement& other);

  bool            isEmpty(void)   const;
  size_t          getNumber(void) const;
  const MElement& get(size_t i)   const;

  const std::vector<const MElement*>& getAll(void) const;
  const Mesh& getMesh(void) const;

  void getAllVertex(std::set<const MVertex*, VertexComparator>& vertex) const;
  void getAllVertexCoordinate(fullMatrix<double>& coord) const;

  const std::vector<size_t>& getOrientationStats(size_t elementType) const;
  const std::vector<size_t>& getTypeStats(void) const;

  std::pair<bool, size_t> isUniform(void) const;

  std::string toString(void) const;

 private:
  void init(void);
  static bool sortPredicate(const MElement* a, const MElement* b);
};


/**
   @fn GroupOfElement::GroupOfElement(std::list<const MElement*>&, const Mesh&)
   @param element A list of element
   @param mesh A Mesh

   Instantiates a new GroupOfElement,
   with the given elements and associated to the given Mesh
   **

   @fn GroupOfElement::GroupOfElement(const Mesh&)
   @param mesh A Mesh

   Instantiates a new empty GroupOfElement associated to the given Mesh
   **

   @fn GroupOfElement::~GroupOfElement
   Deletes this GroupOfElement
   **

   @fn GroupOfElement::add
   @param other An other GroupOfElement

   Adds the elements of the given GroupOfElement in this GroupOfElement
   **

   @fn GroupOfElement::isEmpty
   @return Returns true if this GroupOfElement has no elements
   and false otherwise
   **

   @fn GroupOfElement::getNumber
   @return Returns the number of elements in this GroupOfElement
   **

   @fn GroupOfElement::get
   @param i An interger ranging from 0 to GroupOfElement::getNumber() - 1
   @return Returns the ith element of the GroupOfElement
   **

   @fn GroupOfElement::getAll
   @return Returns all the elements of this GroupOfElement

   The elements are sorted by geomtrical types and by orientations
   **

   @fn GroupOfElement::getMesh
   @return Returns the associated Mesh
   **

   @fn GroupOfElement::getAllVertex
   @param vertex A set of MVertex

   Populates the given set with the MVertex%s of this GroupOfElement
   **

   @fn GroupOfElement::getAllVertexCoordinate
   @param coord A matrix

   Populates the given matrix with the coordinates
   of the MVertex%s of this GroupOfElement:
   @li The ith row is the ith MVertex
   @li The jth column is the jth dimension
   **

   @fn GroupOfElement::getOrientationStats
   @param elementType A geomtrical type
   @return A vector where the i-th entry is the number
   of element in GroupOfElement::getAll() (of the given geomtrical type)
   with a ReferenceSpaceManager::getOrientation() equal to i

   @see
   See <a href="http://www.geuz.org/gmsh">gmsh</a>
   documentation for geomtrical types
   **

   @fn GroupOfElement::getTypeStats
   @return A vector where the i-th entry is the number
   of element in GroupOfElement::getAll() with a geomtrical type equal to i

   @see
   See <a href="http://www.geuz.org/gmsh">gmsh</a>
   documentation for geomtrical types
   **

   @fn GroupOfElement::isUniform
   @return Returns a pair such that:
   @li If this GroupOfElement is composed only
   of element with the same geomtrical type,
   the first entry is true and the second is set to the common geomtrical type
   @li Otherwise the first entry is set to false and the second one is undefined
   **

   @fn GroupOfElement::toString
   @return Returns a string discribing this Group
*/


//////////////////////
// Inline Functions //
//////////////////////

inline bool GroupOfElement::sortPredicate(const MElement* a, const MElement* b){
  return
    ((a->getType()  < b->getType())) ||
    ((a->getType() == b->getType()) &&
     (ReferenceSpaceManager::getOrientation(*a) <
      ReferenceSpaceManager::getOrientation(*b)));
}

inline bool GroupOfElement::isEmpty(void) const{
  return getNumber() == 0;
}

inline size_t GroupOfElement::getNumber(void) const{
  return element.size();
}

inline const MElement& GroupOfElement::get(size_t i) const{
  return *element[i];
}

inline const std::vector<const MElement*>& GroupOfElement::getAll(void) const{
  return element;
}

inline const Mesh& GroupOfElement::getMesh(void) const{
  return *mesh;
}

inline const std::vector<size_t>&
GroupOfElement::getOrientationStats(size_t elementType) const{
  return orientationStat[elementType];
}

inline const std::vector<size_t>& GroupOfElement::getTypeStats(void) const{
  return typeStat;
}

#endif
