#ifndef _MAPENTITY_H_
#define _MAPENTITY_H_

#include <cstring>
#include <vector>
#include <map>

/**
   @class MapEntity
   @brief Map handling geometrical entities

   Map linking geometrical entities to integers,
   and indentifying these entities without taking orientation into account.

   For exemple, the two following edges, [0, 1] and [1, 0],
   will be treated as the same edge.

   An entity is defined by its global vertex numbers.

   @note
   MapEntity::iterator->first   is a std::vector<int>
   MapEntity::iterator->second  is a positive integer
*/

class MapEntity{
 private:
  class VectorSort{
    public:
      bool operator()(const std::vector<int>& a,
                      const std::vector<int>& b) const;
  };

 private:
  typedef std::map<std::vector<int>, size_t, VectorSort> MEntity;

 public:
  typedef MEntity::iterator             iterator;
  typedef MEntity::const_iterator const_iterator;

 private:
  MEntity myMap;

 public:
   MapEntity(void);
  ~MapEntity(void);

  size_t size(void) const;
  void   clear(void);
  bool   insert(const std::vector<int>& v);

  iterator       find(const std::vector<int>& v);
  const_iterator find(const std::vector<int>& v) const;

  iterator       begin(void);
  const_iterator begin(void) const;
  iterator       end(void);
  const_iterator end(void)   const;
};


/**
   @fn MapEntity::MapEntity
   Instanciates an empty MapEntity
   **

   @fn MapEntity::~MapEntity
   Deletes this MapEntity
   **

   @fn MapEntity::size
   @return Returns the size of this MapEntity
   **

   @fn MapEntity::clear
   Clears this MapEntity
   **

   @fn MapEntity::insert
   @param v A vector of vertex global ID

   Inserts the given vector in this MapEntity.
   Its mapped value is defaulted to 0.

   @return Returns true if a new vector was inserted
   and false if the vector existed already
   **

   @fn MapEntity::iterator MapEntity::find(const std::vector<int>&)
   @param v A vector of vertex global ID

   Finds the given vector in this MapEntity

   @return If the given vector is found,
   this method returns an iterator to it.
   Otherwise an iterator pointing to MapEntity::end() is returned.
   **

   @fn MapEntity::const_iterator MapEntity::find(const std::vector<int>&) const
   @param v A vector of vertex global ID

   Finds the given vector in this MapEntity

   @return If the given vector is found,
   this method returns a const_iterator to it.
   Otherwise a const_iterator pointing to MapEntity::end() is returned.
   **

   @fn MapEntity::iterator MapEntity::begin(void)
   @return Returns an iterator to the beginning of this MapEntity
   **

   @fn MapEntity::const_iterator MapEntity::begin(void) const
   @return Returns a const_iterator to the beginning of this MapEntity
   **

   @fn MapEntity::iterator MapEntity::end(void)
   @return Returns an iterator to the end of this MapEntity
   **

   @fn MapEntity::const_iterator MapEntity::end(void) const
   @return Returns a const_iterator to the end of this MapEntity
*/


#endif
