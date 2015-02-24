#ifndef _MAPVERTEX_H_
#define _MAPVERTEX_H_

#include <cstring>
#include <map>
#include "MVertex.h"

/**
   @class MapVertex
   @brief Map linking vertices to integers.

   Map linking vertices to integers.

   @note
   MapVertex::iterator->first   is a const MVertex*
   MapVertex::iterator->second  is a positive integer
*/

class MapVertex{
 private:
  class VertexSort{
    public:
      bool operator()(const MVertex* a, const MVertex* b) const;
  };

 private:
  typedef std::map<const MVertex*, size_t, VertexSort> MpVertex;

 public:
  typedef MpVertex::iterator             iterator;
  typedef MpVertex::const_iterator const_iterator;

 private:
  MpVertex myMap;

 public:
   MapVertex(void);
  ~MapVertex(void);

  size_t size(void) const;
  void   clear(void);
  bool   insert(const MVertex* v);

  iterator       find(const MVertex* v);
  const_iterator find(const MVertex* v) const;

  iterator       begin(void);
  const_iterator begin(void) const;
  iterator       end(void);
  const_iterator end(void)   const;
};


/**
   @fn MapVertex::MapVertex
   Instanciates an empty MapVertex
   **

   @fn MapVertex::~MapVertex
   Deletes this MapVertex
   **

   @fn MapVertex::size
   @return Returns the size of this MapVertex
   **

   @fn MapVertex::clear
   Clears this MapVertex
   **

   @fn MapVertex::insert
   @param v A pointer to a MVertex

   Inserts the given MVertex in this MapVertex.
   Its mapped value is defaulted to 0.

   @return Returns true if a new MVertex was inserted
   and false if the MVertex existed already
   **

   @fn MapVertex::iterator MapVertex::find(const MVertex*)
   @param v A pointer to a MVertex

   Finds the given MVertex in this MapVertex

   @return If the given MVertex is found,
   this method returns an iterator to it.
   Otherwise an iterator pointing to MapVertex::end() is returned.
   **

   @fn MapVertex::const_iterator MapVertex::find(const MVertex*)
   @param v A pointer to a MVertex

   Finds the given MVertex in this MapVertex

   @return If the given MVertex is found,
   this method returns a const_iterator to it.
   Otherwise a const_iterator pointing to MapVertex::end() is returned.
   **

   @fn MapVertex::iterator MapVertex::begin(void)
   @return Returns an iterator to the beginning of this MapVertex
   **

   @fn MapVertex::const_iterator MapVertex::begin(void) const
   @return Returns a const_iterator to the beginning of this MapVertex
   **

   @fn MapVertex::iterator MapVertex::end(void)
   @return Returns an iterator to the end of this MapVertex
   **

   @fn MapVertex::const_iterator MapVertex::end(void) const
   @return Returns a const_iterator to the end of this MapVertex
*/


#endif
