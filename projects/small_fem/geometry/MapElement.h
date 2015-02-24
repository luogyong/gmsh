#ifndef _MAPELEMENT_H_
#define _MAPELEMENT_H_

#include <cstring>
#include <map>
#include "MElement.h"

/**
   @class MapElement
   @brief Map linking elements to integers.

   Map linking elements to integers.

   @note
   MapElement::iterator->first   is a const MElement*
   MapElement::iterator->second  is a positive integer
*/

class MapElement{
 private:
  class ElementSort{
    public:
      bool operator()(const MElement* a, const MElement* b) const;
  };

 private:
  typedef std::map<const MElement*, size_t, ElementSort> MpElement;

 public:
  typedef MpElement::iterator             iterator;
  typedef MpElement::const_iterator const_iterator;

 private:
  MpElement myMap;

 public:
   MapElement(void);
  ~MapElement(void);

  size_t size(void) const;
  void   clear(void);
  bool   insert(const MElement* e);

  iterator       find(const MElement* e);
  const_iterator find(const MElement* e) const;

  iterator       begin(void);
  const_iterator begin(void) const;
  iterator       end(void);
  const_iterator end(void)   const;
};


/**
   @fn MapElement::MapElement
   Instanciates an empty MapElement
   **

   @fn MapElement::~MapElement
   Deletes this MapElement
   **

   @fn MapElement::size
   @return Returns the size of this MapElement
   **

   @fn MapElement::clear
   Clears this MapElement
   **

   @fn MapElement::insert
   @param e A pointer to a MElement

   Inserts the given MElement in this MapElement.
   Its mapped value is defaulted to 0.

   @return Returns true if a new MElement was inserted
   and false if the MElement existed already
   **

   @fn MapElement::iterator MapElement::find(const MElement*)
   @param e A pointer to a MElement

   Finds the given MElement in this MapElement

   @return If the given MElement is found,
   this method returns an iterator to it.
   Otherwise an iterator pointing to MapElement::end() is returned.
   **

   @fn MapElement::const_iterator MapElement::find(const MElement*)
   @param e A pointer to a MElement

   Finds the given MElement in this MapElement

   @return If the given MElement is found,
   this method returns a const_iterator to it.
   Otherwise a const_iterator pointing to MapElement::end() is returned.
   **

   @fn MapElement::iterator MapElement::begin(void)
   @return Returns an iterator to the beginning of this MapElement
   **

   @fn MapElement::const_iterator MapElement::begin(void) const
   @return Returns a const_iterator to the beginning of this MapElement
   **

   @fn MapElement::iterator MapElement::end(void)
   @return Returns an iterator to the end of this MapElement
   **

   @fn MapElement::const_iterator MapElement::end(void) const
   @return Returns a const_iterator to the end of this MapElement
*/


#endif
