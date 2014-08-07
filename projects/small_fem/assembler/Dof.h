#ifndef _DOF_H_
#define _DOF_H_

#include <string>

/**
   @class Dof
   @brief This class represent a degree of freedom

   This class represents a degree of freedom.

   A Dof is defined by a pair of two integers called (entity, type).
   By themselfs, these integers have no meaning, they just define a Dof.
*/


class Dof{
 private:
  int entity;
  int type;

 public:
   Dof(void);
   Dof(const Dof& other);
   Dof(int entity, int type);
  ~Dof(void);

  int getEntity(void) const;
  int getType(void)   const;

  void setEntity(int entity);
  void setType(int type);
  void setDof(int entity, int type);

  bool operator<(const Dof& other) const;
  bool operator>(const Dof& other) const;
  bool operator==(const Dof& other) const;

  std::string toString(void) const;
};


/**
   @fn Dof::Dof(void)
   Instanciates a new Dof with:
   @li entity = 0
   @li type = 0
   **

   @fn Dof::Dof(const Dof&)
   @param other An other Dof

   Instanciates a copy of the given Dof
   **

   @fn Dof::Dof(int, int)
   @param entity A natural number
   @param type A natural number

   Instanciates a new Dof with the given pair (entity, type)
   **

   @fn Dof::~Dof
   Deletes this Dof
   **

   @fn Dof::getEntity
   @return Returns the associated entity of this Dof
   **

   @fn Dof::getType
   @return Returns the associated type of this Dof
   **

   @fn Dof::setEntity
   @param entity A natural number

   Sets this Dof entity to the given value
   **

   @fn Dof::setType
   @param type A natural number

   Sets this Dof type to the given value
   **

   @fn Dof::setDof
   @param entity A natural number
   @param type A natural number

   Sets this Dof to the given pair (entity, type)
   **

   @fn bool Dof::operator<(const Dof& other) const
   @param other An other Dof to compare with the current one
   @return Returns :
   @li true if the current Dof is smaller than the other one
   @li false otherwise
   **

   @fn bool Dof::operator>(const Dof& other) const
   @param other An other Dof to compare with the current one
   @return Returns :
   @li true if the current Dof is greater than the other one
   @li false otherwise
   **

   @fn bool Dof::operator==(const Dof& other) const
   @param other An other Dof to compare with the current one
   @return Returns :
   @li true if the current Dof is equal to the other one
   @li false otherwise
   **

   @fn Dof::toString
   @return Returns the Dof string
   **
*/

//////////////////////
// Inline Functions //
//////////////////////

inline int Dof::getEntity(void) const{
  return entity;
}

inline int Dof::getType(void) const{
  return type;
}

inline void Dof::setEntity(int entity){
  this->entity = entity;
}

inline void Dof::setType(int type){
  this->type = type;
}

inline void Dof::setDof(int entity, int type){
  this->entity = entity;
  this->type   = type;
}

inline bool Dof::operator<(const Dof& other) const{
  return (entity < other.entity) ||
    ((entity == other.entity) && (type < other.type));
}

inline bool Dof::operator>(const Dof& other) const{
  return (entity > other.entity) ||
    ((entity == other.entity) && (type > other.type));
}

inline bool Dof::operator==(const Dof& other) const{
  return (entity == other.entity) && (type == other.type);
}

#endif
