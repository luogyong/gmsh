#include "MapElement.h"

using namespace std;

MapElement::MapElement(void){
  clear();
}

MapElement::~MapElement(void){
}

void MapElement::clear(void){
  myMap.clear();
}

size_t MapElement::size(void) const{
  return myMap.size();
}

bool MapElement::insert(const MElement* e){
  pair<iterator, bool> isOK = myMap.insert(pair<const MElement*, size_t>(e, 0));
  return isOK.second;
}

MapElement::iterator MapElement::find(const MElement* e){
  return myMap.find(e);
}

MapElement::const_iterator MapElement::find(const MElement* e) const{
  return myMap.find(e);
}

MapElement::iterator MapElement::begin(void){
  return myMap.begin();
}

MapElement::const_iterator MapElement::begin(void) const{
  return myMap.begin();
}

MapElement::iterator MapElement::end(void){
  return myMap.end();
}

MapElement::const_iterator MapElement::end(void) const{
  return myMap.end();
}

bool MapElement::ElementSort::operator()(const MElement* a,
                                         const MElement* b) const{
  return a->getNum() < b->getNum();
}
