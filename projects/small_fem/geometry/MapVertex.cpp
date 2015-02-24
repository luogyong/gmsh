#include "MapVertex.h"

using namespace std;

MapVertex::MapVertex(void){
  clear();
}

MapVertex::~MapVertex(void){
}

void MapVertex::clear(void){
  myMap.clear();
}

size_t MapVertex::size(void) const{
  return myMap.size();
}

bool MapVertex::insert(const MVertex* v){
  pair<iterator, bool> isOK = myMap.insert(pair<const MVertex*, size_t>(v, 0));
  return isOK.second;
}

MapVertex::iterator MapVertex::find(const MVertex* v){
  return myMap.find(v);
}

MapVertex::const_iterator MapVertex::find(const MVertex* v) const{
  return myMap.find(v);
}

MapVertex::iterator MapVertex::begin(void){
  return myMap.begin();
}

MapVertex::const_iterator MapVertex::begin(void) const{
  return myMap.begin();
}

MapVertex::iterator MapVertex::end(void){
  return myMap.end();
}

MapVertex::const_iterator MapVertex::end(void) const{
  return myMap.end();
}

bool MapVertex::VertexSort::operator()(const MVertex* a,const MVertex* b) const{
  return a->getNum() < b->getNum();
}
