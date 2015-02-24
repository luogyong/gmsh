#include "algorithm"
#include "MapEntity.h"

using namespace std;

MapEntity::MapEntity(void){
  clear();
}

MapEntity::~MapEntity(void){
}

void MapEntity::clear(void){
  myMap.clear();
}

size_t MapEntity::size(void) const{
  return myMap.size();
}

bool MapEntity::insert(const vector<int>& v){
  // Copy, sort, and insert //
  vector<int> copy(v);
  sort(copy.begin(), copy.end());

  pair<iterator, bool> isOK = myMap.insert(pair<vector<int>, size_t>(copy, 0));
  return isOK.second;
}

MapEntity::iterator MapEntity::find(const vector<int>& v){
  // Copy, sort, and find //
  vector<int> copy(v);
  sort(copy.begin(), copy.end());

  return myMap.find(copy);
}

MapEntity::const_iterator MapEntity::find(const vector<int>& v) const{
  // Copy, sort, and find //
  vector<int> copy(v);
  sort(copy.begin(), copy.end());

  return myMap.find(copy);
}

MapEntity::iterator MapEntity::begin(void){
  return myMap.begin();
}

MapEntity::const_iterator MapEntity::begin(void) const{
  return myMap.begin();
}

MapEntity::iterator MapEntity::end(void){
  return myMap.end();
}

MapEntity::const_iterator MapEntity::end(void) const{
  return myMap.end();
}

bool MapEntity::VectorSort::operator()(const vector<int>& a,
                                       const vector<int>& b) const{
  // Entity with less node are bigger //
  if(a.size() != b.size())
    return a.size() < b.size();

  // Compare, if same size //
  const int size = a.size();
  for(int i = 0; i < size; i++)
    if(a[i] != b[i])
      return a[i] < b[i];

  // If the two vector are the same, then 'a' is NOT smaller then 'b' //
  return false;
}
