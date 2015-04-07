#include <fstream>
#include <vector>
#include "TextSolution.h"

using namespace std;

TextSolution::TextSolution(void){
  data.clear();
}

TextSolution::~TextSolution(void){
  data.clear();
}

void TextSolution::clear(void){
  data.clear();
}

void TextSolution::addValues(size_t step, std::string text){
  data.insert(pair<size_t, string>(step, text));
}

void TextSolution::write(std::string fileName) const{
  // Serialize data //
  // Iterators
  map<size_t, string>::const_iterator  it = data.begin();
  map<size_t, string>::const_iterator end = data.end();

  // Full size
  size_t size;
  size_t sizeMinus;

  if(data.empty()){
    size      = 0;
    sizeMinus = 0;
  }

  else{
    end--;
    size      = end->first + 1;
    sizeMinus = size - 1;
  }

  // Allocate and populate
  vector<string> serial(size, "");

  it  = data.begin();
  end = data.end();
  for(; it != end; it++)
    serial[it->first] = it->second;

  // Write Stream //
  string fullName = fileName + ".pos";
  ofstream stream(fullName.c_str(), ofstream::out);

  stream << "View" "\"" << fileName << "\"" << " {" << endl
         << "  T2(100000,30,66574){";

  if(!data.empty()){
    for(size_t i = 0; i < sizeMinus; i++){
      stream << "\"" << serial[i] << "\"" << ",";
    }
    stream << "\"" << serial[sizeMinus] << "\"";
  }

  else{
    stream << "\"Empty\"";
  }

  stream << "};" << endl
         << "};" << endl;

  // Done !//
  stream.close();
}
