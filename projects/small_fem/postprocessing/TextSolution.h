#ifndef _TEXTSOLUTION_H_
#define _TEXTSOLUTION_H_

#include <string>
#include <map>

/**
   @class TextSolution
   @brief Handles text post-processing field

   Handles text post-processing field

   A TextSolution can write a post-processing map in the
   <a href="http://www.geuz.org/gmsh">gmsh</a> .pos file format.

   The same instance of TextSolution may handle multiple text.
   Each text is associated to an integer called a 'step'.
   If two texts have the same step, the last one override to first one.
 */

class TextSolution{
 private:
  std::map<size_t, std::string> data;

 public:
   TextSolution(void);
  ~TextSolution(void);

  void clear(void);
  void addValues(size_t step, std::string text);

  void write(std::string fileName) const;
};


/**
   @fn TextSolution::TextSolution
   Instanciates a new TextSolution which is empty
   **

   @fn TextSolution::~TextSolution
   Deletes this TextSolution
   **

   @fn TextSolution::clear
   This TextSolution is now empty
   **

   @fn TextSolution::addValues
   @param step An integer value
   @param text A string

   Adds the given text to this TextSolution at the given step.
   **

   @fn TextSolution::write
   @param fileName A file name (without extension)

   Writes this TextSolution in
   <a href="http://www.geuz.org/gmsh">gmsh</a> .pos file format
   into the given file
 */

#endif
