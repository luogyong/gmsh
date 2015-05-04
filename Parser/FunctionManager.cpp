// Gmsh - Copyright (C) 1997-2015 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@geuz.org>.

#include <map>
#include <stack>
#include <string>
#include "FunctionManager.h"

class File_Position
{
 public:
  int lineno;
  gmshfpos_t position;
  gmshFILE file;
  std::string filename;
};

class mystack
{
 public:
  std::stack<File_Position> s;
};

class mymap
{
 public:
  std::map<std::string, File_Position> m;
};

FunctionManager *FunctionManager::instance = 0;

FunctionManager::FunctionManager()
{
  functions = new mymap;
  calls = new mystack;
}

FunctionManager *FunctionManager::Instance()
{
  if(!instance) {
    instance = new FunctionManager;
  }
  return instance;
}

int FunctionManager::enterFunction(const std::string &name, gmshFILE * f,
                                   std::string &filename, int &lno) const
{
  if(functions->m.find(name) == functions->m.end())
    return 0;
  File_Position fpold;
  fpold.lineno = lno;
  fpold.filename = filename;
  fpold.file = *f;
  gmshgetpos(fpold.file, &fpold.position);
  calls->s.push(fpold);
  File_Position fp = (functions->m)[name];
  gmshsetpos(fp.file, &fp.position);
  *f = fp.file;
  filename = fp.filename;
  lno = fp.lineno;
  return 1;
}

int FunctionManager::leaveFunction(gmshFILE * f, std::string &filename, int &lno)
{
  if(!calls->s.size())
    return 0;
  File_Position fp;
  fp = calls->s.top();
  calls->s.pop();
  gmshsetpos(fp.file, &fp.position);
  *f = fp.file;
  filename = fp.filename;
  lno = fp.lineno;
  return 1;
}

int FunctionManager::createFunction(const std::string &name, gmshFILE f,
                                    const std::string &filename, int lno)
{
  if(functions->m.find(name) != functions->m.end())
    return 0;
  File_Position fp;
  fp.file = f;
  fp.filename = filename;
  fp.lineno = lno;
  gmshgetpos(fp.file, &fp.position);
  (functions->m)[name] = fp;
  return 1;
}
