#include <iostream>

#include "SmallFem.h"

#include "LineReferenceSpace.h"
#include "TriReferenceSpace.h"
#include "QuadReferenceSpace.h"
#include "TetReferenceSpace.h"
#include "HexReferenceSpace.h"
#include "PyrReferenceSpace.h"
#include "PriReferenceSpace.h"

using namespace std;

void compute(const Options& option){
  cout << "Line:" << flush;
  LineReferenceSpace line;
  cout << line.getNOrientation() << endl;

  cout << "Triangle:" << flush;
  TriReferenceSpace tri;
  cout << tri.getNOrientation() << endl;

  cout << "Quadrangle:" << flush;
  QuadReferenceSpace quad;
  cout << quad.getNOrientation() << endl;

  cout << "Tetrahedron:" << flush;
  TetReferenceSpace tet;
  cout << tet.getNOrientation() << endl;

  cout << "Hexahedron:" << flush;
  HexReferenceSpace hex;
  cout << hex.getNOrientation() << endl;

  cout << "Pyramid:" << flush;
  PyrReferenceSpace pyr;
  cout << pyr.getNOrientation() << endl;

  cout << "Prism:" << flush;
  PriReferenceSpace  pri;
  cout << pri.getNOrientation() << endl;
}

int main(int argc, char** argv){
  // SmallFEM //
  SmallFem::Keywords("");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
