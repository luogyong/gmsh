#include <iostream>

#include "LineReferenceSpace.h"
#include "TriReferenceSpace.h"
#include "QuadReferenceSpace.h"
#include "TetReferenceSpace.h"
#include "HexReferenceSpace.h"
#include "PyrReferenceSpace.h"
#include "PriReferenceSpace.h"

#include "Timer.h"
#include "SmallFem.h"

using namespace std;

void compute(const Options& option){
  Timer time;

  cout << "Line: " << flush;
  time.start();
  LineReferenceSpace line;
  time.stop();
  cout << line.getNOrientation()
       << " (" << time.time() << " " << time.unit() << ")" << endl << flush;

  cout << "Triangle: " << flush;
  time.start();
  TriReferenceSpace tri;
  time.stop();
  cout << tri.getNOrientation()
       << " (" << time.time() << " " << time.unit() << ")" << endl << flush;

  cout << "Quadrangle: " << flush;
  time.start();
  QuadReferenceSpace quad;
  time.stop();
  cout << quad.getNOrientation()
       << " (" << time.time() << " " << time.unit() << ")" << endl << flush;

  cout << "Tetrahedron: " << flush;
  time.start();
  TetReferenceSpace tet;
  time.stop();
  cout << tet.getNOrientation()
       << " (" << time.time() << " " << time.unit() << ")" << endl << flush;

  cout << "Hexahedron: " << flush;
  time.start();
  HexReferenceSpace hex;
  time.stop();
  cout << hex.getNOrientation()
       << " (" << time.time() << " " << time.unit() << ")" << endl << flush;

  cout << "Pyramid: " << flush;
  time.start();
  PyrReferenceSpace pyr;
  time.stop();
  cout << pyr.getNOrientation()
       << " (" << time.time() << " " << time.unit() << ")" << endl << flush;

  cout << "Prism: " << flush;
  time.start();
  PriReferenceSpace  pri;
  time.stop();
  cout << pri.getNOrientation()
       << " (" << time.time() << " " << time.unit() << ")" << endl << flush;
}

int main(int argc, char** argv){
  // SmallFEM //
  SmallFem::Keywords("");
  SmallFem::Initialize(argc, argv);

  compute(SmallFem::getOptions());

  SmallFem::Finalize();
  return 0;
}
