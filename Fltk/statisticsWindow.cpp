// Gmsh - Copyright (C) 1997-2015 C. Geuzaine, J.-F. Remacle
//
// See the LICENSE.txt file for license information. Please report all
// bugs and problems to the public mailing list <gmsh@geuz.org>.

#include <FL/Fl_Tabs.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Return_Button.H>
#include "FlGui.h"
#include "drawContext.h"
#include "statisticsWindow.h"
#include "paletteWindow.h"
#include "GModel.h"
#include "MElement.h"
#include "PView.h"
#include "Generator.h"
#include "Context.h"
#include "OS.h"
#include "Field.h"

enum QM_HISTO {QMH_SICN_XY, QMH_SICN_3D, QMH_GAMMA_XY, QMH_GAMMA_3D, QMH_RHO_XY, QMH_RHO_3D};

void statistics_cb(Fl_Widget *w, void *data)
{
  FlGui::instance()->stats->show();
}

static void statistics_update_cb(Fl_Widget *w, void *data)
{
  FlGui::instance()->stats->compute(true);
}

static void statistics_histogram_cb(Fl_Widget *w, void *data)
{
  QM_HISTO qmh = *(QM_HISTO*)data;

  std::vector<double> x, y;

  if (qmh == QMH_SICN_XY) {
    for(int i = 0; i < 100; i++){
      x.push_back((double)(2*i-99) / 99);
      y.push_back(FlGui::instance()->stats->quality[0][i]);
    }
    new PView("SICN", "# Elements", x, y);
  }
  else if (qmh == QMH_GAMMA_XY) {
    for(int i = 0; i < 100; i++){
      x.push_back((double)i / 99);
      y.push_back(FlGui::instance()->stats->quality[1][i]);
    }
    new PView("Gamma", "# Elements", x, y);
  }
  else if (qmh == QMH_RHO_XY) {
    for(int i = 0; i < 100; i++){
      x.push_back((double)i / 99);
      y.push_back(FlGui::instance()->stats->quality[2][i]);
    }
    new PView("Rho", "# Elements", x, y);
  }
  else {
    std::vector<GEntity*> entities_;
    GModel::current()->getEntities(entities_);
    std::map<int, std::vector<double> > d;
    for (unsigned int i = 0; i < entities_.size(); i++){
      if (entities_[i]->dim() < 2) continue;
      for (unsigned int j = 0; j < entities_[i]->getNumMeshElements(); j++) {
        MElement *e = entities_[i]->getMeshElement(j);
        if (qmh == QMH_SICN_3D)
          d[e->getNum()].push_back(e->minSICNShapeMeasure());
        else if (qmh == QMH_GAMMA_3D)
          d[e->getNum()].push_back(e->gammaShapeMeasure());
        else if (qmh == QMH_RHO_3D)
          d[e->getNum()].push_back(e->rhoShapeMeasure());
      }
    }
    std::string name = (qmh == QMH_SICN_3D) ? "SICN" :
                       (qmh == QMH_GAMMA_3D) ? "Gamma" :
                       (qmh == QMH_RHO_3D) ? "Rho" : "";
    new PView(name, "ElementData", GModel::current(), d);
  }

  FlGui::instance()->updateViews(true, true);
  drawContext::global()->draw();
}

statisticsWindow::statisticsWindow(int deltaFontSize)
{
  FL_NORMAL_SIZE -= deltaFontSize;

  int num = 0;
  int width = 26 * FL_NORMAL_SIZE;
  int height = 5 * WB + 17 * BH;

  win = new paletteWindow
    (width, height, CTX::instance()->nonModalWindows ? true : false, "Statistics");
  win->box(GMSH_WINDOW_BOX);
  {
    Fl_Tabs *o = new Fl_Tabs(WB, WB, width - 2 * WB, height - 3 * WB - BH);
    {
      group[0] = new Fl_Group
        (WB, WB + BH, width - 2 * WB, height - 3 * WB - 2 * BH, "Geometry");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 1 * BH, IW, BH, "Points");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 2 * BH, IW, BH, "Lines");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 3 * BH, IW, BH, "Surfaces");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 4 * BH, IW, BH, "Volumes");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 5 * BH, IW, BH, "Physical groups");
      group[0]->end();
    }
    {
      group[1] = new Fl_Group
        (WB, WB + BH, width - 2 * WB, height - 3 * WB - 2 * BH, "Mesh");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 1 * BH, IW, BH, "Nodes on Lines");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 2 * BH, IW, BH, "Nodes on surfaces");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 3 * BH, IW, BH, "Nodes in volumes");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 4 * BH, IW, BH, "Triangles");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 5 * BH, IW, BH, "Quadrangles");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 6 * BH, IW, BH, "Tetrahedra");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 7 * BH, IW, BH, "Hexahedra");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 8 * BH, IW, BH, "Prisms");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 9 * BH, IW, BH, "Pyramids");

      value[num++] = new Fl_Output(2 * WB, 2 * WB + 10 * BH, IW, BH, "Time for 1D mesh");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 11 * BH, IW, BH, "Time for 2D mesh");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 12 * BH, IW, BH, "Time for 3D mesh");

      value[num] = new Fl_Output(2 * WB, 2 * WB + 13 * BH, IW, BH, "SICN");
      value[num]->tooltip("~ signed inverse condition number"); num++;
      value[num] = new Fl_Output(2 * WB, 2 * WB + 14 * BH, IW, BH, "Gamma");
      value[num]->tooltip("~ inscribed_radius / circumscribed_radius (simplices)"); num++;
      value[num] = new Fl_Output(2 * WB, 2 * WB + 15 * BH, IW, BH, "Rho");
      value[num]->tooltip("~ min_edge_length / max_edge_length"); num++;

      for(int i = 0; i < 3; i++){
        int ww = 3 * FL_NORMAL_SIZE;
        new Fl_Box
          (FL_NO_BOX, width - 3 * ww - 2 * WB, 2 * WB + (13 + i) * BH, ww, BH, "Plot");
        butt[2 * i] = new Fl_Button
          (width - 2 * ww - 2 * WB, 2 * WB + (13 + i) * BH, ww, BH, "X-Y");
        butt[2 * i + 1] = new Fl_Button
          (width - ww - 2 * WB, 2 * WB + (13 + i) * BH, ww, BH, "3D");
      }
      static const QM_HISTO qmh0 = QMH_SICN_XY, qmh1 = QMH_SICN_3D, qmh2 = QMH_GAMMA_XY,
                            qmh3 = QMH_GAMMA_3D, qmh4 = QMH_RHO_XY, qmh5 = QMH_RHO_3D;
      butt[0]->callback(statistics_histogram_cb, (void*) &qmh0);
      butt[1]->callback(statistics_histogram_cb, (void*) &qmh1);
      butt[2]->callback(statistics_histogram_cb, (void*) &qmh2);
      butt[3]->callback(statistics_histogram_cb, (void*) &qmh3);
      butt[4]->callback(statistics_histogram_cb, (void*) &qmh4);
      butt[5]->callback(statistics_histogram_cb, (void*) &qmh5);

      group[1]->end();
    }
    {
      group[2] = new Fl_Group
        (WB, WB + BH, width - 2 * WB, height - 3 * WB - 2 * BH, "Post-processing");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 1 * BH, IW, BH, "Views");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 2 * BH, IW, BH, "Points");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 3 * BH, IW, BH, "Lines");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 4 * BH, IW, BH, "Triangles");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 5 * BH, IW, BH, "Quadrangles");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 6 * BH, IW, BH, "Tetrahedra");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 7 * BH, IW, BH, "Hexahedra");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 8 * BH, IW, BH, "Prisms");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 9 * BH, IW, BH, "Pyramids");
      value[num++] = new Fl_Output(2 * WB, 2 * WB + 10 * BH, IW, BH, "Strings");
      group[2]->end();
    }
    o->end();
  }

  for(int i = 0; i < num; i++) {
    value[i]->align(FL_ALIGN_RIGHT);
    value[i]->value(0);
  }

  {
    memUsage = new Fl_Box(WB, height - BH - WB, width / 2, BH, "");
    memUsage->align(FL_ALIGN_INSIDE);

    Fl_Return_Button *o = new Fl_Return_Button
      (width - BB - WB, height - BH - WB, BB, BH, "Update");
    o->callback(statistics_update_cb);
  }

  win->position(CTX::instance()->statPosition[0], CTX::instance()->statPosition[1]);
  win->end();

  FL_NORMAL_SIZE += deltaFontSize;
}

void statisticsWindow::compute(bool elementQuality)
{
  int num = 0;
  static double s[50];
  static char label[50][256];

  if(elementQuality)
    GetStatistics(s, quality);
  else
    GetStatistics(s);

  // geom
  sprintf(label[num], "%g", s[0]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[1]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[2]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[3]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[45]); value[num]->value(label[num]); num++;

  // mesh
  sprintf(label[num], "%g", s[4]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[5]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[6]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[7]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[8]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[9]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[10]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[11]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[12]); value[num]->value(label[num]); num++;

  sprintf(label[num], "%g", s[13]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[14]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[15]); value[num]->value(label[num]); num++;

  if(!elementQuality){
    for(int i = 0; i < 6; i += 2) butt[i]->deactivate();
    sprintf(label[num], "Press Update");
    value[num]->deactivate();
    value[num]->value(label[num]); num++;
    sprintf(label[num], "Press Update");
    value[num]->deactivate();
    value[num]->value(label[num]); num++;
    sprintf(label[num], "Press Update");
    value[num]->deactivate();
    value[num]->value(label[num]); num++;
  }
  else{
    for(int i = 0; i < 6; i += 2) butt[i]->activate();
    sprintf(label[num], "%.4g (%.4g->%.4g)", s[17], s[18], s[19]);
    value[num]->activate();
    value[num]->value(label[num]); num++;
    sprintf(label[num], "%.4g (%.4g->%.4g)", s[20], s[21], s[22]);
    value[num]->activate();
    value[num]->value(label[num]); num++;
    sprintf(label[num], "%.4g (%.4g->%.4g)", s[23], s[24], s[25]);
    value[num]->activate();
    value[num]->value(label[num]); num++;
  }

  // post
  sprintf(label[num], "%g", s[26]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[27]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[28]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[29]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[30]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[31]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[32]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[33]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[34]); value[num]->value(label[num]); num++;
  sprintf(label[num], "%g", s[35]); value[num]->value(label[num]); num++;

  static char mem[256];
  long m = GetMemoryUsage();
  if(m){
    sprintf(mem, "Memory usage: %gMb", GetMemoryUsage()/1024./1024.);
    memUsage->label(mem);
  }
}

void statisticsWindow::show()
{
  if(!win->shown())
    compute(false);

  for(int i = 0; i < 3; i++)
    group[i]->hide();

  if(GModel::current()->getMeshStatus(true) > 0)
    group[1]->show();
  else if(PView::list.size())
    group[2]->show();
  else
    group[0]->show();

  win->show();
}
