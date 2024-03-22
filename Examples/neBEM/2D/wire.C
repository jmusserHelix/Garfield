#include <iostream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/ComponentNeBem2d.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);

  MediumMagboltz gas;
  gas.SetComposition("ar", 100.);
  
  ComponentNeBem2d cmp;
  constexpr double r = 2.;
  constexpr unsigned int n = 6;
  std::vector<double> xv(n, 0.);
  std::vector<double> yv(n, 0.);
  for (unsigned int i = 0; i < n; ++i) {
    const double phi = i * TwoPi / n; 
    xv[i] = r * cos(phi);
    yv[i] = r * sin(phi);
  }
  cmp.AddRegion(xv, yv, &gas, 1, 0.);
  cmp.AddWire(0, 0, 100.e-4, 5000.);

  cmp.SetNumberOfDivisions(100);
  cmp.Initialise();

  TCanvas canvas("c", "", 600, 600);
  ViewField fieldView;
  fieldView.SetCanvas(&canvas);
  fieldView.SetComponent(&cmp);
  fieldView.SetArea(-1.1 * r, -1.1 * r, 1.1 * r, 1.1 * r);
  fieldView.PlotContour();
  ViewCell cellView;
  cellView.SetCanvas(&canvas);
  cellView.SetComponent(&cmp);
  cellView.SetArea(-1.1 * r, -1.1 * r, 1.1 * r, 1.1 * r);
  cellView.Plot2d();
  app.Run(true);
}
