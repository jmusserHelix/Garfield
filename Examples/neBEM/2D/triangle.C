#include <iostream>
#include <iomanip>
#include <fstream>
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
  cmp.SetMedium(&gas);
  constexpr unsigned int nDiv = 200;
  cmp.SetNumberOfDivisions(nDiv);
  cmp.AddSegment(-1., 0., 0., 1., 1.);
  cmp.AddSegment( 1., 0., 0., 1., 1.);
  cmp.AddSegment(-1., 0., 1., 0., 0.);
  cmp.Initialise();

  std::ofstream outfile;
  outfile.open("triangle.txt", std::ios::out);
  outfile << "# " << nDiv << " divisions\n";
  for (int i = 1; i < 10; ++i) {
    const double y = i * 0.1;
    const double v = cmp.ElectricPotential(0., y, 0.);
    outfile << y << "  " << std::setprecision(9) << v << "\n";
  }
  outfile.close();

  TCanvas canvas("c", "", 600, 600);
  ViewField fieldView;
  fieldView.SetCanvas(&canvas);
  fieldView.SetComponent(&cmp);
  fieldView.SetArea(-1.1, -0.6, 1.1, 1.6);
  fieldView.PlotContour();
  ViewCell cellView;
  cellView.SetCanvas(&canvas);
  cellView.SetComponent(&cmp);
  cellView.SetArea(-1.1, -0.6, 1.1, 1.6);
  cellView.Plot2d();
  app.Run(true);
}
