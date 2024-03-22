#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ViewField.hh"
#include "Garfield/ViewCell.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  SetDefaultStyle();

  // Setup the gas.
  MediumMagboltz gas("ar");

  // Setup the cell layout.
  ComponentAnalyticField cmp;
  cmp.SetMedium(&gas);
  cmp.AddPlaneY(-1.,    0.);
  cmp.AddPlaneY(+1., 2000.);
  cmp.AddWire(0., 0., 1., 1000.);
  cmp.EnableDipoleTerms();
  cmp.PrintCell();

  // Plot the potential.
  TCanvas canvas("c", "", 600, 600);
  ViewCell cellView;
  cellView.SetCanvas(&canvas);
  cellView.SetComponent(&cmp);
  cellView.EnableWireMarkers(false);
  cellView.SetArea(-1.1, -1.1, 1.1, 1.1);
  ViewField fieldView;
  fieldView.SetCanvas(&canvas);
  fieldView.SetComponent(&cmp);
  fieldView.SetArea(-1.1, -1.1, 1.1, 1.1);
  fieldView.Plot("v", "cont1z");
  cellView.Plot2d();
  gSystem->ProcessEvents();
  app.Run(kTRUE);
}
