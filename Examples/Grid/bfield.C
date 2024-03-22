#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentGrid.hh"
#include "Garfield/Sensor.hh"

#include "Garfield/ViewField.hh"

using namespace Garfield;

int main(int argc, char *argv[]) {

  TApplication app("app", &argc, argv);

  // Load the field map.
  ComponentGrid cmp;
  cmp.SetCylindricalCoordinates();
  cmp.LoadMagneticField("solenoid.txt", "XZ");

  ViewField view;
  view.SetComponent(&cmp);
  view.SetMagneticFieldRange(5.1, 5.3);
  const double r = 5.;
  const double theta = 0.;
  const double x0 = r * cos(theta);
  const double y0 = r * sin(theta);
  view.PlotProfile(x0, y0, -20., x0, y0, 20., "bz"); 
  app.Run(true);
}
