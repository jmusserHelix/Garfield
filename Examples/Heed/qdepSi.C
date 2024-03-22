#include <iostream>
#include <fstream>
#include <cmath>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TAxis.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  MediumSilicon si;

  // Make a component
  ComponentConstant cmp;
  constexpr double width = 50.e-4;
  cmp.SetArea(-2, -2, 0, 2, 2, width);
  cmp.SetMedium(&si);
  cmp.SetElectricField(0., 0., 20.);

  // Make a sensor
  Sensor sensor;
  sensor.AddComponent(&cmp);

  TH1::StatOverflows(true);
  TH1F hNe("hNe", ";deposited charge [electrons];entries", 150, 0., 15000.);
  TH1F hNc("hNc", ";number of clusters;entries", 350, -0.5, 349.5);

  TrackHeed track(&sensor);
  track.SetParticle("pion");
  track.SetBetaGamma(10.);
  track.Initialise(&si, true);
  const unsigned int nTracks = 10000;
  for (unsigned int i = 0; i < nTracks; ++i) {
    if (i % 1000 == 0) std::cout << "Track " << i << "\n";
    track.NewTrack(0., 0., 0., 0., 0., 0., 1.);
    unsigned int nsum = 0;
    unsigned int ncls = 0;
    for (const auto& cluster : track.GetClusters()) {
      nsum += cluster.electrons.size();
      ++ncls;
    }
    hNe.Fill(nsum);
    hNc.Fill(ncls);
  }

  TCanvas cNe("cNe", "", 600, 600);
  hNe.Draw();
  cNe.Update();
  TCanvas cNc("cNc", "", 600, 600);
  hNc.Draw();
  cNc.Update();
  app.Run(true);

}
