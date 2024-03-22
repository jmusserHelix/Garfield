#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackDegrade.hh"
#include "Garfield/Plotting.hh"
#include "Garfield/Random.hh"

using namespace Garfield;

int main(int argc, char * argv[]) {

  randomEngine.Seed(123456);
  TApplication app("app", &argc, argv);
  SetDefaultStyle();

  TH1F hElectrons("hElectrons", "Number of electrons;number of electrons", 
                  200, 0, 200);
  TH1F hClusterSize("hClusterSize", "Cluster size;electrons / cluster", 
                    100, 0.5, 100.5);

  MediumMagboltz gas("ar", 90., "co2", 10.);

  constexpr double width = 1.;
  // Make a component
  ComponentConstant cmp;
  cmp.SetArea(0., -10., -10., width, 10., 10.);
  cmp.SetMedium(&gas);
  cmp.SetElectricField(100., 0., 0.);

  // Make a sensor
  Sensor sensor;
  sensor.AddComponent(&cmp);

  TrackDegrade track(&sensor);
  track.SetBetaGamma(3.);
  track.Initialise(&gas, true);
  for (unsigned int i = 0; i < 1000; ++i) {
    if (i % 10 == 0) std::cout << "Track " << i << "...\n";
    track.NewTrack(0., 0., 0., 0., 1., 0., 0.);
    unsigned int nsum = 0; 
    for (const auto& cluster : track.GetClusters()) {
      nsum += cluster.electrons.size();
      hClusterSize.Fill(cluster.electrons.size());
    }
    hElectrons.Fill(nsum); 
  } 

  TCanvas c1;
  hElectrons.Draw();
  c1.Update();

  TCanvas c2;
  hClusterSize.Draw();
  c2.SetLogy();
  c2.Update();

  app.Run(true); 

}
