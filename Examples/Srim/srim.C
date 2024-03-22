#include <iostream>

#include <TApplication.h>
#include <TH1F.h>
#include <TCanvas.h>

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentConstant.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackSrim.hh"
#include "Garfield/Plotting.hh"

using namespace Garfield;

int main(int argc, char *argv[]) {

  // Application
  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();

  // Define the medium.
  MediumMagboltz gas;
  gas.SetComposition("ar");
  // Set temperature [K] and pressure [Torr].
  gas.SetPressure(760.0);
  gas.SetTemperature(293.15);

  // Make a component (with constant electric field).
  ComponentConstant cmp;
  // Define the active area and medium.
  constexpr double width = 10.;
  constexpr double length = 1.5;
  cmp.SetArea(0., -0.5 * width, -0.5 * width, length, 0.5 * width, 0.5 * width);
  cmp.SetMedium(&gas);
  cmp.SetElectricField(1000., 0., 0.);
  
  // Make a sensor.
  Sensor sensor;
  sensor.AddComponent(&cmp);
  
  // Create a track class and connect it to a sensor.
  TrackSrim tr(&sensor);
  // Read SRIM output from file.
  const std::string file = "Alpha_in_Ar.txt";
  if (!tr.ReadFile(file)) {
    std::cerr << "Reading SRIM file failed.\n";
    return 0;
  }
  // Set the initial kinetic energy of the particle (in eV).
  tr.SetKineticEnergy(1.47e6);
  // Specify how many electrons we want to be grouped to a cluster.
  tr.SetTargetClusterSize(500);
  // tr.SetClustersMaximum(1000);

  // Make some plots of the SRIM data.
  tr.PlotEnergyLoss();
  tr.PlotRange();
  tr.PlotStraggling();
  // Print a table of the SRIM data.
  tr.Print();

  // Setup histograms.
  TH1F* hX = new TH1F("hX", "x-end;x [cm];entries", 100, 0, 1.);
  TH1F* hY = new TH1F("hY", "y-end;y [cm];entries", 100, -0.5, 0.5);
  TH1F* hZ = new TH1F("hZ", "z-end;z [cm];entries", 100, -0.5, 0.5);
  TH1F* hNe = new TH1F("hNe", ";electrons;entries", 100, 54000, 56000);

  // Generate tracks.
  const unsigned int nTracks = 1000;
  for (unsigned int i = 0; i < nTracks; ++i) {
    if (!tr.NewTrack(0., 0., 0., 0., 1., 0., 0.)) {
      std::cerr << "Generating clusters failed; skipping this track.\n";
      continue;
    }
    // Retrieve the clusters.
    const auto& clusters = tr.GetClusters();
    if (clusters.empty()) continue;
    // Count the total number of electrons.
    unsigned int netot = 0;
    for (const auto& cluster : clusters) netot += cluster.n;
    const auto& last = clusters.back();
    hX->Fill(last.x);
    hY->Fill(last.y);
    hZ->Fill(last.z);
    hNe->Fill(netot);
  }

  // Plot the histograms.
  TCanvas* c = new TCanvas("c", "SRIM", 100, 100, 800, 800);
  c->Divide(2, 2);
  c->cd(1); hX->Draw();
  c->cd(2); hY->Draw();
  c->cd(3); hZ->Draw();
  c->cd(4); hNe->Draw();
  c->Update();
  std::cout << "Done.\n";

  // Start loop
  app.Run();
  return 0;
}
