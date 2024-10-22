#include <iostream>
#include <fstream>
#include <cstdlib>

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>

#include "Garfield/ViewDrift.hh"

#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/MediumMagboltz.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/DriftLineRKF.hh"
#include "Garfield/TrackHeed.hh"

using namespace Garfield;

bool readTransferFunction(Sensor& sensor) {

  std::ifstream infile;
  infile.open("mdt_elx_delta.txt", std::ios::in);
  if (!infile) {
    std::cerr << "Could not read delta response function.\n";
    return false;
  }
  std::vector<double> times;
  std::vector<double> values;
  while (!infile.eof()) {
    double t = 0., f = 0.;
    infile >> t >> f;
    if (infile.eof() || infile.fail()) break;
    times.push_back(1.e3 * t);
    values.push_back(f);
  }
  infile.close();
  sensor.SetTransferFunction(times, values);
  return true;
}

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
 
  // Make a gas medium.
  MediumMagboltz gas;
  gas.LoadGasFile("ar_93_co2_7_3bar.gas");
  auto installdir = std::getenv("GARFIELD_INSTALL");
  if (!installdir) {
    std::cerr << "GARFIELD_INSTALL variable not set.\n";
    return 1;
  }
  const std::string path = installdir;
  gas.LoadIonMobility(path + "/share/Garfield/Data/IonMobility_Ar+_Ar.txt");

  // Make a component with analytic electric field.
  ComponentAnalyticField cmp;
  cmp.SetMedium(&gas);
  // Wire radius [cm]
  const double rWire = 25.e-4;
  // Outer radius of the tube [cm]
  const double rTube = 0.71;
  // Voltages
  const double vWire = 2730.;
  const double vTube = 0.;
  // Add the wire in the centre.
  cmp.AddWire(0, 0, 2 * rWire, vWire, "s");
  // Add the tube.
  cmp.AddTube(rTube, vTube, 0);

  // Make a sensor.
  Sensor sensor;
  sensor.AddComponent(&cmp);
  sensor.AddElectrode(&cmp, "s");
  // Set the signal time window.
  const double tstep = 0.5;
  const double tmin = -0.5 * tstep;
  const unsigned int nbins = 1000;
  sensor.SetTimeWindow(tmin, tstep, nbins);
  // Set the delta reponse function.
  if (!readTransferFunction(sensor)) return 0;
  sensor.ClearSignal();

  // Set up Heed.
  TrackHeed track;
  track.SetParticle("muon");
  track.SetEnergy(170.e9);
  track.SetSensor(&sensor);
 
  TCanvas* cD = nullptr;
  ViewDrift driftView;
  constexpr bool plotDrift = true;
  if (plotDrift) {
    cD = new TCanvas("cD", "", 600, 600);
    driftView.SetCanvas(cD);
    track.EnablePlotting(&driftView);
  }
 
  TCanvas* cS = nullptr;
  constexpr bool plotSignal = true;
  if (plotSignal) cS = new TCanvas("cS", "", 600, 600);

  const double rTrack = 0.3;
  const double x0 = rTrack;
  const double y0 = -sqrt(rTube * rTube - rTrack * rTrack);
  const unsigned int nTracks = 1;
  for (unsigned int j = 0; j < nTracks; ++j) {
    sensor.ClearSignal();
    std::vector<std::array<double, 4> > electrons;
    track.NewTrack(x0, y0, 0, 0, 0, 1, 0);
    for (const auto& cluster : track.GetClusters()) {
      for (const auto& electron : cluster.electrons) {
        electrons.push_back({electron.x, electron.y, electron.z, 
                             electron.t});      
      }
    }
    // Loop over the primary electrons along the track.
    const std::size_t ne = electrons.size();
    #pragma omp parallel for
    for (size_t k = 0; k < ne; ++k) {
      DriftLineRKF drift(&sensor);
      drift.SetGainFluctuationsPolya(0., 20000., true);
      // drift.EnableIonTail();
      if (plotDrift) drift.EnablePlotting(&driftView);
      const double xe = electrons[k][0];
      const double ye = electrons[k][1];
      const double ze = electrons[k][2];
      const double te = electrons[k][3];
      drift.DriftElectron(xe, ye, ze, te);
    }

    if (plotDrift) {
      cD->Clear();
      cmp.PlotCell(cD);
      constexpr bool twod = true;
      constexpr bool drawaxis = false;
      driftView.Plot(twod, drawaxis);
    }
    sensor.ConvoluteSignals();
    int nt = 0;
    if (!sensor.ComputeThresholdCrossings(-2., "s", nt)) continue;
    if (plotSignal) sensor.PlotSignal("s", cS);
  }

  app.Run(kTRUE);

}
