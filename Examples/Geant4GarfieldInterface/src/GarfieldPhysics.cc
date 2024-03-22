//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file GarfieldPhysics.cc
/// \brief Implementation of the GarfieldPhysics class
#include "GarfieldPhysics.hh"
#include "GarfieldAnalysis.hh"

#include "Garfield/AvalancheMC.hh"
#include "Garfield/AvalancheMicroscopic.hh"

GarfieldPhysics* GarfieldPhysics::fGarfieldPhysics = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

GarfieldPhysics* GarfieldPhysics::GetInstance() {
  if (!fGarfieldPhysics) fGarfieldPhysics = new GarfieldPhysics();
  return fGarfieldPhysics;
}

void GarfieldPhysics::Dispose() {
  delete fGarfieldPhysics;
  fGarfieldPhysics = nullptr;
}

GarfieldPhysics::~GarfieldPhysics() {
  delete fMediumMagboltz;
  delete fSensor;
  delete fComponentAnalyticField;
  delete fTrackHeed;

  std::cout << "Deconstructor GarfieldPhysics" << std::endl;
}

std::string GarfieldPhysics::GetIonizationModel() { return fIonizationModel; }

void GarfieldPhysics::SetIonizationModel(std::string model, bool useDefaults) {
  if (model != "PAIPhot" && model != "PAI" && model != "Heed") {
    std::cout << "Unknown ionization model " << model << std::endl;
    std::cout << "Using PAIPhot as default model!" << std::endl;
    model = "PAIPhot";
  }
  fIonizationModel = model;

  if (fIonizationModel == "PAIPhot" || fIonizationModel == "PAI") {
    if (useDefaults) {
      // Particle types and energies for which the G4FastSimulationModel with
      // Garfield++ is valid
      this->AddParticleName("e-", 1e-6, 1e-3, "garfield");
      this->AddParticleName("gamma", 1e-6, 1e+8, "garfield");

      // Particle types and energies for which the PAI or PAIPhot model is valid
      this->AddParticleName("e-", 0, 1e+8, "geant4");
      this->AddParticleName("e+", 0, 1e+8, "geant4");
      this->AddParticleName("mu-", 0, 1e+8, "geant4");
      this->AddParticleName("mu+", 0, 1e+8, "geant4");
      this->AddParticleName("proton", 0, 1e+8, "geant4");
      this->AddParticleName("pi+", 0, 1e+8, "geant4");
      this->AddParticleName("pi-", 0, 1e+8, "geant4");
      this->AddParticleName("alpha", 0, 1e+8, "geant4");
      this->AddParticleName("He3", 0, 1e+8, "geant4");
      this->AddParticleName("GenericIon", 0, 1e+8, "geant4");
    }

  } else if (fIonizationModel == "Heed") {
    if (useDefaults) {
      // Particle types and energies for which the G4FastSimulationModel with
      // Garfield++ is valid
      this->AddParticleName("gamma", 1e-6, 1e+8, "garfield");
      this->AddParticleName("e-", 6e-2, 1e+7, "garfield");
      this->AddParticleName("e+", 6e-2, 1e+7, "garfield");
      this->AddParticleName("mu-", 1e+1, 1e+8, "garfield");
      this->AddParticleName("mu+", 1e+1, 1e+8, "garfield");
      this->AddParticleName("pi-", 2e+1, 1e+8, "garfield");
      this->AddParticleName("pi+", 2e+1, 1e+8, "garfield");
      this->AddParticleName("kaon-", 1e+1, 1e+8, "garfield");
      this->AddParticleName("kaon+", 1e+1, 1e+8, "garfield");
      this->AddParticleName("proton", 9.e+1, 1e+8, "garfield");
      this->AddParticleName("anti_proton", 9.e+1, 1e+8, "garfield");
      this->AddParticleName("deuteron", 2.e+2, 1e+8, "garfield");
      this->AddParticleName("alpha", 4.e+2, 1e+8, "garfield");
    }
  }
}

void GarfieldPhysics::AddParticleName(const std::string particleName,
                                      double ekin_min_MeV, double ekin_max_MeV,
                                      std::string program) {
  if (ekin_min_MeV >= ekin_max_MeV) {
    std::cout << "Ekin_min=" << ekin_min_MeV
              << " keV is larger than Ekin_max=" << ekin_max_MeV << " keV"
              << std::endl;
    return;
  }

  if (program == "garfield") {
    std::cout << "Garfield model (Heed) is applicable for G4Particle "
              << particleName << " between " << ekin_min_MeV << " MeV and "
              << ekin_max_MeV << " MeV" << std::endl;

    fMapParticlesEnergyGarfield.insert(std::make_pair(
        particleName, std::make_pair(ekin_min_MeV, ekin_max_MeV)));
  } else {
    std::cout << fIonizationModel << " is applicable for G4Particle "
              << particleName << " between " << ekin_min_MeV << " MeV and "
              << ekin_max_MeV << " MeV" << std::endl;
    fMapParticlesEnergyGeant4.insert(std::make_pair(
        particleName, std::make_pair(ekin_min_MeV, ekin_max_MeV)));
  }
}

bool GarfieldPhysics::FindParticleName(std::string name, std::string program) {
  if (program == "garfield") {
    auto it = fMapParticlesEnergyGarfield.find(name);
    if (it != fMapParticlesEnergyGarfield.end()) return true;
  } else {
    auto it = fMapParticlesEnergyGeant4.find(name);
    if (it != fMapParticlesEnergyGeant4.end()) return true;
  }
  return false;
}

bool GarfieldPhysics::FindParticleNameEnergy(std::string name, double ekin_MeV,
                                             std::string program) {
  if (program == "garfield") {
    auto it = fMapParticlesEnergyGarfield.find(name);
    if (it != fMapParticlesEnergyGarfield.end()) {
      EnergyRange_MeV range = it->second;
      if (range.first <= ekin_MeV && range.second >= ekin_MeV) {
        return true;
      }
    }
  } else {
    auto it = fMapParticlesEnergyGeant4.find(name);
    if (it != fMapParticlesEnergyGeant4.end()) {
      EnergyRange_MeV range = it->second;
      if (range.first <= ekin_MeV && range.second >= ekin_MeV) {
        return true;
      }
    }
  }
  return false;
}

double GarfieldPhysics::GetMinEnergyMeVParticle(std::string name,
                                                std::string program) {
  if (program == "garfield") {
    auto it = fMapParticlesEnergyGarfield.find(name);
    if (it != fMapParticlesEnergyGarfield.end()) {
      EnergyRange_MeV range = it->second;
      return range.first;
    }
  } else {
    auto it = fMapParticlesEnergyGeant4.find(name);
    if (it != fMapParticlesEnergyGeant4.end()) {
      EnergyRange_MeV range = it->second;
      return range.first;
    }
  }
  return -1;
}

double GarfieldPhysics::GetMaxEnergyMeVParticle(std::string name,
                                                std::string program) {
  if (program == "garfield") {
    auto it = fMapParticlesEnergyGarfield.find(name);
    if (it != fMapParticlesEnergyGarfield.end()) {
      EnergyRange_MeV range = it->second;
      return range.second;
    }
  } else {
    auto it = fMapParticlesEnergyGeant4.find(name);
    if (it != fMapParticlesEnergyGeant4.end()) {
      EnergyRange_MeV range = it->second;
      return range.second;
    }
  }
  return -1;
}

void GarfieldPhysics::InitializePhysics() {
  // Define the gas mixture.
  fMediumMagboltz = new Garfield::MediumMagboltz();
  fMediumMagboltz->SetComposition("ar", 70., "co2", 30.);
  fMediumMagboltz->SetTemperature(293.15);
  fMediumMagboltz->SetPressure(760.);
  fMediumMagboltz->Initialise(true);
  // Set the Penning transfer efficiency.
  const double rPenning = 0.57;
  const double lambdaPenning = 0.;
  fMediumMagboltz->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
  fMediumMagboltz->LoadGasFile("ar_70_co2_30_1000mbar.gas");

  fComponentAnalyticField = new Garfield::ComponentAnalyticField();
  fComponentAnalyticField->SetMedium(fMediumMagboltz);
  // Wire radius [cm]
  constexpr double rWire = 25.e-4;
  // Tube radius [cm]
  constexpr double rTube = 1.451;
  // Voltages
  constexpr double vWire = 1000.;
  constexpr double vTube = 0.;
  // Add the wire in the center.
  fComponentAnalyticField->AddWire(0., 0., 2 * rWire, vWire, "w");
  // Add the tube.
  fComponentAnalyticField->AddTube(rTube, vTube, 0, "t");

  fSensor = new Garfield::Sensor();
  fSensor->AddComponent(fComponentAnalyticField);

  fTrackHeed = new Garfield::TrackHeed();
  fTrackHeed->SetSensor(fSensor);
  fTrackHeed->EnableDeltaElectronTransport();
}

void GarfieldPhysics::DoIt(std::string particleName, double ekin_MeV,
                           double time, double x_cm, double y_cm, double z_cm,
                           double dx, double dy, double dz) {
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  fEnergyDeposit = 0;
  fSecondaryParticles.clear();

  Garfield::AvalancheMC drift;
  drift.SetSensor(fSensor);
  drift.SetDistanceSteps(1.e-4);
  Garfield::AvalancheMicroscopic avalanche;
  avalanche.SetSensor(fSensor);

  // Wire radius [cm]
  constexpr double rWire = 25.e-4;
  // Outer radius of the tube [cm]
  constexpr double rTube = 1.45;
  // Half-length of the tube [cm]
  constexpr double lTube = 10.;

  double eKin_eV = ekin_MeV * 1e+6;

  fEnergyDeposit = 0;
  if (fIonizationModel != "Heed" || particleName == "gamma") {
    // Number of electrons produced
    int nc = 0;
    if (particleName == "gamma") {
      fTrackHeed->TransportPhoton(x_cm, y_cm, z_cm, time, eKin_eV, dx, dy, dz,
                                  nc);
    } else {
      fTrackHeed->TransportDeltaElectron(x_cm, y_cm, z_cm, time, eKin_eV, dx,
                                         dy, dz, nc);
      fEnergyDeposit = eKin_eV;
    }

    for (int cl = 0; cl < nc; cl++) {
      double xe, ye, ze, te;
      double ee, dxe, dye, dze;
      fTrackHeed->GetElectron(cl, xe, ye, ze, te, ee, dxe, dye, dze);
      if (fabs(ze) > lTube || sqrt(xe * xe + ye * ye) > rTube) continue;
      nsum++;
      if (particleName == "gamma") {
        fEnergyDeposit += fTrackHeed->GetW();
      }
      analysisManager->FillH3(1, ze * 10, xe * 10, ye * 10);
      if (createSecondariesInGeant4) {
        double newTime = te;
        if (newTime < time) {
          newTime += time;
        }
        fSecondaryParticles.emplace_back(GarfieldParticle(
            "e-", ee, newTime, xe, ye, ze, dxe, dye, dze));
      }

      drift.DriftElectron(xe, ye, ze, te);

      double xe1, ye1, ze1, te1;
      double xe2, ye2, ze2, te2;

      int status;
      drift.GetElectronEndpoint(0, xe1, ye1, ze1, te1, xe2, ye2, ze2, te2,
                                status);

      if (0 < xe2 && xe2 < rWire) {
        xe2 += 2 * rWire;
      } else if (0 > xe2 && xe2 > -rWire) {
        xe2 += -2 * rWire;
      } 
      if (0 < ye2 && ye2 < rWire) {
        ye2 += 2 * rWire;
      } else if (0 > ye2 && ye2 > -rWire) {
        ye2 += -2 * rWire;
      }

      double e2 = 0.1;
      avalanche.AvalancheElectron(xe2, ye2, ze2, te2, e2, 0, 0, 0);

      int ne = 0, ni = 0;
      avalanche.GetAvalancheSize(ne, ni);
      fAvalancheSize += ne;
    }
  } else {
    fTrackHeed->SetParticle(particleName);
    fTrackHeed->SetKineticEnergy(eKin_eV);
    fTrackHeed->NewTrack(x_cm, y_cm, z_cm, time, dx, dy, dz);
    for (const auto& cluster : fTrackHeed->GetClusters()) {
      if (fabs(cluster.z) > lTube) continue;
      if (sqrt(cluster.x * cluster.x + cluster.y * cluster.y) >= rTube) {
        continue;
      }
      nsum += cluster.electrons.size();
      fEnergyDeposit += cluster.energy;
      for (const auto& electron : cluster.electrons) {
        if (fabs(electron.z) > lTube) continue;
        if (sqrt(electron.x * electron.x + electron.y * electron.y) >= rTube) {
          continue;
        }
        analysisManager->FillH3(1, electron.z * 10, electron.x * 10, electron.y * 10);
        if (createSecondariesInGeant4) {
          double newTime = electron.t;
          if (newTime < time) {
            newTime += time;
          }
          fSecondaryParticles.emplace_back(GarfieldParticle(
              "e-", electron.e, newTime, electron.x, electron.y, electron.z, 
              electron.dx, electron.dy, electron.dz));
        }

        drift.DriftElectron(electron.x, electron.y, electron.z, electron.t);

        double xe1, ye1, ze1, te1;
        double xe2, ye2, ze2, te2;

        int status;
        drift.GetElectronEndpoint(0, xe1, ye1, ze1, te1, xe2, ye2, ze2,
                                  te2, status);

        if (0 < xe2 && xe2 < rWire) {
          xe2 += 2 * rWire;
        } else if (0 > xe2 && xe2 > -rWire) {
          xe2 += -2 * rWire;
        }
        if (0 < ye2 && ye2 < rWire) {
          ye2 += 2 * rWire;
        } else if (0 > ye2 && ye2 > -rWire) {
          ye2 += -2 * rWire;
        }

        double e2 = 0.1;
        avalanche.AvalancheElectron(xe2, ye2, ze2, te2, e2, 0, 0, 0);

        int ne = 0, ni = 0;
        avalanche.GetAvalancheSize(ne, ni);
        fAvalancheSize += ne;
      }
    }
  }
  fGain = fAvalancheSize / nsum;
}

