#ifndef G_TRACK_BICHSEL_H
#define G_TRACK_BICHSEL_H

#include <array>

#include "FundamentalConstants.hh"
#include "Track.hh"

namespace Garfield {

/// Generate tracks using differential cross-sections
/// for silicon computed by Hans Bichsel.
/// References:
///   - H. Bichsel, Rev. Mod. Phys. 60 (1988), 663-699
///   - https://faculty.washington.edu/hbichsel/

class TrackBichsel : public Track {
 public:
  struct Cluster {
    double x, y, z, t;
    double energy;
  };

  /// Default constructor
  TrackBichsel() : TrackBichsel(nullptr) {}
  /// Constructor
  TrackBichsel(Sensor* sensor);
  /// Destructor
  virtual ~TrackBichsel() {}

  bool NewTrack(const double x0, const double y0, const double z0,
                const double t0, const double dx0, const double dy0,
                const double dz0) override;
  bool GetCluster(double& xc, double& yc, double& zc, double& tc, int& nc, 
                  double& ec, double& extra) override;
  const std::vector<Cluster>& GetClusters() const { return m_clusters; }

  double GetClusterDensity() override;
  double GetStoppingPower() override;

  bool Initialise();
  bool ComputeCrossSection();

 private:
  constexpr static size_t NEnergyBins = 1250;
  std::array<double, NEnergyBins + 1> m_E;

  /// Optical oscillator strength density.
  std::array<double, NEnergyBins> m_dfdE;
  /// Real part of the dielectric function. 
  std::array<double, NEnergyBins> m_eps1;
  /// Imaginary part of the dielectric function. 
  std::array<double, NEnergyBins> m_eps2;
  /// Integral over the generalised oscillator strength density.
  std::array<double, NEnergyBins> m_int;
  /// Lower limit of the integral over the GOS.
  std::array<double, NEnergyBins> m_k1;

  constexpr static size_t NCdfBins = 10000;
  std::array<double, NCdfBins> m_tab;

  /// Density of silicon.
  double m_density;
  /// Conversion from optical loss function to oscillator strength density.
  double m_conv = 0.0092456;
  
  bool m_initialised = false;

  /// Inverse mean free path [cm-1].
  double m_imfp = 0.;
  /// Stopping power [eV/cm].
  double m_dEdx = 0.;
 
  /// Particle speed
  double m_speed = SpeedOfLight;

  std::vector<Cluster> m_clusters;
  size_t m_cluster = 0;
};
}

#endif
