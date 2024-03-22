#ifndef G_TRACK_ELECTRON
#define G_TRACK_ELECTRON

#include <string>
#include <vector>

#include "Track.hh"

namespace Garfield {

/// [WIP] Ionization calculation based on MIP program (S. Biagi). 

class TrackElectron : public Track {
 public:
  struct Cluster {
    double x, y, z, t;
    double esec;
  };

  /// Constructor
  TrackElectron();
  /// Destructor
  virtual ~TrackElectron() {}

  void SetParticle(const std::string& particle) override;

  bool NewTrack(const double x0, const double y0, const double z0,
                const double t0, const double dx0, const double dy0,
                const double dz0) override;

  bool GetCluster(double& xc, double& yc, double& zc, double& tc, int& nc,
                  double& ec, double& extra) override;
  const std::vector<Cluster>& GetClusters() const { return m_clusters; }

  double GetClusterDensity() override;
  double GetStoppingPower() override;

 private:

  struct Parameters {
    // Dipole moment
    double m2;
    // Constant in ionisation cross-section
    double cIon;
    // Density correction term
    double x0;
    double x1;
    double cDens;
    double aDens;
    double mDens;
    // Opal-Beaty-Peterson splitting factor
    double wSplit;
    // Ionisation threshold
    double ethr;
  };

  std::vector<Cluster> m_clusters;
  size_t m_cluster = 0;

  // Mean free path
  double m_mfp = 0.;
  // Stopping power
  double m_dedx = 0.;

  static bool Setup(Medium* gas, std::vector<Parameters>& par,
                    std::vector<double>& frac);
  static bool Update(const double density, const double beta2,
                     const std::vector<Parameters>& par,
                     const std::vector<double>& frac, 
                     std::vector<double>& prob, double& mfp, double& dedx);
  static double Delta(const double x, const Parameters& par);
  static double Esec(const double e0, const Parameters& par);
};
}

#endif
