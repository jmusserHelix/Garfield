#include <iostream>

#include "Garfield/FundamentalConstants.hh"
#include "Garfield/GarfieldConstants.hh"
#include "Garfield/Random.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackSimple.hh"

namespace Garfield {

TrackSimple::TrackSimple(Sensor* sensor) : Track("Simple") {
  m_sensor = sensor;
}

void TrackSimple::SetClusterDensity(const double d) {
  if (d < Small) {
    std::cerr << m_className << "::SetClusterDensity:\n"
              << "    Cluster density (number of clusters per cm)"
              << " must be positive.\n";
    return;
  }
  m_mfp = 1. / d;
}

double TrackSimple::GetClusterDensity() { return 1. / m_mfp; }

void TrackSimple::SetStoppingPower(const double dedx) {
  if (dedx < Small) {
    std::cerr << m_className << "::SetStoppingPower:\n"
              << "    Stopping power (average energy loss [eV] per cm)"
              << " must be positive.\n";
    return;
  }
  m_eloss = dedx;
}

double TrackSimple::GetStoppingPower() { return m_eloss; }

bool TrackSimple::NewTrack(const double x0, const double y0, const double z0,
                           const double t0, const double dx0, const double dy0,
                           const double dz0) {

  m_clusters.clear();
  m_cluster = 0;
  // Make sure the sensor is defined.
  if (!m_sensor) {
    std::cerr << m_className << "::NewTrack: Sensor is not defined.\n";
    return false;
  }

  // Make sure we are inside a medium.
  Medium* medium = m_sensor->GetMedium(x0, y0, z0);
  if (!medium) {
    std::cerr << m_className << "::NewTrack: No medium at initial position.\n";
    return false;
  }

  double x = x0;
  double y = y0;
  double z = z0;
  double t = t0;

  // Normalise the direction.
  double dx = dx0;
  double dy = dy0;
  double dz = dz0;
  const double d = sqrt(dx * dx + dy * dy + dz * dz);
  if (d < Small) {
    // Choose a random direction.
    RndmDirection(dx, dy, dz);
  } else {
    const double scale = 1. / d;
    dx *= scale;
    dy *= scale;
    dz *= scale;
  }
  bool ok = true;
  while (ok) {
    if (m_useEqualSpacing) {
      x += dx * m_mfp;
      y += dy * m_mfp;
      z += dz * m_mfp;
    } else {
      const double step = -m_mfp * log(RndmUniformPos());
      x += dx * step;
      y += dy * step;
      z += dz * step;
    }

    medium = m_sensor->GetMedium(x, y, z);
    if (!medium) {
      if (m_debug) {
        std::cout << m_className << "::NewTrack: Particle left the medium.\n";
      }
      break;
    }
    Cluster cluster;
    cluster.x = x;
    cluster.y = y;
    cluster.z = z;
    cluster.t = t;
    cluster.energy = m_eloss * m_mfp;
    m_clusters.push_back(std::move(cluster));
  }
  m_cluster = 0;
  return true;
}

bool TrackSimple::GetCluster(double& xc, double& yc, double& zc, double& tc,
                             int& ne, double& ec, double& extra) {
  xc = yc = zc = tc = ec = extra = 0.;
  ne = 0;
  if (m_clusters.empty() || m_cluster >= m_clusters.size()) return false;
  const auto& cluster = m_clusters[m_cluster];
  xc = cluster.x;
  yc = cluster.y;
  zc = cluster.z;
  tc = cluster.t;
  ec = cluster.energy;
  ne = 1; 
  ++m_cluster;
  return true;
}

}
