#ifndef G_COMPONENT_ANSYS121_H
#define G_COMPONENT_ANSYS121_H

#include "ComponentFieldMap.hh"

namespace Garfield {

/// Component for importing and interpolating two-dimensional ANSYS field maps.

class ComponentAnsys121 : public ComponentFieldMap {
 public:
  /// Constructor
  ComponentAnsys121();
  /// Destructor
  ~ComponentAnsys121() {}

  /** Import a field map.
    * \param elist name of the file containing the list of elements
    * \param nlist name of the file containing the list of nodes
    * \param mplist name of the file containing the list of materials
    * \param prnsol name of the file containing the nodal solutions
    * \param unit length unit
    */ 
  bool Initialise(const std::string& elist = "ELIST.lis",
                  const std::string& nlist = "NLIST.lis",
                  const std::string& mplist = "MPLIST.lis",
                  const std::string& prnsol = "PRNSOL.lis", 
                  const std::string& unit = "cm");
  /// Import weighting potentials.
  bool SetWeightingPotential(const std::string& prnsol, 
                             const std::string& label) {
    return SetWeightingField(prnsol, label);
  } 
  bool SetWeightingField(const std::string& prnsol, const std::string& label);
  /// Set the limits of the active region along z.
  void SetRangeZ(const double zmin, const double zmax);

 private:
  bool LoadPotentials(const std::string& prnsol, std::vector<double>& pot);
};
}

#endif
