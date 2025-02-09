#ifndef G_COMPONENT_CST_H
#define G_COMPONENT_CST_H

#include <map>

#include "ComponentFieldMap.hh"

namespace Garfield {

/// Component for importing and interpolating field maps from CST.
/// This interface assumes a certain format of the ascii files
/// Please find the tools to extract the field data from CST
/// in the correct way here: http://www.desy.de/~zenker/garfieldpp.html

class ComponentCST : public ComponentFieldMap {
 public:
  /// Constructor
  ComponentCST();
  /// Destructor
  ~ComponentCST() {}

  void ShiftComponent(const double xShift, const double yShift,
                      const double zShift);

  void GetNumberOfMeshLines(unsigned int& nx, unsigned int& ny,
                            unsigned int& nz) const;
  size_t GetNumberOfElements() const override { return m_nElements; }
  bool GetElement(const size_t i, size_t& mat, bool& drift,
                  std::vector<size_t>& nodes) const override;
  size_t GetNumberOfNodes() const override { return m_nNodes; }
  bool GetNode(const size_t i, double& x, double& y, double& z) const override; 

  void GetElementBoundaries(unsigned int element, double& xmin, double& xmax,
                            double& ymin, double& ymax, double& zmin,
                            double& zmax) const;

  Medium* GetMedium(const double x, const double y, const double z) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, Medium*& m, int& status) override;
  void ElectricField(const double x, const double y, const double z, double& ex,
                     double& ey, double& ez, double& v, Medium*& m,
                     int& status) override;
  void WeightingField(const double x, const double y, const double z,
                      double& wx, double& wy, double& wz,
                      const std::string& label) override;

  double WeightingPotential(const double x, const double y, const double z,
                            const std::string& label) override;
  /**
   * Deprecated version of the interface based on text file import of field
   * data.
   * \param elist Information about the element material of mesh cells. Each
   * line contains the element number
   * and the material index:
   * \code
   * 0  3
   * ...
   * \endcode
   * \param nlist Information about the mesh like this:
   * \code
   *  xmax 136 ymax 79 zmax 425
   *  x−l i n e s
   * 0
   * 8 . 9 2 8 5 7 e −07
   * 1 . 7 8 5 7 1 e −06
   * ...
   * y−l i n e s
   * 0
   * 8 . 9 2 8 5 7 e −07
   * 1 . 7 8 5 7 1 e −06
   * ...
   * z−l i n e s
   * 0.0027
   * 0.00270674
   * ...
   * \endcode
   * \param mplist Information about material properties used in the simulation:
   * \code
   *  Materials 4
   *  Material 1 PERX 1 . 0 0 0 0 0 0
   *  Material 2 RSVX 0 . 0 0 0 0 0 0 PERX 0 . 1 0 0 0 0 0 0E+11
   *  Material 3 PERX 3 . 5 0 0 0 0 0
   *  Material 4 PERX 4 . 8 0 0 0 0 0
   *  \endcode
   *  \param prnsol Information about the node potentials. Each line contains
   * the node number and the potential:
   *  \code
   *  0 1000.00
   *  ...
   *  \endcode
   * \param unit The units used in the nlist input file
   */
  bool Initialise(std::string elist, std::string nlist, std::string mplist,
                  std::string prnsol, std::string unit = "cm");
  /**
   * Import of field data based on binary files.
   * See http://www.desy.de/~zenker/garfieldpp.html to get information about the
   * binary files export from CST.
   * \param dataFile The binary file containing the field data exported from
   * CST.
   * \param unit The units used in the binary file. They are not necessarily
   * equal to CST units.
   */

  bool Initialise(std::string dataFile, std::string unit = "cm");
  /**
   * Initialise a weighting field.
   * This function can handle the deprecated text based file format (see
   * Initialise( std::string elist...) for
   * the expected file format, which is similar to prnsol.
   * It also also handles binary files including the weighting field.
   * \param prnsol The input file (binary/text file)
   * \param label The name of the weighting field to be added. If a weighting
   * field with same name already exist it is replaced by the new one.
   * \param isBinary Depending on the file type you use, adapt this switch.
   *
   */
  bool SetWeightingField(std::string prnsol, std::string label,
                         bool isBinary = true);

  // Range

  void SetRangeZ(const double zmin, const double zmax);
  /**
   * Use these functions to disable a certain field component.
   * Is a field component is disabled ElectricField and WeightingField will
   * return 0 for this component.
   * This is useful if you want to have calculated global field distortions and
   * you want to
   * add the field of a GEM. If you would simply add both components the field
   * component
   * in drift direction would be added twice!
   */
  void DisableXField() { disableFieldComponent[0] = true; }
  void DisableYField() { disableFieldComponent[1] = true; }
  void DisableZField() { disableFieldComponent[2] = true; }
  /**
   * If you calculate the electric field component in \f$x\f$ direction along a
   * line in x direction
   * this field component will be constant inside mesh elements (by
   * construction). This can be observed
   * by plotting \f$E_x\f$ in \f$x\f$ direction. If you plot \f$E_x\f$  in y
   * direction the field will
   * be smooth (also by construction). Yuri Piadyk proposed also to shape the
   * electric field.
   * This is done as follows. The field component calculated as described above
   * is assumed to appear
   * in the center of the mesh element.
   * ________________________
   * |    M1  P |     M2     | x direction
   * |----x---x-|-----x------|-->
   * |          |            |
   * |          |            |
   * |__________|____________|
   *
   *  element 1   element 2
   *
   * Lets consider only the \f$x\f$ direction and we want to calculate
   * \f$E_x(P)\f$. The field in the
   * center of the element containing \f$P\f$ is \f$E_x(M_1) = E_1\f$. Without
   * shaping it is \f$E_1\f$ along the
   * \f$x\f$ direction in everywhere in element 1.
   * The idea of the shaping is to do a linear interpolation of the \f$E_x\f$
   * between the field \f$E_1\f$
   * and \f$E_x(M_2)=E_2\f$. This results in a smooth electric field \f$E_x\f$
   * also in \f$x\f$ direction.
   * If P would be left from \f$M_1\f$ the field in the left neighboring element
   * would be considered.
   * In addition it is also checked if the material in both elements used for
   * the interpolation is the same.
   * Else no interpolation is done.
   * \remark This shaping gives you a nice and smooth field, but you introduce
   * additional information.
   * This information is not coming from the CST simulation, but from the
   * assumption that the field between
   * elements changes in a linear way, which might be wrong! So you might
   * consider to increase the number
   * of mesh cells used in the simulation rather than using this smoothing.
   */
  void EnableShaping() { doShaping = true; }
  void DisableShaping() { doShaping = false; }
  /**
   * Calculate the element index from the position in the x/y/z position vectors
   * (m_xlines, m_ylines, m_zlines).
   * This is public since it is used in ViewFEMesh::DrawCST.
   */
  int Index2Element(const unsigned int i, const unsigned int j,
                    const unsigned int k) const;
  /**
   * Find the positions in the x/y/z position vectors (m_xlines, m_ylines,
   * m_zlines) for a given point.
   * The operator used for the comparison is <=. Therefore, the last entry in
   * the vector will never be
   * returned for a point inside the mesh.
   */
  bool Coordinate2Index(const double x, const double y, const double z,
                        unsigned int& i, unsigned int& j, unsigned int& k) const;

 protected:
  void SetRange() override;

  double GetElementVolume(const size_t i) const override;
  void GetAspectRatio(const size_t i, 
                      double& dmin, double& dmax) const override;

  /**
   * Calculate the index in the vectors m_xlines, m_ylines, m_zlines, which is
   * before the given coordinate.
   * \remark x, y, z need to be mapped before using this function!
   * \param x The x coordinate mapped to the basic cell.
   * \param y The y coordinate mapped to the basic cell.
   * \param z The z coordinate mapped to the basic cell.
   * \param i The index of the m_xlines vector, where m_xlines.at(i) < x <
   * m_xlines.at(i+1).
   * \param j The index of the m_ylines vector, where m_ylines.at(j) < y <
   * m_ylines.at(j+1).
   * \param k The index of the m_zlines vector, where m_zlines.at(k) < z <
   * m_zlines.at(k+1).
   * \param position_mapped The calculated mapped position (x,y,z) -> (x_mapped,
   * y_mapped, z_mapped)
   * \param mirrored Information if x, y, or z direction is mirrored.
   */
  bool Coordinate2Index(const double x, const double y, const double z,
                        unsigned int& i, unsigned int& j, unsigned int& k,
                        double* position_mapped, bool* mirrored) const;

 private:
  std::vector<double> m_xlines;  ///< x positions used in the CST mesh
  std::vector<double> m_ylines;  ///< y positions used in the CST mesh
  std::vector<double> m_zlines;  ///< z positions used in the CST mesh
  /// Potentials resulting from the CST simulation.
  std::vector<float> m_potential;  
  /// Map of weighting field potentials
  std::map<std::string, std::vector<float> > m_weightingFields;  
  /// Material id for each element (unsigned char since it uses only 1 byte)
  std::vector<unsigned char> m_elementMaterial;  

  unsigned int m_nx = 0;  ///< Number of mesh lines in x direction
  unsigned int m_ny = 0;  ///< Number of mesh lines in y direction
  unsigned int m_nz = 0;  ///< Number of mesh lines in z direction
  size_t m_nElements = 0; ///< Number of elements
  size_t m_nNodes = 0;    ///< Number of nodes
  // If true x,y,z fields of this component are disabled (e=0 V/cm).
  bool disableFieldComponent[3] = {false, false, false};
  bool doShaping = false;

  void ElectricFieldBinary(const double x, const double y, const double z,
                           double& ex, double& ey, double& ez, double& v,
                           Medium*& m, int& status,
                           const bool calculatePotential = false) const;
  float GetFieldComponent(const unsigned int i, const unsigned int j,
                          const unsigned int k, const double rx,
                          const double ry, const double rz,
                          const char component,
                          const std::vector<float>& potentials) const;

  float GetPotential(const unsigned int i, const unsigned int j,
                     const unsigned int k, 
                     const double rx, const double ry, const double rz, 
                     const std::vector<float>& potentials) const;

  void ShapeField(float& ex, float& ey, float& ez, const double rx,
                  const double ry, const double rz, const unsigned int i,
                  const unsigned int j, const unsigned int k,
                  const std::vector<float>& potentials) const;

  /* Calculate the index (i,j,k) along x,y,z direction of the given element.
   * i,j,k start at 0 and reach at maximum
   * m_xlines-1,m_ylines-1,m_zlines-1
   */
  void Element2Index(const size_t element, unsigned int& i, unsigned int& j,
                     unsigned int& k) const;

  int Index2Node(const unsigned int i, const unsigned int j,
                 const unsigned int k) const;

  void Node2Index(const size_t node, unsigned int& i, unsigned int& j,
                  unsigned int& k) const;
};

}

#endif
