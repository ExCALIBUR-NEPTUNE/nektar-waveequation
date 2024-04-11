#ifndef NEKTAR_SOLVERS_EQUATIONSYSTEMS_WAVE_EQUATION_SYSTEM_H
#define NEKTAR_SOLVERS_EQUATIONSYSTEMS_WAVE_EQUATION_SYSTEM_H

#include <map>
#include <string>

#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/EquationSystem.h>

using namespace Nektar::SolverUtils;

namespace Nektar {
class WaveEquationSystem : public EquationSystem {
public:
  std::map<std::string, int> field_to_index;

  friend class MemoryManager<WaveEquationSystem>;

  WaveEquationSystem(const LibUtilities::SessionReaderSharedPtr &pSession,
                     const SpatialDomains::MeshGraphSharedPtr &pGraph);

  /// Creates an instance of this class
  static EquationSystemSharedPtr
  create(const LibUtilities::SessionReaderSharedPtr &pSession,
         const SpatialDomains::MeshGraphSharedPtr &pGraph) {
    double timestep;
    pSession->LoadParameter("TimeStep", timestep);
    EquationSystemSharedPtr p =
        MemoryManager<WaveEquationSystem>::AllocateSharedPtr(pSession, pGraph);
    p->InitObject();
    return p;
  }
  /// Name of class
  static std::string className;

  virtual ~WaveEquationSystem();
  /**
   *  Helper function to map from field name to field indices.
   *
   *  @param name Field name
   *  @returns index (probably 0 or 1).
   */
  int GetFieldIndex(const std::string name);

  void setDtMultiplier(const double dtMultiplier);

  void setTheta(const double theta);

  double timeStep();

  void SetVolume(const double volume);

  void LorenzGaugeSolve(const int field_t_index, const int field_t_minus1_index,
                        const int source_index);

protected:
  StdRegions::ConstFactorMap m_factors;
  void v_InitObject(bool DeclareFields = true) override;
  void v_DoSolve() override;
  void v_GenerateSummary(SolverUtils::SummaryList &s) override;
  void SubtractMean(int field_index);
  std::vector<Array<OneD, NekDouble>> m_physarrays;
  std::vector<Array<OneD, NekDouble>> m_coeffarrays;

private:

  double m_timeStep;
  double m_theta;

  Array<OneD, bool> v_GetSystemSingularChecks() override;
  Array<OneD, NekDouble> m_laplacetmp;
  Array<OneD, NekDouble> m_implicittmp;

  //  std::map<int, Array<OneD, NekDouble>> m_mapIntToArray;

  // Diffusion object
  DiffusionSharedPtr m_diffusion;

  void GetDiffusionFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscous_tensor);
};
} // namespace Nektar

#endif
