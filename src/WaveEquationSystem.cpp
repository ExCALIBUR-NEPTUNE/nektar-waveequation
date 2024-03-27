#include "WaveEquationSystem.hpp"
#include <cstddef>

namespace Nektar {
std::string WaveEquationSystem::className =
    GetEquationSystemFactory().RegisterCreatorFunction("WaveEquationSystem",
                                                       WaveEquationSystem::create);
WaveEquationSystem::WaveEquationSystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : EquationSystem(pSession, pGraph), m_factors(),m_diff_in_arr(1), m_diff_out_arr(1), m_diff_fields(1) {

  pSession->LoadParameter("TimeStep", m_timestep);
  ASSERTL1(m_timestep > 0,
           "TimeStep must be set and be > 0 in the xml config file.");

  // m_factors[StdRegions::eFactorLambda] = 0.0;
  m_factors[StdRegions::eFactorTau] = 1.0;
  auto variables = pSession->GetVariables();
  int index = 0;
  for (auto vx : variables) {
    this->field_to_index[vx] = index;
    index++;
  }

  this->m_session->LoadParameter("Theta", this->m_theta, 0.0);
  this->m_session->LoadParameter("TimeStep", this->m_timeStep);

  ASSERTL1(this->GetFieldIndex("u0") > -1, "Could not get index for u1.");
  ASSERTL1(this->GetFieldIndex("u1") > -1, "Could not get index for u0.");
  ASSERTL1(this->GetFieldIndex("u_1") > -1, "Could not get index for u_1.");
  ASSERTL1(this->GetFieldIndex("s") > -1, "Could not get index for s.");
}


void WaveEquationSystem::v_InitObject(bool DeclareField)
{
  this->EquationSystem::v_InitObject(true); // call EquationSystem v_InitObject

  int nPts = GetNpoints();
  int nCfs = GetNcoeffs();
  int nVar = m_session->GetVariables().size();

  for (int i = 0; i < nVar; ++i) {
    this->m_physarrays.push_back(Array<OneD, NekDouble>(nPts));
    this->m_coeffarrays.push_back(Array<OneD, NekDouble>(nCfs));
    this->m_fields[i]->SetPhysArray(this->m_physarrays[i]);
    this->m_fields[i]->SetCoeffsArray(this->m_coeffarrays[i]);
  }

  // Read ICs from the file
  const int domain = 0; // if this is different to the DOMAIN in the mesh it segfaults.
  this->v_SetInitialConditions(0.0, true, domain);

  for (auto f : m_fields) {
    ASSERTL1(f->GetNpoints() > 0, "GetNpoints must return > 0");
    ASSERTL1(f->GetNcoeffs() > 0, "GetNcoeffs must return > 0");
  }
  m_laplacetmp = Array<OneD, NekDouble>(nPts);
  m_implicittmp = Array<OneD, NekDouble>(nPts);

  // Set up diffusion object
  std::string diff_type;
  m_session->LoadSolverInfo("DiffusionType", diff_type, "LDG");
  m_diffusion = GetDiffusionFactory().CreateInstance(diff_type, diff_type);
  m_diffusion->SetFluxVector(&WaveEquationSystem::GetDiffusionFluxVector, this);
  m_diffusion->InitObject(m_session, m_fields);

  // Allocate temporary arrays used in the diffusion calc
  m_diff_in_arr[0] = Array<OneD, NekDouble>(nPts);
  m_diff_out_arr[0] = Array<OneD, NekDouble>(nPts, 0.0);
  m_diff_fields[0] = m_fields[this->GetFieldIndex("u0")];
}


int WaveEquationSystem::GetFieldIndex(const std::string name) {
  ASSERTL1(this->field_to_index.count(name) > 0,
           "Could not map field name to index.");
  return this->field_to_index[name];
}

WaveEquationSystem::~WaveEquationSystem() {}

void WaveEquationSystem::v_GenerateSummary(SolverUtils::SummaryList &s) {
  EquationSystem::SessionSummary(s);
}

Array<OneD, bool> WaveEquationSystem::v_GetSystemSingularChecks() {
  auto singular_bools =
      Array<OneD, bool>(m_session->GetVariables().size(), false);
  return singular_bools;
}

void WaveEquationSystem::v_DoSolve() {
  const int u1 = this->GetFieldIndex("u1");
  const int u0 = this->GetFieldIndex("u0");
  const int u_1 = this->GetFieldIndex("u_1");
  const int s = this->GetFieldIndex("s");

  LorenzGaugeSolve(u0, u_1, s);
}

void WaveEquationSystem::SubtractMean(const int field_index) {
  auto field = this->m_fields[field_index];
  // Nektar reduces the integral across all ranks
  double integral = field->Integral();
  double volume = 1.0; // TODO this properly
  const double mean = integral / volume;
  const int nPts = GetNpoints();
  Vmath::Sadd(nPts, -mean, field->GetPhys(), 1, field->UpdatePhys(), 1);
  field->FwdTrans(field->GetPhys(), field->UpdateCoeffs());
}

void WaveEquationSystem::setTheta(const double theta) {
  ASSERTL1(0 <= theta,
           "Theta (0 = explicit, 1=implicit) must not be negative.");
  ASSERTL1(theta <= 1,
           "Theta (0 = explicit, 1=implicit) must not be greater than 1.");
  m_theta = theta;
}



/**
 * @brief Return the flux vector for the unsteady diffusion problem.
 */
void WaveEquationSystem::GetDiffusionFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscous_tensor) {
  boost::ignore_unused(in_arr);

  unsigned int nDim = q_field.size();
  unsigned int nConvectiveFields = q_field[0].size();
  unsigned int nPts = q_field[0][0].size();

  // Hard-code diff coeffs
  NekDouble d[2] = {1.0, 1.0};

  for (unsigned int j = 0; j < nDim; ++j) {
    for (unsigned int i = 0; i < nConvectiveFields; ++i) {
      Vmath::Smul(nPts, d[j], q_field[j][i], 1, viscous_tensor[j][i], 1);
    }
  }
}

void WaveEquationSystem::Laplace(
                                Array<OneD, NekDouble>& rhs,
                                const int index) {
  const int nPts = GetNpoints();

  // Use diffusion object to calculate second deriv
  // - Copy target field into temp array to work around diffusion API restrictions
  Vmath::Zero(nPts, m_diff_out_arr[0], 1);
  Vmath::Vcopy(nPts, m_fields[index]->GetPhys(), 1, m_diff_in_arr[0], 1);
  m_diffusion->Diffuse(1, m_diff_fields, m_diff_in_arr, m_diff_out_arr);
  // - Add result to RHS
  Vmath::Vadd(nPts, m_diff_out_arr[0], 1, rhs, 1, rhs, 1);

}

void WaveEquationSystem::LorenzGaugeSolve(const int field_t_index,
                                          const int field_t_minus1_index,
                                          const int source_index) {
  // copy across into shorter variable names to make sure code fits
  // on one line, more readable that way.
  const int f0     = field_t_index;
  const int f_1    = field_t_minus1_index;
  const int s      = source_index;
  const int nPts   = GetNpoints();
  const int nCfs   = GetNcoeffs();
  const double dt2 = std::pow(m_timestep, 2);

  auto f0phys  = m_fields[f0]->UpdatePhys();
  auto f_1phys = m_fields[f_1]->UpdatePhys();
  auto sphys   = m_fields[s]->GetPhys();

  Array<OneD, NekDouble> rhs(nPts, 0.0);
  Laplace(rhs, f0); // rhs = ∇² f0

  if (m_theta == 0.0) {
    // f⁺ = (2 + Δt^2 ∇²) f⁰ - f⁻ + Δt^2 s
    Vmath::Smul(nPts, dt2, rhs, 1, rhs, 1);
    // rhs = Δt^2 ∇² f0
    // Svtvp (n, a, x, _, y, _, z, _) -> z = a * x + y
    Vmath::Svtvp(nPts, 2.0, f0phys, 1, rhs, 1, rhs, 1);
    // rhs = Δt^2 ∇² f0 + 2f0
    Vmath::Svtvp(nPts, -1.0, f_1phys, 1, rhs, 1, rhs, 1);
    // rhs = Δt^2 ∇² f0 + 2f0 - f_1
    Vmath::Svtvp(nPts, dt2, sphys, 1, rhs, 1, rhs, 1);
    // rhs = Δt^2 ∇² f0 + 2f0 - f_1 + dt^2 s

    Vmath::Vcopy(nPts, f0phys, 1, f_1phys, 1);
    // f_1 -> f0 // f_1 now holds f0 (phys values)
    // Copy f_1 coefficients to f0 (no need to solve again!) ((N.B. phys values
    // copied across above)) N.B. phys values were copied above
    Vmath::Vcopy(nCfs, m_fields[f0]->GetCoeffs(), 1,
                 m_fields[f_1]->UpdateCoeffs(), 1);

    Vmath::Vcopy(nPts, rhs, 1, f0phys, 1);

    m_fields[f0]->FwdTrans(f0phys, m_fields[f0]->UpdateCoeffs());
  } else {
    // need in the form (∇² - lambda)f⁺ = rhs, where
    double lambda = 2.0 / dt2 / m_theta;

    // and currently rhs = ∇² f0
    Vmath::Smul(nPts, -2 * (1 - m_theta) / m_theta, rhs, 1, rhs, 1);
    // Svtvp (n, a, x, _, y, _, z, _) -> z = a * x + y
    Vmath::Svtvp(nPts, -2 * lambda, f0phys, 1, rhs, 1, rhs,
                 1); // rhs now holds the f0 rhs values

    // and currently rhs = -2 (lambda + (1-theta)/theta ∇²) f0

    Laplace(m_implicittmp, f_1); // m_implicittmp = ∇² f_1

    // Svtvp (n, a, x, _, y, _, z, _) -> z = a * x + y
    Vmath::Svtvp(nPts, -lambda, f_1phys, 1, m_implicittmp, 1, m_implicittmp, 1);
    // m_implicittmp = (∇² - lambda) f_1
    Vmath::Vsub(nPts, rhs, 1, m_implicittmp, 1, rhs, 1); // rhs now holds the f0 and f_1 rhs values
    // rhs = rhs - m_implicittmp
    // rhs = -2 (lambda + (1-theta)/theta ∇²) f0 - (∇² - lambda) f_1

    // Svtvp (n, a, x, _, y, _, z, _) -> z = a * x + y
    Vmath::Svtvp(nPts, -2.0 / m_theta, sphys, 1, rhs, 1, rhs, 1);
    // rhs now has the source term too
    // rhs = -2 (lambda + (1-theta)/theta ∇²) f0 - (∇² - lambda) f_1 - 2/theta * s

    // copy f0 coefficients to f_1 (no need to solve again!)
    Vmath::Vcopy(nPts, m_fields[f0]->GetPhys(), 1,
        m_fields[f_1]->UpdatePhys(), 1);
    Vmath::Vcopy(nCfs, m_fields[f0]->GetCoeffs(), 1,
        m_fields[f_1]->UpdateCoeffs(), 1);

    bool rhsAllZero = true;
    for (auto i : rhs) {
      if (i != 0.0) {
        rhsAllZero = false;
        break;
      }
    }
    // now solve f1 but store in f0
    if (!rhsAllZero) {
      m_factors[StdRegions::eFactorLambda] = lambda;
      m_fields[f0]->HelmSolve(rhs, m_fields[f0]->UpdateCoeffs(), m_factors);
      m_fields[f0]->BwdTrans(m_fields[f0]->GetCoeffs(), m_fields[f0]->UpdatePhys());
    } else {
      Vmath::Zero(nCfs, m_fields[f0]->UpdateCoeffs(), 1);
      Vmath::Zero(nPts, m_fields[f0]->UpdatePhys(), 1);
    }
  }
}

} // namespace Nektar
