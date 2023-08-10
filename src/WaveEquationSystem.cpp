#include "WaveEquationSystem.hpp"
#include <cstddef>

namespace Nektar {
std::string WaveEquationSystem::className =
    GetEquationSystemFactory().RegisterCreatorFunction("WaveEquationSystem",
                                                       WaveEquationSystem::create);
WaveEquationSystem::WaveEquationSystem(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : EquationSystem(pSession, pGraph), m_factors() {

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
  this->EquationSystem::v_InitObject(DeclareField); // call EquationSystem v_InitObject

  std::vector<std::string> variables = m_session->GetVariables();
  // Check variables are defined.
  for (int i = 0; i < variables.size(); ++i)
  {
      ASSERTL0(variables[i] == m_session->GetVariable(i),
               "Variable '" + variables[i] + "' not defined.");
  }

  int nPts = GetNpoints();
  int nVar = variables.size();
  for (int i = 0; i < nVar; ++i) {
    this->m_physarrays.push_back(Array<OneD, NekDouble>(nPts));
    this->m_coeffarrays.push_back(Array<OneD, NekDouble>(nPts));
    this->m_fields[i]->SetPhysArray(this->m_physarrays[i]);
    this->m_fields[i]->SetCoeffsArray(this->m_coeffarrays[i]);
  }

  // Read ICs from the file
  this->v_SetInitialConditions(0.0, true, 0);

  for (auto f : m_fields) {
    ASSERTL1(f->GetNpoints() > 0, "GetNpoints must return > 0");
  }

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

void WaveEquationSystem::Laplace(Array<OneD, NekDouble>& tmp,
                                Array<OneD, NekDouble>& rhs,
                                const int index) {
  const int nPts = tmp.GetCount();
  Vmath::Zero(nPts, tmp, 1);
  auto foo = m_fields[index]->GetPhys();
  m_fields[index]->PhysDeriv(MultiRegions::eX, foo, tmp);
  m_fields[index]->PhysDeriv(MultiRegions::eX, m_fields[index]->GetPhys(),
                          tmp);
  m_fields[index]->PhysDeriv(MultiRegions::eX, tmp, tmp);// tmp = ∇x² f
  Vmath::Vadd(nPts, tmp, 1, rhs, 1, rhs, 1); // rhs = rhs + tmp = rhs + ∇x² f
  Vmath::Zero(nPts, tmp, 1);
  m_fields[index]->PhysDeriv(MultiRegions::eY, m_fields[index]->GetPhys(),
                          tmp);
  m_fields[index]->PhysDeriv(MultiRegions::eY, tmp, tmp);// tmp = ∇y² f
  Vmath::Vadd(nPts, tmp, 1, tmp, 1, rhs, 1); // rhs = rhs + tmp = rhs + ∇y² f
  // rhs = ∇² f
}

void WaveEquationSystem::LorenzGaugeSolve(const int field_t_index,
                                          const int field_t_minus1_index,
                                          const int source_index) {
  // copy across into shorter variable names to make sure code fits
  // on one line, more readable that way.
  const int f0 = field_t_index;
  const int f_1 = field_t_minus1_index;
  const int s = source_index;
  const int nPts = GetNpoints();
  const double dt2 = std::pow(m_timestep, 2);

  auto f0phys = m_fields[f0]->UpdatePhys();
  auto f_1phys = m_fields[f_1]->UpdatePhys();
  auto sphys = m_fields[s]->GetPhys();

  Array<OneD, NekDouble> tmp1(nPts, 0.0);
  Array<OneD, NekDouble> rhs(nPts, 0.0);
  Laplace(tmp1, rhs, f0); // rhs = ∇² f0

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
    Vmath::Vcopy(nPts, m_fields[f0]->GetCoeffs(), 1,
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

    Array<OneD, NekDouble> tmp2(nPts, 0.0);
    Laplace(tmp1, tmp2, f_1); // tmp2 = ∇² f_1

    // Svtvp (n, a, x, _, y, _, z, _) -> z = a * x + y
    Vmath::Svtvp(nPts, -lambda, f_1phys, 1, tmp2, 1, tmp2, 1);
    // tmp2 = (∇² - lambda) f_1
    Vmath::Vsub(nPts, rhs, 1, tmp2, 1, rhs, 1); // rhs now holds the f0 and f_1 rhs values
    // rhs = rhs - tmp2
    // rhs = -2 (lambda + (1-theta)/theta ∇²) f0 - (∇² - lambda) f_1

    // Svtvp (n, a, x, _, y, _, z, _) -> z = a * x + y
    Vmath::Svtvp(nPts, -2.0 / m_theta, sphys, 1, rhs, 1, rhs, 1);
    // rhs now has the source term too
    // rhs = -2 (lambda + (1-theta)/theta ∇²) f0 - (∇² - lambda) f_1 - 2/theta * s

    // copy f_1 coefficients to f0 (no need to solve again!)
    Vmath::Vcopy(nPts, m_fields[f0]->GetPhys(), 1,
        m_fields[f_1]->UpdatePhys(), 1);
    Vmath::Vcopy(nPts, m_fields[f0]->GetCoeffs(), 1,
        m_fields[f_1]->UpdateCoeffs(), 1);

    bool rhsAllZero = true;
    for (auto i : rhs) {
      if (i != 0.0) {
        rhsAllZero = false;
        break;
      }
    }

    if (!rhsAllZero) {
      m_factors[StdRegions::eFactorLambda] = lambda;
      m_fields[f0]->HelmSolve(rhs, m_fields[f0]->UpdateCoeffs(), m_factors);
      m_fields[f0]->BwdTrans(m_fields[f0]->GetCoeffs(), m_fields[f0]->UpdatePhys());
    } else {
      Vmath::Zero(nPts, m_fields[f0]->UpdateCoeffs(), 1);
      Vmath::Zero(nPts, m_fields[f0]->UpdatePhys(), 1);
    }
  }

  m_fields[f0]->SetPhysState(true); // don't need this
}

} // namespace Nektar
