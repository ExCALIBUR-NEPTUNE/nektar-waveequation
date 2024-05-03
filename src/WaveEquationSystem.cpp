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
  this->EquationSystem::v_InitObject(true); // call EquationSystem v_InitObject

  int nPts = GetNpoints();
  int nCfs = GetNcoeffs();
  int nVar = m_session->GetVariables().size();

  /*
  for (int i = 0; i < nVar; ++i) {
    this->m_physarrays.push_back(Array<OneD, NekDouble>(nPts));
    this->m_coeffarrays.push_back(Array<OneD, NekDouble>(nCfs));
    this->m_fields[i]->SetPhysArray(this->m_physarrays[i]);
    this->m_fields[i]->SetCoeffsArray(this->m_coeffarrays[i]);
  }
  */

  // Read ICs from the file
  const int domain = 0; // if this is different to the DOMAIN in the mesh it segfaults.
  this->SetInitialConditions(0.0, true, domain);

  for (auto f : m_fields) {
    ASSERTL1(f->GetNpoints() > 0, "GetNpoints must return > 0");
    ASSERTL1(f->GetNcoeffs() > 0, "GetNcoeffs must return > 0");
  }

  // Set up diffusion object
  std::string diff_type;
  m_session->LoadSolverInfo("DiffusionType", diff_type, "LDG");
  m_diffusion = GetDiffusionFactory().CreateInstance(diff_type, diff_type);
  m_diffusion->SetFluxVector(&WaveEquationSystem::GetDiffusionFluxVector, this);
  m_diffusion->InitObject(m_session, m_fields);
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
  auto &f0coeff  = m_fields[f0]->UpdateCoeffs();
  auto &f_1coeff = m_fields[f_1]->UpdateCoeffs();
  auto scoeff   = m_fields[s]->GetCoeffs();

  Array<OneD, NekDouble> rhs(nCfs, 0.0), tmp(nCfs, 0.0), tmp2(nCfs, 0.0);

  // Apply Laplacian matrix op -> tmp
  MultiRegions::GlobalMatrixKey laplacianKey(StdRegions::eLaplacian);
  MultiRegions::GlobalMatrixKey massKey(StdRegions::eMass);

  if (m_theta == 0.0) {
    // Evaluate M^{-1} * L * u
    m_fields[f0]->GeneralMatrixOp(laplacianKey, f0coeff, tmp);
    m_fields[f0]->MultiplyByInvMassMatrix(tmp, tmp2);

    // Temporary copy for f_0 to transfer to f_{-1}
    Vmath::Vcopy(nCfs, f0coeff, 1, tmp, 1);

    // Central difference timestepping
    for (int i = 0; i < nCfs; ++i)
    {
      f0coeff[i] = 2 * f0coeff[i] - dt2 * tmp2[i] - f_1coeff[i] + dt2 * scoeff[i];
    }

    // Update f_{-1}
    Vmath::Vcopy(nCfs, tmp, 1, f_1coeff, 1);

    // backward transform -- not really necessary
    m_fields[f0]->BwdTrans(f0coeff, f0phys);
    m_fields[f_1]->BwdTrans(f_1coeff, f_1phys);
  } else if (m_theta == 1.0) {
    double lambda = 1.0 / dt2;

    for (int i = 0; i < nCfs; ++i)
    {
        // This is negative, because HelmSolve will negate the input to be
        // consistent with the Helmholtz equation definition.
        tmp2[i] = -lambda * (2 * f0coeff[i] - f_1coeff[i]);
    }

    // Evaluate mass matrix action
    m_fields[f0]->GeneralMatrixOp(massKey, tmp2, tmp);

    // Zero storage
    Vmath::Zero(nCfs, tmp2, 1);

    for (int i = 0; i < nCfs; ++i)
    {
        // This is negative, because HelmSolve will negate the input to be
        // consistent with the Helmholtz equation definition.
        tmp[i] -= scoeff[i];
    }

    m_factors[StdRegions::eFactorLambda] = lambda;

    m_fields[f0]->HelmSolve(tmp, tmp2, m_factors, StdRegions::NullVarCoeffMap,
                            MultiRegions::NullVarFactorsMap,
                            NullNekDouble1DArray, false);

    // Rotate storage
    Vmath::Vcopy(nCfs, f0coeff, 1, f_1coeff, 1);
    Vmath::Vcopy(nCfs, tmp2, 1, f0coeff, 1);

    m_fields[f0]->BwdTrans(f0coeff, f0phys);
    m_fields[f_1]->BwdTrans(f_1coeff, f_1phys);

  } else {
    // need in the form (∇² - lambda)f⁺ = rhs, where
    double lambda = 2.0 / dt2 / m_theta;

    for (int i = 0; i < nCfs; ++i)
    {
        // This is negative, because HelmSolve will negate the input to be
        // consistent with the Helmholtz equation definition.
        tmp2[i] = -lambda * (2 * f0coeff[i] - f_1coeff[i]);
    }

    m_fields[f0]->GeneralMatrixOp(massKey, tmp2, rhs);

    for (int i = 0; i < nCfs; ++i)
    {
        tmp2[i] = -(2 * (1 - m_theta) / m_theta * f0coeff[i] + f_1coeff[i]);
    }

    // zero tmp
    Vmath::Zero(nCfs, tmp, 1);

    // now do diffusion operator
    m_fields[f0]->GeneralMatrixOp(laplacianKey, tmp2, tmp);

    for (int i = 0; i < nCfs; ++i)
    {
        // copy the second term and sources into rhs
        rhs[i] += tmp[i] + 2.0 / m_theta * scoeff[i]; // should have minus sign?
    }

    // Zero storage
    Vmath::Zero(nCfs, tmp2, 1);

    m_fields[f0]->HelmSolve(rhs, tmp2, m_factors, StdRegions::NullVarCoeffMap,
                            MultiRegions::NullVarFactorsMap,
                            NullNekDouble1DArray, false);

    // Rotate storage
    Vmath::Vcopy(nCfs, f0coeff, 1, f_1coeff, 1);
    Vmath::Vcopy(nCfs, tmp2, 1, f0coeff, 1);

    m_fields[f0]->BwdTrans(f0coeff, f0phys);
    m_fields[f_1]->BwdTrans(f_1coeff, f_1phys);
  }
}

} // namespace Nektar
