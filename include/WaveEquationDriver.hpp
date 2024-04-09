#ifndef __LORENZBORIS_H_
#define __LORENZBORIS_H_

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/Timer.h>

#include <SolverUtils/Driver.h>
#include <SolverUtils/EquationSystem.h>

#include <WaveEquationSystem.hpp>

#include <functional>
#include <memory>
#include <vector>

using namespace Nektar;
using namespace Nektar::SolverUtils;

/// Forward declaration
template <typename T> class WaveEquationDriver;

/**
 *  This is the class that sets up the 2D3V EM PIC
 *  simulation and contains the main loop.
 *
 */
template <typename T> class WaveEquationDriver : public WaveEquationSystem {
private:
  LibUtilities::SessionReaderSharedPtr session;
  SpatialDomains::MeshGraphSharedPtr graph;

  int rank;

public:
  /// The number of time steps in the main loop.
  int m_numTimeSteps;
  /// The current time step of the simulation.
  int m_timeStep;
  /// The parameter that controls implicitness (0 = explicit, 1 = implicit)
  double m_theta;
  /// the current time step number
  int m_currentTimeStep;

  /**
   *  Create new simulation instance using a nektar++ session. The parameters
   *  for the simulation are read from the nektar+ input file.
   *
   *  @param session Nektar++ session object.
   *  @param graph Nektar++ MeshGraph instance.
   */
  WaveEquationDriver(LibUtilities::SessionReaderSharedPtr session,
                     SpatialDomains::MeshGraphSharedPtr graph)
      : WaveEquationSystem(session, graph), session(session), graph(graph) {
    this->session->LoadParameter("NumTimeSteps", this->m_numTimeSteps, 1024);
  };


  void v_InitObject(bool DeclareField) override
  {
    WaveEquationSystem::v_InitObject(DeclareField);
  }


  /**
   *  Run the simulation.
   */
  inline void run() {
    std::cout << "Calling WaveEquationDriver::run" << std::endl;
    std::cout << " ... for " << m_numTimeSteps << " steps." << std::endl;
    int u0 = this->GetFieldIndex("u0");
    int u_1 = this->GetFieldIndex("u_1");
    int s = this->GetFieldIndex("s");
    int checkpoint = 0;
    for (int stepx = 1; stepx <= this->m_numTimeSteps; stepx++) {
      this->m_currentTimeStep = stepx;
      // These 5 lines perform the simulation timestep.
      this->LorenzGaugeSolve(u0, u_1, s);
      if (stepx % m_checksteps == 0) {
        checkpoint++;
        this->Checkpoint_Output(checkpoint);
      }
    } // MAIN LOOP END
    std::cout << "Called WaveEquationDriver::run" << std::endl;
  }

  /**
   * Finalise the simulation, i.e. close output files and free objects.
   */
  inline void finalise() {
  }
};

#endif
