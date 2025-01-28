/**
   \file PTOF/Steppers.h
   \author Tomás Aquino
   \date 09/03/2022
   \brief Objects and utilities to handle time stepping along advective and
   diffusive trajectories.
*/

#ifndef PTOF_STEPPERS_H
#define PTOF_STEPPERS_H

#include "CTRW/JumpGenerator.h"
#include "CTRW/Meta.h"
#include "CTRW/TimeGenerator.h"
#include "General/Meta.h"
#include <cstddef>
#include <type_traits>
#include <utility>

namespace ptof {

/**
   \struct CTRWStepper PTOF/Steppers.h
   "PTOF/Steppers.h"
   \brief Types to choose options of CTRW stepping mode.
*/
struct CTRWStepper {
  /**
     \struct Asynchronous PTOF/Steppers.h
     "PTOF/Steppers.h"
     \brief Step particles independtly
  */
  struct Asynchronous {
    template <typename CTRW, typename Transitions, typename SolverParameters,
              typename Time, typename Update = meta::DoNothing>
    static void evolve(CTRW &ctrw, Transitions &&transitions,
                       SolverParameters const &parameters_solvers,
                       Time new_time, Time old_time, Update &&update = {}) {
      ctrw.evolve(
          [new_time](typename CTRW::Particle const &part) {
            return part.state_new().time < new_time &&
                   !part.state_new().info.absorbed;
          },
          transitions);
      update(ctrw, new_time, old_time);
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "CTRW Stepper\n"
             "--------------------------------------------------------------\n"
             "Asynchronous stepping\n"
             "--------------------------------------------------------------\n";
      return output;
    }
  };

  /**
     \struct ParticleStep PTOF/Steppers.h
     "PTOF/Steppers.h"
     \brief Evolve CTRW particles with synchronization over particle steps.
  */
  struct ParticleStep {
    template <typename CTRW, typename Transitions, typename SolverParameters,
              typename Time, typename Update = meta::DoNothing>
    static void evolve(CTRW &ctrw, Transitions &&transitions,
                       SolverParameters const &parameters_solvers,
                       Time new_time, Time old_time, Update &&update = {}) {
      std::size_t nr_stepped = 1;
      while (nr_stepped != 0) {
        nr_stepped = ctrw.step(
            [new_time](typename CTRW::Particle const &part) {
              return part.state_new().time < new_time &&
                     !part.state_new().info.absorbed;
            },
            transitions);
        update(ctrw, new_time, old_time);
      }
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "CTRW Stepper\n"
             "--------------------------------------------------------------\n"
             "Particle-based stepping\n"
             "--------------------------------------------------------------\n";
      return output;
    }
  };

  /**
     \struct TimeStep PTOF/Steppers.h
     "PTOF/Steppers.h"
     \brief Evolve CTRW particles with synchronization over prescribed time
     step.
  */
  struct TimeStep {
    template <typename CTRW, typename Transitions, typename SolverParameters,
              typename Time, typename Update = meta::DoNothing>
    static void evolve(CTRW &ctrw, Transitions &&transitions,
                       SolverParameters const &parameters_solvers,
                       Time new_time, Time old_time, Update &&update = {}) {
      double current_time = old_time;
      while (current_time < new_time) {
        current_time += parameters_solvers.ctrw_time_step;
        ctrw.evolve(
            [current_time](typename CTRW::Particle const &part) {
              return part.state_new().time < current_time &&
                     !part.state_new().info.absorbed;
            },
            transitions);
        update(ctrw, new_time, old_time);
      }
    }
  };

  /**
       \brief Output generic information about object.
       \param output Output stream.
    */
  inline static std::ostream &info(std::ostream &output) {
    output
        << "--------------------------------------------------------------\n"
           "CTRW Stepper\n"
           "--------------------------------------------------------------\n"
           "Time-based stepping\n"
           "--------------------------------------------------------------\n";
    return output;
  }
};

struct Stepper {
  struct Euler {};
  struct RK4 {};
};

/**
   \struct Steppers_Advection_RK4_Diffusion_Euler PTOF/Steppers.h
   "PTOF/Steppers.h"
   \brief Time steppers for RK4 advection and stochastic forward Euler
   diffusion.
*/
struct Steppers_Advection_RK4_Diffusion_Euler {
  Steppers_Advection_RK4_Diffusion_Euler() = delete;

  /**
     \param params_solvers Solver parameters.
     \return Deterministic time step TimeGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename SolverParameters>
  static auto makeTimeGenerator(SolverParameters const &params_solvers) {
    if constexpr (meta::has_time_step_v<SolverParameters>)
      return ctrw::TimeGenerator_Step{params_solvers.time_step};
    return ctrw::TimeGenerator_Step{0.};
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters.
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Advection--diffusion JumpGenerator.
     \note
     - \p params_transport must define:
         - \c diff_coeff
     - If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
   */
  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto
  makeJumpGenerator(VelocityField &&velocity_field, Boundary &&boundary,
                    TransportParameters const &params_transport,
                    SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>)
      return ctrw::JumpGenerator_Add{
          ctrw::JumpGenerator_Velocity_RK4{
              std::forward<VelocityField>(velocity_field),
              params_solvers.time_step, std::forward<Boundary>(boundary)},
          ctrw::JumpGenerator_Diffusion<ParallelOption>{
              params_transport.diff_coeff, params_solvers.time_step, dim}};
    return ctrw::JumpGenerator_Add{
        ctrw::JumpGenerator_Velocity_RK4{
            std::forward<VelocityField>(velocity_field), 0.,
            std::forward<Boundary>(boundary)},
        ctrw::JumpGenerator_Diffusion<ParallelOption>{
            params_transport.diff_coeff, 0., dim}};
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure advection JumpGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto makeJumpGenerator_Advection(
      VelocityField &&velocity_field, Boundary &&boundary,
      TransportParameters const &params_transport,
      SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>)
      return ctrw::JumpGenerator_Velocity_RK4{
          std::forward<VelocityField>(velocity_field), params_solvers.time_step,
          std::forward<Boundary>(boundary)};
    return ctrw::JumpGenerator_Velocity_RK4{
        std::forward<VelocityField>(velocity_field), 0.,
        std::forward<Boundary>(boundary)};
  }

  /**
       \brief Output generic information about object.
       \param output Output stream.
    */
  inline static std::ostream &info(std::ostream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Steppers\n"
           "--------------------------------------------------------------\n"
           "Advection: Euler\n"
           "Diffusion: Stochastic Euler\n"
           "--------------------------------------------------------------\n";
    return output;
  }
};

/**
   \struct Steppers_Advection_Euler_Diffusion_Euler PTOF/Steppers.h
   "PTOF/Steppers.h"
   \brief Time steppers for forward Euler advection and stochastic forward
   Euler diffusion.
*/
struct Steppers_Advection_Euler_Diffusion_Euler {
  Steppers_Advection_Euler_Diffusion_Euler() = delete;

  /**
     \param params_solvers Solver parameters.
     \return Deterministic time step TimeGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename SolverParameters>
  static auto makeTimeGenerator(SolverParameters const &params_solvers) {
    if constexpr (meta::has_time_step_v<SolverParameters>)
      return ctrw::TimeGenerator_Step{params_solvers.time_step};
    return ctrw::TimeGenerator_Step{0.};
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters.
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Advection--diffusion JumpGenerator.
     \note
     -\p params_transport must define:
         - \c diff_coeff
     -If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto
  makeJumpGenerator(VelocityField &&velocity_field, Boundary &&boundary,
                    TransportParameters const &params_transport,
                    SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>)
      return ctrw::JumpGenerator_Add{
          ctrw::JumpGenerator_Velocity{
              std::forward<VelocityField>(velocity_field),
              params_solvers.time_step},
          ctrw::JumpGenerator_Diffusion<ParallelOption>{
              params_transport.diff_coeff, params_solvers.time_step, dim}};
    return ctrw::JumpGenerator_Add{
        ctrw::JumpGenerator_Velocity{
            std::forward<VelocityField>(velocity_field), 0.},
        ctrw::JumpGenerator_Diffusion<ParallelOption>{
            params_transport.diff_coeff, 0., dim}};
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure advection JumpGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto makeJumpGenerator_Advection(
      VelocityField &&velocity_field, Boundary &&boundary,
      TransportParameters const &params_transport,
      SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>)
      return ctrw::JumpGenerator_Velocity{
          std::forward<VelocityField>(velocity_field),
          params_solvers.time_step};
    return ctrw::JumpGenerator_Velocity{
        std::forward<VelocityField>(velocity_field), 0.};
  }

  /**
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure diffusion JumpGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename ParallelOption, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto makeJumpGenerator_Diffusion(
      Boundary &&boundary, TransportParameters const &params_transport,
      SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>)
      return ctrw::JumpGenerator_Diffusion<ParallelOption>{
          params_transport.diff_coeff, params_solvers.time_step, dim};
    return ctrw::JumpGenerator_Diffusion<ParallelOption>{
        params_transport.diff_coeff, 0., dim};
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Steppers\n"
           "--------------------------------------------------------------\n"
           "Advection: Euler\n"
           "Diffusion: Stochastic Euler\n"
           "--------------------------------------------------------------\n";
    return output;
  }
};
} // namespace ptof

#endif /* PTOF_STEPPERS_H */
