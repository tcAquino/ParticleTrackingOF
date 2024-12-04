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
#include <cstddef>
#include <type_traits>
#include <utility>

namespace ptof {
/**
   \struct Steppers_Advection_RK4_Diffusion_Euler PTOF/Steppers.h
   "PTOF/Steppers.h"
   \brief Time steppers for RK4 advection and stochastic forward Euler
   diffusion.
*/
struct Steppers_Advection_RK4_Diffusion_Euler {
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
     -# \c diff_coeff
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
};

/**
   \struct Steppers_Advection_Euler_Diffusion_Euler PTOF/Steppers.h
   "PTOF/Steppers.h"
   \brief Time steppers for forward Euler advection and stochastic forward Euler
   diffusion.
*/
struct Steppers_Advection_Euler_Diffusion_Euler {
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
     -# \c diff_coeff
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
};
} // namespace ptof

#endif /* PTOF_STEPPERS_H */
