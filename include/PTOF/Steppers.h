/**
 \file PTOF/Steppers.h
 \author Tomás Aquino
 \date 09/03/2022
*/

#ifndef PTOF_STEPPERS_H
#define PTOF_STEPPERS_H

#include <cstddef>
#include <type_traits>
#include <utility>
#include "CTRW/JumpGenerator.h"
#include "CTRW/TimeGenerator.h"
#include "General/Useful.h"

namespace ptof
{
  /** \struct Steppers_Advection_RK4_Diffusion_Euler PTOF/Steppers.h "PTOF/Steppers.h"
   * \brief Time steppers for RK4 advection and stochastic forward Euler diffusion. */
  struct Steppers_Advection_RK4_Diffusion_Euler
  {
    /**
     \param solver_params Solver parameters.
     \return Deterministic time step TimeGenerator.
     \note
     \p solver_params must define:
     - time_step
    */
    template <typename SolverParameters>
    static auto makeTimeGenerator(SolverParameters const& solver_params)
    {
      if constexpr (useful::has_time_step<SolverParameters>::value)
        return ctrw::TimeGenerator_Step{ solver_params.time_step };
      return ctrw::TimeGenerator_Step{ 0. };
    }
    
    /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters.
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Advection--diffusion JumpGenerator.
     \note
     \p params_transport must define:
     - diff_coeff
     
     \p params_solvers must define:
     - time_step
     */
    template
    <typename VelocityField,
    typename Boundary,
    typename TransportParameters,
    typename SolverParameters>
    static auto makeJumpGenerator
    (VelocityField&& velocity_field,
     Boundary&& boundary,
     TransportParameters const& params_transport,
     SolverParameters const& params_solvers,
     std::size_t dim)
    {
      if constexpr (useful::has_time_step<SolverParameters>::value)
        return
          ctrw::JumpGenerator_Add{
            ctrw::JumpGenerator_Velocity_RK4{
              std::forward<VelocityField>(velocity_field),
              params_solvers.time_step,
              std::forward<Boundary>(boundary) },
            ctrw::JumpGenerator_Diffusion{
              params_transport.diff_coeff,
              params_solvers.time_step,
              dim } };
      return
        ctrw::JumpGenerator_Add{
          ctrw::JumpGenerator_Velocity_RK4{
            std::forward<VelocityField>(velocity_field),
            0.,
            std::forward<Boundary>(boundary) },
          ctrw::JumpGenerator_Diffusion{
            params_transport.diff_coeff,
            0.,
            dim } };
    }
    
    /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure advection JumpGenerator.
     \note
     \p params_solvers must define:
     - time_step
    */
    template
    <typename VelocityField,
    typename Boundary,
    typename TransportParameters,
    typename SolverParameters>
    static auto makeJumpGenerator_Advection
    (VelocityField&& velocity_field,
     Boundary&& boundary,
     TransportParameters const& params_transport,
     SolverParameters const& params_solvers,
     std::size_t dim)
    {
      if constexpr (useful::has_time_step<SolverParameters>::value)
        return
          ctrw::JumpGenerator_Velocity_RK4{
            std::forward<VelocityField>(velocity_field),
            params_solvers.time_step,
            std::forward<Boundary>(boundary) };
      return
        ctrw::JumpGenerator_Velocity_RK4{
          std::forward<VelocityField>(velocity_field),
          0.,
          std::forward<Boundary>(boundary) };
    }
  };
  
  /** \struct Steppers_Advection_Euler_Diffusion_Euler PTOF/Steppers.h "PTOF/Steppers.h"
   \brief Time steppers for forward Euler advection and stochastic forward Euler diffusion. */
  struct Steppers_Advection_Euler_Diffusion_Euler
  {
    /**
     \param solver_params Solver parameters.
     \return Deterministic time step TimeGenerator.
     \note
     \p solver_params must define:
     - time_step
    */
    template <typename SolverParameters>
    static auto makeTimeGenerator(SolverParameters const& solver_params)
    {
      if constexpr (useful::has_time_step<SolverParameters>::value)
        return ctrw::TimeGenerator_Step{ solver_params.time_step };
      return ctrw::TimeGenerator_Step{ 0. };
    }
    
    /**
      \param velocity_field Velocity field as a function of state.
      \param boundary Boundary condition enforcer.
      \param params_transport Transport parameters.
      \param params_solvers Solver parameters.
      \param dim Spatial dimension.
      \return Advection--diffusion JumpGenerator.
      \note
      \p params_transport must define:
      - diff_coeff
     
      \p params_solvers must define:
      - time_step
     */
    template
    <typename VelocityField,
    typename Boundary,
    typename TransportParameters,
    typename SolverParameters>
    static auto makeJumpGenerator
    (VelocityField&& velocity_field,
     Boundary&& boundary,
     TransportParameters const& params_transport,
     SolverParameters const& params_solvers,
     std::size_t dim)
    {
      if constexpr (useful::has_time_step<SolverParameters>::value)
        return
          ctrw::JumpGenerator_Add{
            ctrw::JumpGenerator_Velocity{
              std::forward<VelocityField>(velocity_field),
              params_solvers.time_step },
            ctrw::JumpGenerator_Diffusion{
              params_transport.diff_coeff,
              params_solvers.time_step,
              dim } };
      return
        ctrw::JumpGenerator_Add{
          ctrw::JumpGenerator_Velocity{
            std::forward<VelocityField>(velocity_field),
            0. },
          ctrw::JumpGenerator_Diffusion{
            params_transport.diff_coeff,
            0.,
            dim } };
    }
    
    /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure advection JumpGenerator.
     \note
     \p params_solvers must define:
     - time_step
    */
    template
    <typename VelocityField,
    typename Boundary,
    typename TransportParameters,
    typename SolverParameters>
    static auto makeJumpGenerator_Advection
    (VelocityField&& velocity_field,
     Boundary&& boundary,
     TransportParameters const& params_transport,
     SolverParameters const& params_solvers,
     std::size_t dim)
    {
      if constexpr (useful::has_time_step<SolverParameters>::value)
        return
          ctrw::JumpGenerator_Velocity{
            std::forward<VelocityField>(velocity_field),
            params_solvers.time_step };
      return
        ctrw::JumpGenerator_Velocity{
          std::forward<VelocityField>(velocity_field),
          0. };
    }
  };
}

#endif /* PTOF_STEPPERS_H */
