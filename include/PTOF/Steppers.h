/**
* \file PTOF/Steppers.h
* \author Tomás Aquino
* \date 09/03/2022
*/

#ifndef PTOF_STEPPERS_H
#define PTOF_STEPPERS_H

#include <cstddef>
#include <utility>
#include "CTRW/JumpGenerator.h"
#include "CTRW/TimeGenerator.h"

namespace ptof
{
  /** \struct Steppers_Advection_RK4_Diffusion_Euler PTOF/Steppers.h "PTOF/Steppers.h"
   * \brief  Time steppers for RK4 advection
   * and stochastic forward Euler diffusion. */
  struct Steppers_Advection_RK4_Diffusion_Euler
  {
    /** Make constant time step TimeGenerator. */
    template <typename Parameters>
    static auto makeTimeGenerator(Parameters const& params)
    {
      return
        ctrw::TimeGenerator_Step{
          params.time_step
        };
    }
    
    /** Make advection--diffusion JumpGenerator.
     * Given:
     * - velocity field;
     * - boundary enforcer;
     * - transport parameters; 
     * - solver parameters;
     * - spatial dimension. */
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
      return
        ctrw::JumpGenerator_Add{
          ctrw::JumpGenerator_Velocity_State_RK4{
            std::forward<VelocityField>(velocity_field),
            params_solvers.time_step,
            std::forward<Boundary>(boundary) },
          ctrw::JumpGenerator_Diffusion{
            params_transport.diff_coeff,
            params_solvers.time_step,
            dim } };
    }
    
    /** Make purely-advective JumpGenerator.
     * Given:
     * - velocity field; 
     * - boundary enforcer;
     * - transport parameters;
     * - solver parameters;
     * - spatial dimension. */
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
      return
        ctrw::JumpGenerator_Velocity_State_RK4{
          std::forward<VelocityField>(velocity_field),
          params_solvers.time_step,
          std::forward<Boundary>(boundary) };
    }
  };
  
  /** \struct Steppers_Advection_Euler_Diffusion_Euler PTOF/Steppers.h "PTOF/Steppers.h"
   * \brief  Time steppers for forward Euler advection
   * and stochastic forward Euler diffusion. */
  struct Steppers_Advection_Euler_Diffusion_Euler
  {
    // Make constant time step TimeGenerator
    template <typename Parameters>
    static auto makeTimeGenerator(Parameters const& params)
    {
      return
        ctrw::TimeGenerator_Step{
          params.time_step
        };
    }
    
    /** Make advection--diffusion JumpGenerator.
     * Given:
     * - velocity field;
     * - boundary enforcer;
     * - transport parameters; 
     * - solver parameters;
     * - spatial dimension. */
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
      return
        ctrw::JumpGenerator_Add{
          ctrw::JumpGenerator_Velocity_State{
            std::forward<VelocityField>(velocity_field),
            params_solvers.time_step,
            std::forward<Boundary>(boundary)
          },
          ctrw::JumpGenerator_Diffusion{
            params_transport.diff_coeff,
            params_solvers.time_step,
            dim
          }
        };
    }
  };
}


#endif /* PTOF_STEPPERS_H */
