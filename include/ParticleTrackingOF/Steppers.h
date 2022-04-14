//
//  Steppers.h
//
//  Created by Tomás Aquino on 09/03/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Steppers_OF_h
#define Steppers_OF_h

#include <cstddef>
#include <utility>
#include "Stochastic/CTRW/JumpGenerator.h"
#include "Stochastic/CTRW/TimeGenerator.h"

namespace ptof
{
  // Time steppers for RK4 advection
  // and stochastic forward Euler diffusion
  struct Steppers_Advection_RK4_Diffusion_Euler
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
    
    // Make advection--diffusion JumpGenerator
    // Given velocity field, boundary enforcer,
    // transport parameters, solver parameters,
    // and spatial dimension
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
  
  // Time steppers for forward Euler advection
  // and stochastic forward Euler diffusion
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
    
    // Make advection--diffusion JumpGenerator
    // Given velocity field, boundary enforcer,
    // transport parameters, solver parameters,
    // and spatial dimension
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


#endif /* Steppers_OF_h */
