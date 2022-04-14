//
//  Transitions.h
//
//  Created by Tomás Aquino on 09/03/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Transitions_h
#define Transitions_h

#include "Stochastic/CTRW/Transitions_State.h"

namespace ptof
{
  // Make transitions object
  // to handle conservative transport
  // given velocity field,
  // geometry information,
  // boundary condition enforcer,
  // transport parameters,
  // and solver parameters
  template
  <typename Geometry,
  typename Steppers,
  typename VelocityField,
  typename Boundary,
  typename TransportParameters,
  typename SolverParameters>
  auto makeTransportTransitions
  (VelocityField&& velocity_field,
   Geometry const& geometry,
   Boundary const& boundary,
   TransportParameters const& params_transport,
   SolverParameters const& params_solvers)
  {
    return
    ctrw::Transitions_Locate_Time_Position{
      velocity_field.get_locator(),
      Steppers::makeTimeGenerator(params_solvers),
      Steppers::makeJumpGenerator(std::forward<VelocityField>(velocity_field),
                        boundary,
                        params_transport,
                        params_solvers,
                        geometry.dim),
      boundary
    };
  }
}


#endif /* Transitions_h */
