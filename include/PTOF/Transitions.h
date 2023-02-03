/**
* \file PTOF/Transitions.h
* \author Tomás Aquino
* \date 09/03/2022
*/

#ifndef PTOF_TRANSITIONS_H
#define PTOF_TRANSITIONS_H

#include "CTRW/Transitions.h"

namespace ptof
{
  /** Make transitions object
   * to handle conservative transport.
   * Given:
   * - velocity field;
   * - geometry information;
   * - boundary condition enforcer;
   * - transport parameters;
   * - solver parameters. */
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
        boundary };
  }
  
  /** Make transitions object
   * to handle purely-advective transport
   * Given:
   * - velocity field;
   * - geometry information;
   * - boundary condition enforcer;
   * - transport parameters;
   * - solver parameters. */
  template
  <typename Geometry,
  typename Steppers,
  typename VelocityField,
  typename Boundary,
  typename TransportParameters,
  typename SolverParameters>
  auto makeTransportTransitions_Advection
  (VelocityField&& velocity_field,
   Geometry const& geometry,
   Boundary const& boundary,
   TransportParameters const& params_transport,
   SolverParameters const& params_solvers)
  {
    return
      ctrw::Transitions_Locate_Time_Position{
        velocity_field.get_locator(),
        [&velocity_field, &params_solvers](auto state){
          return params_solvers.step_length/Foam::mag(velocity_field(state)); },
        Steppers::makeJumpGenerator_Advection(std::forward<VelocityField>(velocity_field),
                                              boundary,
                                              params_transport,
                                              params_solvers,
                                              geometry.dim),
        boundary };
  }
}

#endif /* PTOF_TRANSITIONS_H */
