/**
* \file PTOF/Transitions.h
* \author Tomás Aquino
* \date 09/03/2022
*/

#ifndef PTOF_TRANSITIONS_H
#define PTOF_TRANSITIONS_H

#include "CTRW/Transitions.h"
#include "PTOF/TimeStepAdaptor.h"

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
   Boundary& boundary,
   TransportParameters const& params_transport,
   SolverParameters const& params_solvers)
  {
    return
      ctrw::Transitions_AdaptiveTimeStep_Time_Position{
        TimeStepAdaptor_CellSize_SurfaceReaction{
          geometry,
          velocity_field,
          boundary.reaction,
          params_transport,
          params_solvers
        },
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
   Boundary& boundary,
   TransportParameters const& params_transport,
   SolverParameters const& params_solvers)
  {
    return
      ctrw::Transitions_AdaptiveTimeStep_Time_Position{
        TimeStepAdaptor_CellSize_SurfaceReaction{
          geometry,
          velocity_field,
          boundary.reaction,
          params_transport,
          params_solvers },
        Steppers::makeTimeGenerator(params_solvers),
        Steppers::makeJumpGenerator_Advection(std::forward<VelocityField>(velocity_field),
                                              boundary,
                                              params_transport,
                                              params_solvers,
                                              geometry.dim),
        boundary };
  }
}

#endif /* PTOF_TRANSITIONS_H */
