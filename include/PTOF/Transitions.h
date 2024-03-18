/**
 \file PTOF/Transitions.h
 \author Tomás Aquino
 \date 09/03/2022
*/

#ifndef PTOF_TRANSITIONS_H
#define PTOF_TRANSITIONS_H

#include "CTRW/Transitions.h"
#include "PTOF/Field.h"
#include "PTOF/TimeStepAdaptor.h"

namespace ptof
{
  /** \brief Make transitions object to handle advective-diffusive transport. */
  template
  <typename Steppers,
  typename Geometry,
  typename VelocityField,
  typename Boundary,
  typename TransportParameters,
  typename ReactionParameters,
  typename SolverParameters>
  auto makeTransportTransitions
  (VelocityField const& velocity_field,
   Geometry const& geometry,
   Boundary& boundary,
   TransportParameters const& params_transport,
   ReactionParameters const& params_reaction,
   SolverParameters const& params_solvers)
  {
    return
      ctrw::Transitions_AdaptiveTimeStep_Time_Position{
        TimeStepAdaptor_CellSize_SurfaceReaction{
          geometry,
          velocity_field,
          boundary.surface_reaction,
          params_transport,
          params_reaction,
          params_solvers,
          CheckOptions::Check{} },
        Steppers::makeTimeGenerator(params_solvers),
        Steppers::makeJumpGenerator(velocity_field,
                                    boundary,
                                    params_transport,
                                    params_solvers,
                                    geometry.dim),
        boundary };
  }
  
  /** \brief Make transitions object to handle purely-advective transport. */
  template
  <typename Steppers,
  typename Geometry,
  typename VelocityField,
  typename Boundary,
  typename TransportParameters,
  typename ReactionParameters,
  typename SolverParameters>
  auto makeTransportTransitions_Advection
  (VelocityField const& velocity_field,
   Geometry const& geometry,
   Boundary& boundary,
   TransportParameters const& params_transport,
   ReactionParameters const& params_reaction,
   SolverParameters const& params_solvers)
  {
    return
      ctrw::Transitions_AdaptiveTimeStep_Time_Position{
        TimeStepAdaptor_CellSize_SurfaceReaction{
          geometry,
          velocity_field,
          boundary.surface_reaction,
          params_transport,
          params_reaction,
          params_solvers,
          CheckOptions::Check{} },
        Steppers::makeTimeGenerator(params_solvers),
        Steppers::makeJumpGenerator_Advection(velocity_field,
                                              boundary,
                                              params_transport,
                                              params_solvers,
                                              geometry.dim),
        boundary };
  }
}

#endif /* PTOF_TRANSITIONS_H */
