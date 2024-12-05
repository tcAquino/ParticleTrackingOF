/**
   \file PTOF/Transitions.h
   \author Tomás Aquino
   \date 09/03/2022
   \brief Objects to handle particle dynamics.
*/

#ifndef PTOF_TRANSITIONS_H
#define PTOF_TRANSITIONS_H

#include "CTRW/Transitions.h"
#include "General/Meta.h"
#include "PTOF/Field.h"
#include "PTOF/TimeStepAdaptor.h"

namespace ptof {
/** \brief Make Transitions object to handle advective-diffusive transport. */
template <typename Steppers, typename Geometry, typename VelocityField,
          typename Boundary, typename TransportParameters,
          typename ReactionParameters, typename SolverParameters>
auto makeTransportTransitions(VelocityField const &velocity_field,
                              Geometry const &geometry, Boundary &boundary,
                              TransportParameters const &params_transport,
                              ReactionParameters const &params_reaction,
                              SolverParameters const &params_solvers) {
  return ctrw::Transitions_AdaptiveTimeStep_Time_Position{
      TimeStepAdaptor_CellSize_SurfaceReaction{
          geometry, velocity_field, boundary.surface_reaction, params_transport,
          params_reaction, params_solvers, meta::Selector_t<CheckOptions::Check>{}},
      Steppers::makeTimeGenerator(params_solvers),
      Steppers::template makeJumpGenerator<typename Geometry::ParallelOption>(
          velocity_field, boundary, params_transport, params_solvers,
          geometry.dim),
      geometry.locator, boundary};
}

/** \brief Make Transitions object to handle purely-advective transport. */
template <typename Steppers, typename Geometry, typename VelocityField,
          typename Boundary, typename TransportParameters,
          typename ReactionParameters, typename SolverParameters>
auto makeTransportTransitions_Advection(
    VelocityField const &velocity_field, Geometry const &geometry,
    Boundary &boundary, TransportParameters const &params_transport,
    ReactionParameters const &params_reaction,
    SolverParameters const &params_solvers) {
  return ctrw::Transitions_AdaptiveTimeStep_Time_Position{
      TimeStepAdaptor_CellSize_SurfaceReaction{
          geometry, velocity_field, boundary.surface_reaction, params_transport,
          params_reaction, params_solvers, meta::Selector_t<CheckOptions::Check>{}},
      Steppers::makeTimeGenerator(params_solvers),
      Steppers::makeJumpGenerator_Advection(velocity_field, boundary,
                                            params_transport, params_solvers,
                                            geometry.dim),
      geometry.locator, boundary};
}

/**
   \brief Make Transitions object to handle transport and bulk reaction.
   \details Conservative transitions and reaction are determined from template
   parameters that must be explicitly specified:
   - \tparam Transport Conservative transport information.
   - \tparam Solvers Numerical solver information.
   - \tparam Reaction Reaction information.
*/
template <typename Transport, typename Solvers, typename Reaction,
          typename VelocityField, typename Geometry, typename Boundary>
auto makeTransitions(VelocityField const &velocity_field,
                     Geometry const &geometry, Boundary &boundary,
                     typename Transport::Parameters const &params_transport,
                     typename Reaction::Parameters const &params_reaction,
                     typename Solvers::Parameters const &params_solvers) {
  return ctrw::Transitions_CTRW_Transport_Reaction{
    Transport::template makeTransitions<Solvers>(
        velocity_field, geometry, boundary, params_transport,
        params_reaction, params_solvers),
    Reaction::makeBulkReaction(geometry, params_reaction, params_transport,
                               params_solvers)};
}

/**
   \brief Make Transitions object to handle transport and bulk reaction.
   \details Conservative transitions are determined from template parameters
   that must be explicitly specified:
   - \tparam Transport Conservative transport information.
   - \tparam Solvers Numerical solver information.
*/
template <typename Transport, typename Solvers, typename VelocityField,
          typename Geometry, typename Boundary, typename ReactionParameters,
          typename BulkReaction>
auto makeTransitions(VelocityField const &velocity_field,
                     Geometry const &geometry, Boundary &boundary,
                     typename Transport::Parameters const &params_transport,
                     ReactionParameters const &params_reaction,
                     typename Solvers::Parameters const &params_solvers,
                     BulkReaction const &bulk_reaction) {
  return ctrw::Transitions_CTRW_Transport_Reaction{
      Transport::template makeTransitions<Solvers>(
          velocity_field, geometry, boundary, params_transport, params_reaction,
          params_solvers),
      bulk_reaction};
}
} // namespace ptof

#endif /* PTOF_TRANSITIONS_H */
