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
#include "PTOF/Meta.h"
#include "PTOF/Reaction.h"
#include "PTOF/TimeGenerator.h"

namespace ptof {
template <typename Transitions_Transport, typename BulkReaction,
          typename SurfaceReaction>
class Transitions_Transport_Reaction {
public:
  Transitions_Transport_Reaction(Transitions_Transport &&transitions_transport,
                                 BulkReaction &&bulk_reaction,
                                 SurfaceReaction &&surface_reaction)
      : _transitions_transport{std::forward<Transitions_Transport>(
            transitions_transport)},
        _bulk_reaction{std::forward<BulkReaction>(bulk_reaction)},
        _surface_reaction{std::forward<SurfaceReaction>(surface_reaction)} {}

  template <typename State> void operator()(State &state) {
    if constexpr (meta::has_adsorbed_v<typename State::Info>) {
      _surface_reaction.desorb(state);
      return;
    }
    auto time_old = state.time;
    _transitions_transport(state);
    _bulk_reaction(state, state.time - time_old);
  }

  auto const &transitions_transport() const { return _transitions_transport; }

  auto const &bulk_reaction() const { return _bulk_reaction; }

  auto const &surface_reaction() const { return _surface_reaction; }

private:
  Transitions_Transport _transitions_transport;
  BulkReaction _bulk_reaction;
  SurfaceReaction _surface_reaction;
};
template <typename Transitions_Transport, typename BulkReaction,
          typename SurfaceReaction>
Transitions_Transport_Reaction(Transitions_Transport &&, BulkReaction &&,
                               SurfaceReaction &&)
    -> Transitions_Transport_Reaction<Transitions_Transport, BulkReaction,
                                      SurfaceReaction>;

/** \brief Make Transitions object to handle transport. */
template <typename Solvers, typename Geometry, typename VelocityField,
          typename Boundary, typename TransportParameters,
          typename ReactionParameters>
auto makeTransportTransitions(
    VelocityField const &velocity_field, Geometry const &geometry,
    Boundary &boundary, TransportParameters const &params_transport,
    ReactionParameters const &params_reaction,
    typename Solvers::Parameters const &params_solvers) {
  return ctrw::Transitions_Time_Position{
      Solvers::makeTimeGenerator(geometry, velocity_field, boundary,
                                 params_transport, params_reaction,
                                 params_solvers),
      Solvers::template makeJumpGenerator<typename Geometry::ParallelOption>(
          velocity_field, boundary, params_transport, params_solvers,
          geometry.dim),
      boundary};
}

/**
   \brief Make Transitions object to handle transport and reaction.
   \details Conservative transitions and reaction are determined from template
   parameters that must be explicitly specified:
   - \tparam Transport Conservative transport information.
   - \tparam ReactionHandler Reaction information.
*/
template <typename Transport, typename ReactionHandler, typename VelocityField,
          typename Geometry, typename Boundary>
auto makeTransitions(
    VelocityField const &velocity_field, Geometry const &geometry,
    Boundary &boundary, typename Transport::Parameters const &params_transport,
    typename ReactionHandler::Parameters const &params_reaction,
    typename Transport::Solvers::Parameters const &params_solvers) {
  return Transitions_Transport_Reaction{
      Transport::makeTransitions(velocity_field, geometry, boundary,
                                 params_transport, params_reaction,
                                 params_solvers),
      boundary.surface_reaction,
      ReactionHandler::makeBulkReaction(geometry, params_reaction,
                                        params_transport, params_solvers)};
}

/**
   \brief Make Transitions object to handle transport and reaction.
   \details Conservative transitions are determined from template parameters
   that must be explicitly specified:
   - \tparam Transport Conservative transport information.
*/
template <typename Transport, typename VelocityField, typename Geometry,
          typename Boundary, typename ReactionParameters, typename BulkReaction>
auto makeTransitions(
    VelocityField const &velocity_field, Geometry const &geometry,
    Boundary &boundary, typename Transport::Parameters const &params_transport,
    ReactionParameters const &params_reaction,
    typename Transport::Solvers::Parameters const &params_solvers,
    BulkReaction &bulk_reaction) {
  return Transitions_Transport_Reaction{
      Transport::makeTransitions(velocity_field, geometry, boundary,
                                 params_transport, params_reaction,
                                 params_solvers),
      bulk_reaction, boundary.surface_reaction};
}
} // namespace ptof

#endif /* PTOF_TRANSITIONS_H */
