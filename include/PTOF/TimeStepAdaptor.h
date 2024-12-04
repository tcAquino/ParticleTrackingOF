/**
   \file PTOF/TimeStepAdaptor.h
   \author Tomás Aquino
   \date 09/26/2017
   \brief Objects for adaptive time stepping.
*/

#ifndef PTOF_TIMESTEPADAPTOR_H
#define PTOF_TIMESTEPADAPTOR_H

#include "CTRW/Meta.h"
#include "General/Meta.h"
#include "PTOF/Field.h"
#include "PTOF/Reaction.h"
#include "PTOF/Useful.h"
#include <MinMax.H>
#include <algorithm>
#include <fieldTypes.H>
#include <zero.H>

namespace ptof {

/**
   \class TimeStepAdaptor_CellSize PTOF/TimeStepAdaptor.h
   "PTOF/TimeStepAdaptor.h"
   \brief Adaptive time step control based on local cell size.
   \details Takes reaction limitations to time step only if close to reactive
   boundary.
   \note See CheckOptions class for bounds checking options.
*/
template <typename Geometry, typename VelocityField, typename SurfaceReaction,
          typename TransportParameters, typename ReactionParameters,
          typename SolverParameters, typename CheckOption = CheckOptions::Check>
class TimeStepAdaptor_CellSize_SurfaceReaction {
public:
  /** Whether to check if requested positions are outside of mesh. */
  static constexpr bool check_if_outside =
      !std::is_same_v<CheckOption, CheckOptions::NoCheck>;

  /** Whether to warn when requested positions are outside of mesh. */
  static constexpr bool warn_if_outside =
      std::is_same_v<CheckOption, CheckOptions::Warn>;

  /**
     \brief Constructor.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field.
     \param surface_reaction Surface reaction handler.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
  */
  TimeStepAdaptor_CellSize_SurfaceReaction(
      Geometry const &geometry, VelocityField const &velocity_field,
      SurfaceReaction &surface_reaction,
      TransportParameters const &params_transport,
      ReactionParameters const &params_reaction,
      SolverParameters const &params_solvers, meta::Selector_t<CheckOption>)
      : TimeStepAdaptor_CellSize_SurfaceReaction{
            geometry,         velocity_field,  surface_reaction,
            params_transport, params_reaction, params_solvers} {}

  /**
     \brief Constructor.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field.
     \param surface_reaction Surface reaction handler.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
  */
  TimeStepAdaptor_CellSize_SurfaceReaction(
      Geometry const &geometry, VelocityField const &velocity_field,
      SurfaceReaction &surface_reaction,
      TransportParameters const &params_transport,
      ReactionParameters const &params_reaction,
      SolverParameters const &params_solvers)
      : _geometry{geometry}, _velocity_field{velocity_field},
        _surface_reaction{surface_reaction}, _volume_factor{1.},
        _diff_coeff{params_transport.diff_coeff},
        _local_time_step_adv{params_solvers.local_time_step_adv},
        _local_time_step_diff{params_solvers.local_time_step_diff},
        _local_time_step_react{params_solvers.local_time_step_react},
        _time_step_global{std::min({params_solvers.global_time_step_adv *
                                        params_transport.advection_time,
                                    params_solvers.global_time_step_diff *
                                        params_transport.diffusion_time,
                                    params_solvers.global_time_step_react *
                                        params_reaction.reaction_time})} {
    for (std::size_t dd = 2; dd-- > Geometry::dim;) {
      Foam::vector direction = Foam::zero{};
      direction[dd] = 1.;
      std::vector<double> max;
      std::vector<double> min;
      for (auto const &patch : _geometry.mesh().boundary()) {
        max.push_back(Foam::max(patch.Cf() & direction));
        min.push_back(Foam::min(patch.Cf() & direction));
      }
      _volume_factor *= *std::max_element(max.begin(), max.end()) -
                        *std::min_element(min.begin(), min.end());
    }
  }

  TimeStepAdaptor_CellSize_SurfaceReaction(
      Geometry &&geometry, VelocityField const &velocity_field,
      SurfaceReaction &surface_reaction,
      TransportParameters const &params_transport,
      ReactionParameters const &params_reaction,
      SolverParameters const &params_solvers) = delete;

  TimeStepAdaptor_CellSize_SurfaceReaction(
      Geometry const &geometry, VelocityField &&velocity_field,
      SurfaceReaction &surface_reaction,
      TransportParameters const &params_transport,
      ReactionParameters const &params_reaction,
      SolverParameters const &params_solvers) = delete;

  TimeStepAdaptor_CellSize_SurfaceReaction(
      Geometry &&geometry, VelocityField &&velocity_field,
      SurfaceReaction &surface_reaction,
      TransportParameters const &params_transport,
      ReactionParameters const &params_reaction,
      SolverParameters const &params_solvers) = delete;

  /**
     \brief Update time step.
     \param state Particle state.
     \param time_generator Time step generator object to be updated.
     \param jump_generator Jump generator object to be updated.
     \return New time step.
  */
  template <typename State, typename TimeGenerator, typename JumpGenerator>
  double operator()(State &state, TimeGenerator &time_generator,
                    JumpGenerator &jump_generator) const {
    auto cell_id = state.cell;
    if constexpr (check_if_outside)
      if (outside<warn_if_outside>(state.cell, make_point(state.position),
                                   "Using nearest cell."))
        cell_id = _geometry.locator.nearest_cell(state);
    double cell_side = std::pow(_geometry.mesh().V()[cell_id] / _volume_factor,
                                1. / Geometry::dim);
    double local_advection_time = cell_side / Foam::mag(_velocity_field(state));
    double local_diffusion_time = cell_side * cell_side / (2. * _diff_coeff);

    double reaction_rate = 0.;
    if constexpr (!std::is_same_v<SurfaceReaction, SurfaceReaction_DoNothing>)
      for (auto face : _geometry.mesh().cells()[cell_id]) {
        double rate = _surface_reaction.rate(face);
        if (rate > reaction_rate)
          reaction_rate = rate;
      }
    double local_reaction_time = _diff_coeff / (reaction_rate * reaction_rate);

    double time_step =
        std::max(std::min({_local_time_step_adv * local_advection_time,
                           _local_time_step_diff * local_diffusion_time,
                           _local_time_step_react * local_reaction_time}),
                 _time_step_global);

    if constexpr (meta::has_time_step_setter_v<TimeGenerator>)
      time_generator.time_step(time_step);
    if constexpr (meta::has_time_step_setter_v<JumpGenerator>)
      jump_generator.time_step(time_step);
    if constexpr (meta::has_time_step_setter_v<SurfaceReaction>)
      _surface_reaction.time_step(time_step);

    return time_step;
  }

private:
  Geometry const &_geometry; /**< Domain geometry info and utilities. */
  VelocityField const &_velocity_field; /**< Velocity field. */
  SurfaceReaction &_surface_reaction;   /**< Surface reaction. */
  double _volume_factor; /**< Volume factor to correct cell volumes in 1D and
                            2D. */
  double _diff_coeff;
  double _local_time_step_adv;
  double _local_time_step_diff;
  double _local_time_step_react;
  double _time_step_global;
};
template <typename Geometry, typename VelocityField, typename SurfaceReaction,
          typename TransportParameters, typename ReactionParameters,
          typename SolverParameters, typename CheckOption>
TimeStepAdaptor_CellSize_SurfaceReaction(
    Geometry const &, VelocityField const &, SurfaceReaction &,
    TransportParameters const &, ReactionParameters const &,
    SolverParameters const &, meta::Selector_t<CheckOption>)
    -> TimeStepAdaptor_CellSize_SurfaceReaction<
        Geometry, VelocityField, SurfaceReaction, TransportParameters,
        ReactionParameters, SolverParameters, CheckOption>;
} // namespace ptof

#endif /* PTOF_TIMESTEPADAPTOR_H */
