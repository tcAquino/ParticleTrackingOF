/**
   \file PTOF/TimeGenerator.h
   \author Tomas Aquino
   \date 09/26/2017
   \brief Objects time increment generation.
*/

#ifndef PTOF_TIMEGENERATOR_H
#define PTOF_TIMEGENERATOR_H

#include "CTRW/Meta.h"
#include "General/Constant.h"
#include "General/Meta.h"
#include "PTOF/CheckOptions.h"
#include "PTOF/SurfaceReaction.h"
#include "PTOF/Useful.h"
#include <MinMax.H>
#include <algorithm>
#include <fieldTypes.H>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <zero.H>

namespace ptof {
/**
   \class TimeGenerator_Adaptive PTOF/TimeGenerator.h
   "PTOF/TimeGenerator.h"
   \brief Adaptive time step control based on local cell size.
   \details Takes reaction limitations to time step only if close to reactive
   boundary.
   \note See CheckOptions class for bounds checking options.
*/
template <typename Geometry, typename VelocityField, typename SurfaceReaction,
          typename CheckOption = CheckOptions::Check>
class TimeGenerator_Adaptive {
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
     \param velocity_field Velocity field as a function of state.
     \param surface_reaction Surface reaction handler.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
  */
  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  TimeGenerator_Adaptive(Geometry const &geometry,
                         VelocityField const &velocity_field,
                         SurfaceReaction &surface_reaction,
                         TransportParameters const &params_transport,
                         ReactionParameters const &params_reaction,
                         SolverParameters const &params_solvers)
      : _geometry{geometry}, _velocity_field{velocity_field},
        _surface_reaction{surface_reaction}, _volume_factor{1.},
        _diff_coeff{params_transport.diff_coeff},
        _min_global_local{params_solvers.min_global_local},
        _local_time_step_adv{params_solvers.local_time_step_adv},
        _local_time_step_diff{params_solvers.local_time_step_diff},
        _local_time_step_surf_react{params_solvers.local_time_step_surf_react},
        _global_time_step{compute_global_time_step(
            params_transport, params_reaction, params_solvers)} {
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

  /**
   \brief Constructor.
   \param geometry Domain geometry info and utilities.
   \param velocity_field Velocity field as a function of state.
   \param surface_reaction Surface reaction handler.
   \param params_transport Transport parameters.
   \param params_reaction Reaction parameters.
   \param params_solvers Solver parameters.
*/
  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  TimeGenerator_Adaptive(Geometry const &geometry,
                         VelocityField const &velocity_field,
                         SurfaceReaction &surface_reaction,
                         TransportParameters const &params_transport,
                         ReactionParameters const &params_reaction,
                         SolverParameters const &params_solvers, CheckOption)
      : TimeGenerator_Adaptive{geometry,         velocity_field,
                               surface_reaction, params_transport,
                               params_reaction,  params_solvers} {}

  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  TimeGenerator_Adaptive(Geometry &&geometry,
                         VelocityField const &velocity_field,
                         SurfaceReaction &surface_reaction,
                         TransportParameters const &params_transport,
                         ReactionParameters const &params_reaction,
                         SolverParameters const &params_solvers) = delete;

  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  TimeGenerator_Adaptive(Geometry const &geometry,
                         VelocityField &&velocity_field,
                         SurfaceReaction &surface_reaction,
                         TransportParameters const &params_transport,
                         ReactionParameters const &params_reaction,
                         SolverParameters const &params_solvers) = delete;

  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  TimeGenerator_Adaptive(Geometry &&geometry, VelocityField &&velocity_field,
                         SurfaceReaction &surface_reaction,
                         TransportParameters const &params_transport,
                         ReactionParameters const &params_reaction,
                         SolverParameters const &params_solvers) = delete;

  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  TimeGenerator_Adaptive(Geometry &&geometry,
                         VelocityField const &velocity_field,
                         SurfaceReaction &surface_reaction,
                         TransportParameters const &params_transport,
                         ReactionParameters const &params_reaction,
                         SolverParameters const &params_solvers,
                         CheckOption) = delete;

  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  TimeGenerator_Adaptive(Geometry const &geometry,
                         VelocityField &&velocity_field,
                         SurfaceReaction &surface_reaction,
                         TransportParameters const &params_transport,
                         ReactionParameters const &params_reaction,
                         SolverParameters const &params_solvers,
                         CheckOption) = delete;

  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  TimeGenerator_Adaptive(Geometry &&geometry, VelocityField &&velocity_field,
                         SurfaceReaction &surface_reaction,
                         TransportParameters const &params_transport,
                         ReactionParameters const &params_reaction,
                         SolverParameters const &params_solvers,
                         CheckOption) = delete;

  /**
     \brief Compute adaptive time step.
     \param state Particle state.
     \return New time step.
  */
  template <typename State> double operator()(State const &state) const {
    auto cell_id = state.cell;
    if constexpr (check_if_outside) {
      if (outside<warn_if_outside>(cell_id, make_point(state.position),
                                   "Using nearest cell.")) {
        cell_id = _geometry.locator.nearest_cell(state);
      }
    }
    double cell_side =
        (_constrain_local_adv || _constrain_local_diff)
            ? std::pow(_geometry.mesh().V()[cell_id] / _volume_factor,
                       1. / Geometry::dim)
            : 0.;

    double local_advection_time = advection_time(state, cell_side);
    double local_diffusion_time = diffusion_time(state, cell_side);
    double local_time_step =
        std::min({_local_time_step_adv * local_advection_time,
                  _local_time_step_diff * local_diffusion_time});
    if constexpr (!std::is_same_v<SurfaceReaction, SurfaceReaction_DoNothing>) {
      if (_constrain_local_surf_react) {
        auto const &faces = _geometry.mesh().cells()[cell_id];
        double max_face_rate = 0.;
        for (auto face : faces) {
          double rate = _surface_reaction.rate(face);
          if (rate > max_face_rate) {
            max_face_rate = rate;
          }
        }
        double local_reaction_time =
            _diff_coeff / (cnst::pi * max_face_rate * max_face_rate);
        local_time_step = std::min(
            local_time_step, _local_time_step_surf_react * local_reaction_time);
      }
    }

    double time_step = _min_global_local
                           ? std::min(local_time_step, _global_time_step)
                           : std::max(local_time_step, _global_time_step);
    if constexpr (meta::has_time_step_setter_v<SurfaceReaction>) {
      _surface_reaction.time_step(time_step);
    }

    return time_step;
  }

private:
  /**
     \brief Compute global time step constraint.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
     \note Deals with edge cases where constraint values are zero or infinity.
  */
  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  double compute_global_time_step(TransportParameters const &params_transport,
                                  ReactionParameters const &params_reaction,
                                  SolverParameters const &params_solvers) {
    bool global_constraints_are_zero =
        (params_solvers.global_time_step_adv == 0. ||
         params_solvers.global_time_step_diff == 0. ||
         params_solvers.global_time_step_surf_react == 0.);
    // bool global_transport_constraints_are_inf =
    //     (params_solvers.global_time_step_adv ==
    //          std::numeric_limits<double>::infinity() &&
    //      params_solvers.global_time_step_diff &&
    //      std::numeric_limits<double>::infinity());

    bool constrain_global_adv = (!global_constraints_are_zero &&
                                 params_solvers.global_time_step_adv !=
                                     std::numeric_limits<double>::infinity());
    double global_time_step_adv = constrain_global_adv
                                      ? params_solvers.global_time_step_adv *
                                            params_transport.advection_time
                                      : params_solvers.global_time_step_adv;

    bool constrain_global_diff = (!global_constraints_are_zero &&
                                  params_solvers.global_time_step_diff !=
                                      std::numeric_limits<double>::infinity());
    double global_time_step_diff = constrain_global_diff
                                       ? params_solvers.global_time_step_diff *
                                             params_transport.diffusion_time
                                       : params_solvers.global_time_step_diff;

    bool constrain_global_react = (!global_constraints_are_zero &&
                                   params_solvers.global_time_step_surf_react !=
                                       std::numeric_limits<double>::infinity());
    double global_time_step_surf_react =
        constrain_global_react ? params_solvers.global_time_step_surf_react *
                                     params_reaction.reaction_time
                               : params_solvers.global_time_step_surf_react;

    return std::min({global_time_step_adv, global_time_step_diff,
                     global_time_step_surf_react});
  }

  /** \brief Compute local advection time. */
  template <typename State>
  double advection_time(State const &state, double cell_side) const {
    if constexpr (!std::is_same_v<VelocityField, meta::Empty>) {
      return _constrain_local_adv
                 ? cell_side / Foam::mag(_velocity_field(state))
                 : 1.;
    }
    return std::numeric_limits<double>::infinity();
  }

  /** \brief Compute local advection time. */
  template <typename State>
  double diffusion_time(State const &state, double cell_side) const {
    return _constrain_local_diff ? cell_side * cell_side / (2. * _diff_coeff)
                                 : 1.;
  }

  Geometry const &_geometry; /**< Domain geometry info and utilities. */
  VelocityField const &_velocity_field; /**< Velocity field. */
  SurfaceReaction &_surface_reaction;   /**< Surface reaction. */
  double _volume_factor;  /**< Volume factor to correct cell volumes in 1D and
                             2D. */
  double _diff_coeff;     /**< Diffusion coefficient. */
  bool _min_global_local; /**< Use minimum of local and global constraints if
                             true, otherwise use maximum. */
  double _local_time_step_adv;  /**< Local advection time step constraint. */
  double _local_time_step_diff; /**< Local diffusion time step constraint. */
  double _local_time_step_surf_react; /**< Local surface reaction time step
                                         constraint. */
  double _global_time_step;           /**< Global time step constraint. */

  bool _local_constraints_are_zero =
      _local_time_step_adv == 0. || _local_time_step_diff == 0. ||
      _local_time_step_surf_react ==
          0.; /**< Whether local constraints are all zero. */
  bool _local_transport_constraints_are_inf =
      _local_time_step_adv == std::numeric_limits<double>::infinity() &&
      _local_time_step_diff &&
      std::numeric_limits<double>::infinity(); /**< Whether local transport
                                                   constraints are all
                                                   infinity. */
  bool _constrain_local_adv =
      !_local_constraints_are_zero &&
      _local_time_step_diff !=
          std::numeric_limits<double>::infinity(); /**< Whether to constrain
                                                      local advection. */
  bool _constrain_local_diff =
      !_local_constraints_are_zero &&
      _local_time_step_diff !=
          std::numeric_limits<double>::infinity(); /**< Whether to constrain
                                 local diffusion. */
  bool _constrain_local_surf_react =
      !_local_constraints_are_zero &&
      _local_time_step_diff !=
          std::numeric_limits<double>::infinity(); /**< Whether to constrain
                                 local surface reaction. */
};
template <typename Geometry, typename VelocityField, typename SurfaceReaction,
          typename TransportParameters, typename ReactionParameters,
          typename SolverParameters>
TimeGenerator_Adaptive(Geometry const &, VelocityField const &,
                       SurfaceReaction &, TransportParameters const &,
                       ReactionParameters const &, SolverParameters const &)
    -> TimeGenerator_Adaptive<Geometry, VelocityField, SurfaceReaction,
                              CheckOptions::Check>;
template <typename Geometry, typename VelocityField, typename SurfaceReaction,
          typename TransportParameters, typename ReactionParameters,
          typename SolverParameters, typename CheckOption>
TimeGenerator_Adaptive(Geometry const &, VelocityField const &,
                       SurfaceReaction &, TransportParameters const &,
                       ReactionParameters const &, SolverParameters const &,
                       CheckOption)
    -> TimeGenerator_Adaptive<Geometry, VelocityField, SurfaceReaction,
                              CheckOption>;
} // namespace ptof

#endif /* PTOF_TIMEGENERATOR_H */
