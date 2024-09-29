/**
 \file PTOF/InitialCondition_Cases.h
 \author Tomás Aquino
 \date 24/02/2022
*/

#ifndef PTOF_INITIALCONDITION_CASES_H
#define PTOF_INITIALCONDITION_CASES_H

#include "CTRW/Meta.h"
#include "General/Useful.h"
#include "PTOF/Boundary.h"
#include "PTOF/InitialCondition.h"
#include "PTOF/InitialConditionList.h"
#include "PTOF/Useful.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fieldTypes.H>
#include <map>
#include <point.H>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace ptof {
/** \class InitialConditionParameters_Cases PTOF/InitialCondition_Cases.h
 * "PTOF/InitialCondition_Cases.h" \brief Initial condition condition parameters
 * to handle all initial condition types in \c InitialCondition_Cases. */
struct InitialConditionParameters_Cases {
public:
  InitialConditionList::Type type;
  double distance_wall;
  std::vector<std::pair<double, double>> region_boundaries;
  std::string filename_data;
  double initial_mass;
  double time_min;
  double time_step_accuracy_adv;
  double time_step_accuracy_diff;
  double time_step_accuracy_react;
  double time_max;
  double time_step;

  template <typename Geometry, typename TransportParameters,
            typename ReactionParameters>
  InitialConditionParameters_Cases(Directories const &directories,
                                   std::string const &name,
                                   Geometry const &geometry,
                                   TransportParameters const &params_transport,
                                   ReactionParameters const &params_reaction) {
    std::string ic_name;
    auto input =
        useful::open_read(directories.dir_parameters +
                          "/parameters_initial_condition_" + name + ".dat");
    useful::read_first_from_line(ic_name, input);
    if (!InitialConditionList{}.contains(ic_name))
      throw std::runtime_error{"Initial condition type " + ic_name +
                               " not supported"};
    type = InitialConditionList::type(ic_name);
    if (type == InitialConditionList::Type::uniform_near_solid)
      useful::read_first_from_line(distance_wall, input);
    if (type == InitialConditionList::Type::uniform_region_cartesian ||
        type == InitialConditionList::Type::fluxweighted_region_cartesian) {
      auto split_line = useful::split_line(input);
      if (split_line.size() % 2 != 0)
        std::runtime_error{
            "Initial condition region boundaries must come in pairs"};
      for (std::size_t ii = 0; ii < split_line.size(); ii += 2) {
        region_boundaries.push_back(
            {std::stod(split_line[ii]), std::stod(split_line[ii + 1])});
      }
    }
    if (type == InitialConditionList::Type::prescribed_positions) {
      useful::read_first_from_line(filename_data, input);
      useful::expand_env_in_place(
          useful::expand_home_dir_in_place(filename_data));
    }
    if (type != InitialConditionList::Type::prescribed_positions_masses &&
        type == InitialConditionList::Type::prescribed_positions_masses_tags)
      useful::read_first_from_line(initial_mass, input);
    useful::read_first_from_line(time_min, input);
    if (type == InitialConditionList::Type::uniform_inlet_continuous ||
        type == InitialConditionList::Type::fluxweighted_inlet_continuous) {
      useful::read_first_from_line(time_max, input);
      useful::read_first_from_line(time_step_accuracy_adv, input);
      useful::read_first_from_line(time_step_accuracy_diff, input);
      useful::read_first_from_line(time_step_accuracy_react, input);
      time_step =
          std::min({time_step_accuracy_adv * params_transport.advection_time,
                    time_step_accuracy_diff * params_transport.diffusion_time,
                    time_step_accuracy_react * params_reaction.reaction_time});
    }
  }

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Initial condition parameters\n"
           "--------------------------------------------------------------\n"
           "- Initial condition type:\n"
           "\tuniform: Homogeneous throughout the domain\n"
           "\tfluxweighted: Flux-weighted throughout the domain\n"
           "\tuniform_inlet: Homogeneous at the inlet\n"
           "\tfluxweighted_inlet: Flux-weighted at the inlet\n"
           "\tuniform_solid: Homogeneous at the solid surface\n"
           "\tuniform_near_solid: Homogeneous at a fixed distance to the\n"
           "\t                    solid interface\n"
           "\tuniform_region_cartesian: Homogeneous in a prescribed\n"
           "\t                          cartesian region\n"
           "\tfluxweighted_region_cartesian: Flux-weighted in a prescribed\n"
           "\t                               cartesian region\n"
           "\tprescribed_positions: Prescribed particle positions\n"
           "\tprescribed_positions_masses: Prescribed particle positions and\n"
           "                               masses\n"
           "\tprescribed_positions_masses_tags: Prescribed particle\n"
           "                                    positions, masses, and tags\n"
           "\tuniform_inlet_continuous: Continuous injection homogeneous at\n"
           "\t                          the inlet\n"
           "\tfluxweighted_inlet_continuous: Continuous injection\n"
           "\t                               flux-weighted at the inlet\n"
           "- Initial distance from solid phase (pass only for type\n"
           "  uniform_near_solid)\n"
           "- Region boundaries (pass only for types uniform_region_cartesian\n"
           "  or fluxweighted_region_cartesian)\n"
           "- Filename (with full path) for prescribed positions, positions\n"
           "  and masses, or positions, masses, and tags (pass only for type\n"
           "  prescribed_positions, prescribed_position_masses, or\n"
           "  prescribed_positions_masses_tags, respectively)\n"
           "- Total injected mass in each injection step (pass for all types\n"
           "  except prescribed_position_masses and\n"
           "  prescribed_positions_masses_tags)\n"
           "- Initial injection time\n"
           "- Final injection time (pass only for types\n"
           "  uniform_inlet_continuous or fluxweighted_inlet_continuous)\n"
           "- Maximum timestep for continuous injection discretization in\n"
           "  units of advection time (pass only for types\n"
           "  uniform_inlet_continuous or fluxweighted_inlet_continuous)\n"
           "- Maximum timestep for continuous injection discretization in\n"
           "  units of diffusion time (pass only for types\n"
           "  uniform_inlet_continuous or fluxweighted_inlet_continuous)\n"
           "- Maximum timestep for continuous injection discretization in\n"
           "  units of reaction times (pass only for types\n"
           "  uniform_inlet_continuous or fluxweighted_inlet_continuous)\n"
           "--------------------------------------------------------------\n";
  }
};

/** \class InitialCondition_Cases PTOF/InitialCondition_Cases.h
 * "PTOF/InitialCondition_Cases.h" \brief InitialCondition object to handle
 * implemented initial condition types. */
template <typename Geometry, typename VelocityField, typename ParticleMaker,
          typename Mask = useful::Empty>
class InitialCondition_Cases {
public:
  using Parameters =
      InitialConditionParameters_Cases; /**< Initial condition parameters.*/
  Parameters const &parameters;         /**< Initial condition parameters. */

private:
  typename Geometry::Locator const &_locator; /**< Object to search mesh. */
  VelocityField const &_velocity_field; /**< Underlying flow field.       */
  ParticleMaker _particle_maker; /**< Makes a particle given a position.*/
  std::size_t _nr_particles;     /**< Number of particles per injection step. */
  std::mt19937 _rng{std::random_device{}()}; /**< Random number generator. */
  Mask _mask{};          /**< Mask to disallow cells/faces        */
  double _threshold{0.}; /**< Threshold for mask        */

  struct PositionMaker {
    double time = 0.;
    auto operator()(Foam::point const &position) { return position; }
  } _position_maker; /** Functor to forward positions and hold time value.*/

public:
  /** Constructor.
   \param geometry Domain geometry info and utilities.
   \param velocity_field Velocity field as a function of state.
   \param particle_maker Functor to make a particle given a position.
   \param nr_particles Number of particles per injection step.
   \param parameters Output parameters. */
  InitialCondition_Cases(Geometry const &geometry,
                         VelocityField const &velocity_field,
                         ParticleMaker particle_maker, std::size_t nr_particles,
                         Parameters const &parameters)
      : parameters{parameters}, _locator{geometry.locator},
        _velocity_field{velocity_field}, _particle_maker{particle_maker},
        _nr_particles{nr_particles} {}

  InitialCondition_Cases(Geometry &&, VelocityField const &, ParticleMaker,
                         std::size_t, Parameters const &) = delete;

  InitialCondition_Cases(Geometry const &, VelocityField &&, ParticleMaker,
                         std::size_t, Parameters const &) = delete;

  InitialCondition_Cases(Geometry const &, VelocityField const &, ParticleMaker,
                         std::size_t, Parameters &&) = delete;

  InitialCondition_Cases(Geometry &&, VelocityField &&, ParticleMaker,
                         std::size_t, Parameters const &) = delete;

  InitialCondition_Cases(Geometry const &, VelocityField &&, ParticleMaker,
                         std::size_t, Parameters &&) = delete;

  InitialCondition_Cases(Geometry &&, VelocityField const &, ParticleMaker,
                         std::size_t, Parameters &&) = delete;

  InitialCondition_Cases(Geometry &&, VelocityField &&, ParticleMaker,
                         std::size_t, Parameters &&) = delete;

  /** Constructor.
  \param geometry Domain geometry info and utilities.
  \param velocity_field Velocity field as a function of state.
  \param particle_maker Functor to make a particle given a position.
  \param nr_particles Number of particles per injection step.
  \param parameters Output parameters.
  \param mask Scalar field assigning values to mesh cell indices through
  operator[]. \param threshold Threshold for the mask, such that cells where the
  mask is above the threshold are considered. */
  InitialCondition_Cases(Geometry const &geometry,
                         VelocityField const &velocity_field,
                         ParticleMaker particle_maker, std::size_t nr_particles,
                         Parameters const &parameters, Mask &&mask,
                         double threshold = 0.)
      : parameters{parameters}, _locator{geometry.locator},
        _velocity_field{velocity_field}, _particle_maker{particle_maker},
        _nr_particles{nr_particles}, _mask{std::forward<Mask>(mask)},
        _threshold{threshold} {}

  InitialCondition_Cases(Geometry &&, VelocityField const &, ParticleMaker,
                         std::size_t, Parameters const &, Mask &&,
                         double = 0.) = delete;

  InitialCondition_Cases(Geometry const &, VelocityField &&_field,
                         ParticleMaker, std::size_t, Parameters const &,
                         Mask &&, double = 0.) = delete;

  InitialCondition_Cases(Geometry const &, VelocityField const &, ParticleMaker,
                         std::size_t, Parameters &&, Mask &&,
                         double = 0.) = delete;

  InitialCondition_Cases(Geometry &&, VelocityField &&, ParticleMaker,
                         std::size_t, Parameters const &, Mask &&,
                         double = 0.) = delete;

  InitialCondition_Cases(Geometry const &, VelocityField &&, ParticleMaker,
                         std::size_t, Parameters &&, Mask &&,
                         double = 0.) = delete;

  InitialCondition_Cases(Geometry &&, VelocityField const &, ParticleMaker,
                         std::size_t, Parameters &&, Mask &&,
                         double = 0.) = delete;

  InitialCondition_Cases(Geometry &&, VelocityField &&, ParticleMaker,
                         std::size_t, Parameters &&, Mask &&,
                         double = 0.) = delete;

  /** \brief Make particles according to prescribed initial condition and
   particle maker. \param nr_particles Number of particles to make \return
   Container with particles.  */
  auto operator()(std::size_t nr_particles) {
    return make_particles(nr_particles, parameters.type, _particle_maker);
  }

  /** \brief Make particles according to prescribed initial condition and
  particle maker, with \c nr_particles particles as specified in \c parameters.
  \return Container with particles.  */
  auto operator()() { return this->operator()(_nr_particles); }

  /** \brief Make a single particle according to prescribed initial condition
  and particle maker. \return One particle.  */
  auto make_particle() { return this->operator()(1)[0]; }

  /** \brief Make positions according to prescribed initial condition.
  \param nr_positions Number of positions to make
  \return Container with position.  */
  auto make_positions(std::size_t nr_positions) {
    // This works by slightly abusing make_particles with a position maker
    // instead of a particle maker, in order to produce only positions instead
    // of particles, but according to the same prescription as for producing
    // particles.
    return make_particles(nr_positions, parameters.type, _position_maker);
  }

  /** \brief Make a single position according to prescribed initial condition
  and particle maker. \return One position.  */
  auto make_position() { return make_positions(1)[0]; }

private:
  /** \brief Make given number of particles.
   \param nr_particles Number of particles to make.
   \param type Initial condition type.
   \param particle_maker Functor to make a particle given a position.
   \return Container with particles.
   */
  template <typename ParticleMakerOther>
  auto make_particles(std::size_t nr_particles, InitialConditionList::Type type,
                      ParticleMakerOther &particle_maker) {
    switch (type) {
    case InitialConditionList::Type::uniform:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return uniform(nr_particles, _locator.mesh(), _rng, particle_maker);
      else
        return uniform(nr_particles, _locator.mesh(), _rng, particle_maker,
                       _mask, _threshold);
    case InitialConditionList::Type::fluxweighted:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return fluxweighted(nr_particles, _velocity_field, _locator.mesh(),
                            _rng, particle_maker);
      else
        return fluxweighted(nr_particles, _velocity_field, _locator.mesh(),
                            _rng, particle_maker, _mask, _threshold);
    case InitialConditionList::Type::uniform_inlet:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return uniform_patches(nr_particles, {"inlet"}, _locator, _rng,
                               particle_maker);
      else
        return uniform_patches(nr_particles, {"inlet"}, _locator, _rng,
                               particle_maker, _mask, _threshold);
    case InitialConditionList::Type::fluxweighted_inlet:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return fluxweighted_patches(nr_particles, {"inlet"}, _velocity_field,
                                    _locator, _rng, particle_maker);
      else
        return fluxweighted_patches(nr_particles, {"inlet"}, _velocity_field,
                                    _locator, _rng, particle_maker, _mask,
                                    _threshold);
    case InitialConditionList::Type::uniform_solid:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return uniform_patches(nr_particles, {"wallFluidSolid"}, _locator, _rng,
                               particle_maker);
      else
        return uniform_patches(nr_particles, {"wallFluidSolid"}, _locator, _rng,
                               particle_maker, _mask, _threshold);
    case InitialConditionList::Type::uniform_near_solid:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return uniform_near_boundary_patches(nr_particles, {"wallFluidSolid"},
                                             parameters.distance_wall, _locator,
                                             _rng, particle_maker);
      else
        return uniform_near_boundary_patches(
            nr_particles, {"wallFluidSolid"}, parameters.distance_wall,
            _locator, _rng, particle_maker, _mask, _threshold);
    case InitialConditionList::Type::uniform_region_cartesian:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return uniform_cells(
            nr_particles,
            cell_ids_region_cartesian(parameters.region_boundaries, _locator),
            _locator.mesh(), _rng, particle_maker);
      else
        return uniform_cells(
            nr_particles,
            apply_mask(cell_ids_region_cartesian(parameters.region_boundaries,
                                                 _locator),
                       _mask, _threshold),
            _locator.mesh(), _rng, particle_maker);
    case InitialConditionList::Type::fluxweighted_region_cartesian:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return fluxweighted_cells(
            nr_particles,
            cell_ids_region_cartesian(parameters.region_boundaries, _locator),
            _velocity_field, _locator.mesh(), _rng, particle_maker);
      else
        return fluxweighted_cells(
            nr_particles,
            apply_mask(cell_ids_region_cartesian(parameters.region_boundaries,
                                                 _locator),
                       _mask, _threshold),
            _velocity_field, _locator.mesh(), _rng, particle_maker);
    case InitialConditionList::Type::prescribed_positions:
      return prescribed_positions(nr_particles, parameters.filename_data,
                                  particle_maker);
    case InitialConditionList::Type::prescribed_positions_masses:
      return prescribed_positions_masses(nr_particles, parameters.filename_data,
                                         particle_maker);
    case InitialConditionList::Type::prescribed_positions_masses_tags:
      return prescribed_positions_masses_tags(
          nr_particles, parameters.filename_data, particle_maker);
    case InitialConditionList::Type::uniform_inlet_continuous:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return continuous_injection(
            [this](std::size_t nr_particles,
                   ParticleMakerOther &particle_maker) {
              return uniform_patches(nr_particles, {"inlet"}, _locator, _rng,
                                     particle_maker);
            },
            nr_particles, parameters.time_min, parameters.time_max,
            parameters.time_step, particle_maker);
      else
        return continuous_injection(
            [this](std::size_t nr_particles,
                   ParticleMakerOther &particle_maker) {
              return uniform_patches(nr_particles, {"inlet"}, _locator, _rng,
                                     particle_maker, _mask, _threshold);
            },
            nr_particles, parameters.time_min, parameters.time_max,
            parameters.time_step, particle_maker);
    case InitialConditionList::Type::fluxweighted_inlet_continuous:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return continuous_injection(
            [this](std::size_t nr_particles,
                   ParticleMakerOther &particle_maker) {
              return fluxweighted_patches(nr_particles, {"inlet"},
                                          _velocity_field, _locator, _rng,
                                          particle_maker);
            },
            nr_particles, parameters.time_min, parameters.time_max,
            parameters.time_step, particle_maker);
      else
        return continuous_injection(
            [this](std::size_t nr_particles,
                   ParticleMakerOther &particle_maker) {
              return fluxweighted_patches(nr_particles, {"inlet"},
                                          _velocity_field, _locator, _rng,
                                          particle_maker, _mask, _threshold);
            },
            nr_particles, parameters.time_min, parameters.time_max,
            parameters.time_step, particle_maker);
    default:
      throw std::runtime_error{std::string{"Initial condition type "} +
                               InitialConditionList::name(parameters.type) +
                               " not supported"};
    }
  }

public:
  /** \brief Output information about current object. */
  template <typename OStream> void info_runtime(OStream &output) const {
    output
        << "--------------------------------------------------------------\n"
           "Initial condition\n"
           "--------------------------------------------------------------\n"
           "Type: " +
               InitialConditionList::name(parameters.type) + "\n"
        << "--------------------------------------------------------------\n";
  }
};
template <typename Geometry, typename VelocityField, typename ParticleMaker,
          typename Parameters, typename Mask>
InitialCondition_Cases(Geometry const &, VelocityField const &, ParticleMaker,
                       std::size_t, Parameters const &, Mask &&, double)
    -> InitialCondition_Cases<Geometry, VelocityField, ParticleMaker, Mask>;
} // namespace ptof

#endif /* PTOF_INITIALCONDITION_CASES_H */
