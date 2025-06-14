/**
   \file PTOF/InitialCondition_Cases.h
   \author Tomás Aquino
   \date 24/02/2022
   \brief Create particles according to different initial conditions.
*/

#ifndef PTOF_INITIALCONDITION_CASES_H
#define PTOF_INITIALCONDITION_CASES_H

#include "CTRW/Meta.h"
#include "General/IO.h"
#include "General/Meta.h"
#include "General/Operation.h"
#include "PTOF/Boundary.h"
#include "PTOF/InitialCondition.h"
#include "PTOF/InitialConditionList.h"
#include "PTOF/InitialConditionParameters_Cases.h"
#include "PTOF/ParticleMaker.h"
#include "PTOF/Useful.h"
#include "Stochastic/Random.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fieldTypes.H>
#include <limits>
#include <map>
#include <memory>
#include <ostream>
#include <point.H>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ptof {
/**
   \class InitialCondition_Cases PTOF/InitialCondition_Cases.h
   "PTOF/InitialCondition_Cases.h"
   \brief InitialCondition object to handle implemented initial condition types.
*/
template <typename Particle, typename Geometry> class InitialCondition_Cases {
public:
  using Parameters =
      InitialConditionParameters_Cases; /**< Initial condition parameters.*/

  std::size_t const
      nr_particles; /**< Number of particles per injection step. */
  InitialConditionList::Type type; /**< Type of initial condition. */

  /**
     \brief Constructor.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field as a function of state.
     \param nr_particles Number of particles to make.
     \param parameters Output parameters.
     \param mask Scalar field.
     \param threshold Threshold for the mask, such that cells where the mask is
     above or equal to the threshold are considered.
  */
  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(Geometry const &geometry,
                         VelocityField const &velocity_field,
                         std::size_t nr_particles, Parameters const &parameters,
                         Mask const &mask = {}, double threshold = 0.)
      : nr_particles{nr_particles}, type{parameters.type},
        _specific_parameters{*parameters.specific_parameters.get()},
        _injection_parameters{*parameters.injection_parameters.get()},
        _particle_maker{make_particle_maker(geometry, parameters)},
        _initial_condition{
            set_type(geometry, velocity_field, parameters, mask, threshold)} {}

  /**
     \brief Constructor.
     \note Specifies \c Particle type.
  */
  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(meta::Selector_t<Particle>, Geometry const &geometry,
                         VelocityField const &velocity_field,
                         std::size_t nr_particles, Parameters const &parameters,
                         Mask const &mask = {}, double threshold = 0.)
      : InitialCondition_Cases{geometry,   velocity_field, nr_particles,
                               parameters, mask,           threshold} {}

  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(Geometry &&, VelocityField const &, std::size_t,
                         Parameters const &, Mask const & = {},
                         double = 0.) = delete;

  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(Geometry const &, VelocityField const &, std::size_t,
                         Parameters &&, Mask const & = {},
                         double = 0.) = delete;

  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(Geometry &&, VelocityField const &, std::size_t,
                         Parameters &&, Mask const & = {},
                         double = 0.) = delete;

  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(meta::Selector_t<Particle>, Geometry &&,
                         VelocityField const &, std::size_t, Parameters const &,
                         Mask const & = {}, double = 0.) = delete;

  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(meta::Selector_t<Particle>, Geometry const &,
                         VelocityField &&, std::size_t, Parameters const &,
                         Mask const & = {}, double = 0.) = delete;

  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(meta::Selector_t<Particle>, Geometry &&,
                         VelocityField const &, std::size_t, Parameters &&,
                         Mask const & = {}, double = 0.) = delete;

  /**
     \brief Make particles according to prescribed initial condition and
     particle maker, with \c nr_particles particles as specified in \c
     parameters.
     \return Container with particles.
  */
  auto operator()() { return make_particles(nr_particles); }

  /**
     \brief Make particles.
     \param nr_particles Number of particles to make.
     \return Container with particles.
  */
  auto operator()(std::size_t nr_particles) {
    return make_particles(nr_particles);
  }

  /**
     \brief Make a single particle.
     \return One particle.
  */
  auto make_particle() { return _initial_condition->make_particle(); }

  /**
     \brief Make a single position.
     \return Position.
  */
  auto make_position() { return _initial_condition->make_position(); }

  /**
     \brief Make a single position along with location cell.
     \return Position and cell.
  */
  auto make_position_and_cell() {
    return _initial_condition->make_position_and_cell();
  }

  /**
     \brief Make particles.
     \param nr_particles Number of particles to make.
     \return Container with particles.
  */
  auto make_particles(std::size_t nr_particles) {
    return _initial_condition->make_particles(nr_particles);
  };

  /**
     \brief Make positions.
     \param nr_positions Number of positions to make.
     \return Container with positions.
  */
  auto make_positions(std::size_t nr_positions) {
    return _initial_condition->make_positions(nr_particles);
  };

  /**
     \brief Output information about current object.
     \param output Output stream.
  */
  std::ostream &info_runtime(std::ostream &output) const {
    output << io::line() << "Initial condition\n"
           << io::line() << "Type: " << InitialConditionList::name(type)
           << "\n";
    _injection_parameters.info_runtime(output);
    _specific_parameters.info_runtime(output) << io::line();
    return output;
  }

private:
  template <typename TT, bool periodic = false> struct BP {
    using type = meta::Empty;
  };

  template <typename TT> struct BP<TT, true> {
    using type = typename TT::BoundaryPeriodic;
  };

  static constexpr bool periodic =
      meta::has_periodicity_v<typename Particle::State>;

  using BoundaryPeriodic =
      typename std::conditional_t<periodic, BP<Geometry, true>,
                                  BP<Geometry, false>>::type;

  using ParticleMaker = std::conditional_t<
      periodic,
      ParticleMaker_Periodic<Particle, typename Geometry::Locator const &,
                             BoundaryPeriodic const &>,
      ParticleMaker_Generic<Particle, typename Geometry::Locator const &>>;
  using IC = InitialCondition<ParticleMaker &, Geometry const &>;

  Parameters::SpecificParameters &_specific_parameters;
  Parameters::InjectionParameters &_injection_parameters;
  ParticleMaker _particle_maker;
  std::unique_ptr<IC>
      _initial_condition; /**< Pointer to initial condition object. */

  ParticleMaker make_particle_maker(Geometry const &geometry,
                                    Parameters const &parameters) {
    double initial_particle_mass = 0.;
    double initial_time = 0.;
    if (parameters.continuity == "pulse" &&
        type != InitialConditionList::Type::prescribed_positions_masses &&
        type != InitialConditionList::Type::prescribed_positions_masses_tags) {
      using InjectionParameters = Parameters::InjectionParameters_Pulse;
      InjectionParameters &injection_parameters =
          downcast_injection_parameters<InjectionParameters>();
      initial_particle_mass = injection_parameters.mass / nr_particles;
      initial_time = injection_parameters.time;
    }
    if constexpr (periodic)
      return ParticleMaker{meta::Selector_t<Particle>{}, geometry.locator,
                           geometry.boundary_periodic, initial_time,
                           initial_particle_mass};
    else
      return ParticleMaker{meta::Selector_t<Particle>{}, geometry.locator,
                           initial_time, initial_particle_mass};
  }

  template <typename VelocityField, typename Mask = meta::Empty>
  std::unique_ptr<IC> set_type(Geometry const &geometry,
                               VelocityField const &velocity_field,
                               Parameters const &parameters,
                               Mask const &mask = {}, double threshold = 0.) {
    switch (type) {
    case InitialConditionList::Type::point: {
      using InitialCondition =
          InitialCondition_Point<ParticleMaker &, Geometry const &>;
      using SpecificParameters = Parameters::SpecificParameters_Position;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      return set_injection(InitialCondition{_particle_maker, geometry,
                                            ParticleMaker::State::make_position(
                                                specific_parameters.position)},
                           parameters);
    }
    case InitialConditionList::Type::uniform_all_cells: {
      using InitialCondition =
          InitialCondition_UniformCellCenters<ParticleMaker &, Geometry const &,
                                              std::vector<Foam::label>>;
      return set_injection(
          InitialCondition{_particle_maker, geometry,
                           cell_ids_masked(all_cell_ids(geometry.mesh()),
                                           geometry, mask, threshold)},
          parameters);
    }
    case InitialConditionList::Type::fluxweighted_all_cells: {
      using InitialCondition = InitialCondition_FluxweightedCellCenters<
          ParticleMaker &, Geometry const &, std::vector<Foam::label>>;
      if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
        return set_injection(
            InitialCondition{_particle_maker, geometry,
                             cell_ids_masked(all_cell_ids(geometry.mesh()),
                                             geometry, mask, threshold),
                             velocity_field},
            parameters);
      else
        throw std::runtime_error{std::string("Initial condition type ") +
                                 InitialConditionList::name(type) + ": " +
                                 "Velocity field not provided"};
    }
    case InitialConditionList::Type::uniform_patch_faces: {
      using InitialCondition =
          InitialCondition_UniformFaceCenters<ParticleMaker &, Geometry const &,
                                              std::vector<Foam::label>>;
      using SpecificParameters = Parameters::SpecificParameters_Patches;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      return set_injection(
          InitialCondition{
              _particle_maker, geometry,
              patch_face_ids_masked(specific_parameters.patch_names, geometry,
                                    mask, threshold)},
          parameters);
    }
    case InitialConditionList::Type::fluxweighted_patch_faces: {
      using InitialCondition = InitialCondition_FluxweightedFaceCenters<
          ParticleMaker &, Geometry const &, std::vector<Foam::label>>;
      using SpecificParameters = Parameters::SpecificParameters_Patches;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
        return set_injection(
            InitialCondition{
                _particle_maker, geometry,
                patch_face_ids_masked(specific_parameters.patch_names,
                                      geometry, mask, threshold),
                velocity_field},
            parameters);
      else
        throw std::runtime_error{std::string("Initial condition type ") +
                                 InitialConditionList::name(type) + ": " +
                                 "Velocity field not provided"};
    }
    case InitialConditionList::Type::uniform_near_patch: {
      using InitialCondition =
          InitialCondition_UniformNearFaces<ParticleMaker &, Geometry const &,
                                            std::vector<Foam::label>>;
      using SpecificParameters =
          Parameters::SpecificParameters_Patches_Distance;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      return set_injection(
          InitialCondition{
              _particle_maker, geometry,
              patch_face_ids_masked(specific_parameters.patch_names, geometry,
                                    mask, threshold),
              specific_parameters.distance, specific_parameters.nr_tries},
          parameters);
    }
    case InitialConditionList::Type::fluxweighted_near_patch: {
      using InitialCondition = InitialCondition_FluxweightedNearFaces<
          ParticleMaker &, Geometry const &, std::vector<Foam::label>>;
      using SpecificParameters =
          Parameters::SpecificParameters_Patches_Distance;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
        return set_injection(
            InitialCondition{
                _particle_maker, geometry,
                patch_face_ids_masked(specific_parameters.patch_names,
                                      geometry, mask, threshold),
                velocity_field, specific_parameters.distance,
                specific_parameters.nr_tries},
            parameters);
      else
        throw std::runtime_error{std::string("Initial condition type ") +
                                 InitialConditionList::name(type) + ": " +
                                 "Velocity field not provided"};
    }
    case InitialConditionList::Type::uniform_cartesian_cells: {
      using InitialCondition =
          InitialCondition_UniformCellCenters<ParticleMaker &, Geometry const &,
                                              std::vector<Foam::label>>;
      using SpecificParameters = Parameters::SpecificParameters_Boundaries;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      return set_injection(
          InitialCondition{
              _particle_maker, geometry,
              cell_ids_masked(
                  cell_ids_region_cartesian(
                      specific_parameters.region_boundaries, geometry.locator),
                  geometry, mask, threshold)},
          parameters);
    }
    case InitialConditionList::Type::fluxweighted_cartesian_cells: {
      using InitialCondition = InitialCondition_FluxweightedCellCenters<
          ParticleMaker &, Geometry const &, std::vector<Foam::label>>;
      using SpecificParameters = Parameters::SpecificParameters_Boundaries;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
        return set_injection(
            InitialCondition{
                _particle_maker, geometry,
                cell_ids_masked(cell_ids_region_cartesian(
                                    specific_parameters.region_boundaries,
                                    geometry.locator),
                                geometry, mask, threshold),
                velocity_field},
            parameters);
      else
        throw std::runtime_error{std::string("Initial condition type ") +
                                 InitialConditionList::name(type) + ": " +
                                 "Velocity field not provided"};
    }
    case InitialConditionList::Type::prescribed_positions: {
      using InitialCondition =
          InitialCondition_PrescribedPositions<ParticleMaker &,
                                               Geometry const &>;
      using SpecificParameters = Parameters::SpecificParameters_Filename;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      return set_injection(InitialCondition{_particle_maker, geometry,
                                            specific_parameters.filename},
                           parameters);
    }
    case InitialConditionList::Type::prescribed_positions_masses: {
      using InitialCondition =
          InitialCondition_PrescribedPositionsMasses<ParticleMaker &,
                                                     Geometry const &>;
      using SpecificParameters = Parameters::SpecificParameters_Filename;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      return set_injection(InitialCondition{_particle_maker, geometry,
                                            specific_parameters.filename},
                           parameters);
    }
    case InitialConditionList::Type::prescribed_positions_masses_tags: {
      using InitialCondition =
          InitialCondition_PrescribedPositionsMassesTags<ParticleMaker &,
                                                         Geometry const &>;
      using SpecificParameters = Parameters::SpecificParameters_Filename;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      return set_injection(InitialCondition{_particle_maker, geometry,
                                            specific_parameters.filename},
                           parameters);
    }
    case InitialConditionList::Type::prescribed_positions_masses_tags_times: {
      using InitialCondition =
          InitialCondition_PrescribedPositionsMassesTagsTimes<ParticleMaker &,
                                                              Geometry const &>;
      using SpecificParameters = Parameters::SpecificParameters_Filename;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      return set_injection(InitialCondition{_particle_maker, geometry,
                                            specific_parameters.filename},
                           parameters);
    }
    case InitialConditionList::Type::uniform_cylinder: {
      using Distribution = stochastic::uniform_cylinder_distribution;
      using InitialCondition =
          InitialCondition_DistributedPosition<ParticleMaker &,
                                               Geometry const &, Distribution>;
      using SpecificParameters = Parameters::SpecificParameters_Cylinder;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      return set_injection(
          InitialCondition{_particle_maker, geometry,
                           Distribution{specific_parameters.center_1,
                                        specific_parameters.center_2,
                                        specific_parameters.radius,
                                        specific_parameters.inner_radius},
                           specific_parameters.nr_tries},
          parameters);
    }
    case InitialConditionList::Type::uniform_parallelipiped: {
      using Distribution = stochastic::uniform_parallelipiped_distribution;
      using InitialCondition =
          InitialCondition_DistributedPosition<ParticleMaker &,
                                               Geometry const &, Distribution>;
      using SpecificParameters = Parameters::SpecificParameters_Parallelipiped;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      return set_injection(
          InitialCondition{_particle_maker, geometry,
                           Distribution{specific_parameters.corner,
                                        specific_parameters.sides},
                           specific_parameters.nr_tries},
          parameters);
    }
    case InitialConditionList::Type::uniform_sphere: {
      using Distribution = stochastic::uniform_sphere_distribution;
      using InitialCondition =
          InitialCondition_DistributedPosition<ParticleMaker &,
                                               Geometry const &, Distribution>;
      using SpecificParameters = Parameters::SpecificParameters_Sphere;
      SpecificParameters &specific_parameters =
          downcast_specific_parameters<SpecificParameters>();
      return set_injection(
          InitialCondition{_particle_maker, geometry,
                           Distribution{specific_parameters.center,
                                        specific_parameters.radius,
                                        specific_parameters.inner_radius},
                           specific_parameters.nr_tries},
          parameters);
    }

    default:
      throw std::runtime_error{std::string{"Initial condition type "} +
                               InitialConditionList::name(parameters.type) +
                               " not supported"};
    }
  }

  template <typename Container, typename Mask = meta::Empty>
  auto cell_ids_masked(Container cell_ids, Geometry const &geometry,
                       Mask const &mask = {}, double threshold = 0.) {
    if constexpr (!std::is_same_v<Mask, meta::Empty>)
      apply_mask_cells_inplace(cell_ids, mask, threshold);
    return cell_ids;
  }

  template <typename Container, typename Mask = meta::Empty>
  auto patch_face_ids_masked(Container face_ids, Geometry const &geometry,
                             Mask const &mask = {}, double threshold = 0.) {
    if constexpr (!std::is_same_v<Mask, meta::Empty>)
      apply_mask_patch_faces_inplace(face_ids, geometry.mesh(), mask,
                                     threshold);
    return face_ids;
  }

  template <typename Mask = meta::Empty>
  auto patch_face_ids_masked(std::vector<std::string> const &patch_names,
                             Geometry const &geometry, Mask const &mask = {},
                             double threshold = 0.) {
    return patch_face_ids_masked(patch_face_ids(patch_names, geometry.mesh()),
                                 geometry);
  }

  template <typename InitialCondition>
  std::unique_ptr<IC> set_injection(InitialCondition ic,
                                    Parameters const &parameters) {
    if (parameters.continuity == "continuous") {
      using InjectionParameters =
          Parameters::InjectionParameters_ContinuousNoMass;
      InjectionParameters &injection_parameters =
          downcast_injection_parameters<InjectionParameters>();
      return std::make_unique<InitialCondition_Continuous<InitialCondition>>(
          std::move(ic), injection_parameters.time,
          injection_parameters.time_end, injection_parameters.time_step);
    }
    return std::make_unique<InitialCondition>(std::move(ic));
  }

  template <typename SpecificParameters>
  SpecificParameters &downcast_specific_parameters() {
    return dynamic_cast<SpecificParameters &>(_specific_parameters);
  }

  template <typename InjectionParameters>
  InjectionParameters &downcast_injection_parameters() {
    return dynamic_cast<InjectionParameters &>(_injection_parameters);
  }
};
template <typename Particle, typename Geometry, typename VelocityField,
          typename Mask>
InitialCondition_Cases(meta::Selector_t<Particle>, Geometry const &,
                       VelocityField const &, std::size_t nr_particles,
                       InitialConditionParameters_Cases const &, Mask const &,
                       double) -> InitialCondition_Cases<Particle, Geometry>;
template <typename Particle, typename Geometry, typename VelocityField,
          typename Mask>
InitialCondition_Cases(meta::Selector_t<Particle>, Geometry const &,
                       VelocityField const &,
                       InitialConditionParameters_Cases const &, Mask const &)
    -> InitialCondition_Cases<Particle, Geometry>;
template <typename Particle, typename Geometry, typename VelocityField>
InitialCondition_Cases(meta::Selector_t<Particle>, Geometry const &,
                       VelocityField const &, std::size_t nr_particles,
                       InitialConditionParameters_Cases const &)
    -> InitialCondition_Cases<Particle, Geometry>;
} // namespace ptof

#endif /* PTOF_INITIALCONDITION_CASES_H */
