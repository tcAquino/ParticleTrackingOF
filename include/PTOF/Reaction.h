/**
   \file PTOF/Reaction.h
   \author Tomás Aquino
   \date 10/03/2022
   \brief Objects and utilities for surface and bulk reactions.
*/

#ifndef PTOF_REACTION_H
#define PTOF_REACTION_H

#include "General/Constants.h"
#include "General/IO.h"
#include "General/Parallel.h"
#include "PTOF/Useful.h"
#include <cmath>
#include <fieldTypes.H>
#include <functional>
#include <list>
#include <ostream>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ptof {
/**
   \brief Append homogeneous surface concentrations over patch faces to map.
   \param face_ids Mesh face indices where to assign concentration.
   \param surface_concentration Homogeneous surface concentration value.
   \param mesh Mesh object.
   \param surface_concentrations Map of surface concentrations over given face
   indices, to append result.
*/
template <typename Container, typename Mesh>
void uniform_solid_reactant(
    Container const &face_ids, double surface_concentration, Mesh const &mesh,
    std::unordered_map<Foam::label, double> &surface_concentrations) {
  for (auto face : face_ids)
    surface_concentrations[face] = surface_concentration;
}

/**
   \param patch_names Mesh patch names where to assign concentration.
   \param surface_concentration_values Surface concentration value in each
   patch.
   \param mesh Mesh object.
   \return Map of homogeneous surface concentrations over patch face indices.
*/
template <typename Mesh>
auto uniform_solid_reactant(
    std::vector<std::string> const &patch_names,
    std::vector<double> const &surface_concentration_values, Mesh const &mesh) {
  std::unordered_map<Foam::label, double> surface_concentrations;
  for (std::size_t ii = 0; ii < patch_names.size(); ++ii)
    uniform_solid_reactant(patch_face_ids({patch_names[ii]}, mesh),
                           surface_concentration_values[ii], mesh,
                           surface_concentrations);
  return surface_concentrations;
}

/**
   \param position Spatial position.
   \param locator Object to locate positions in mesh.
   \return Pair of nearest boundary face index and distance to it.
*/
template <typename Locator>
std::pair<Foam::label, double>
nearest_boundary_face_dist(Foam::vector const &position,
                           Locator const &locator) {
  auto const &mesh = locator.mesh();
  // Do not use hint because this can get stuck in local minima
  auto face_id = locator.mesh_search().findNearestBoundaryFace(position);
  auto dist = Foam::mag(position - mesh.faces()[face_id].centre(mesh.points()));
  return {face_id, dist};
}

/**
   \class SurfaceReaction_AFluidPlusASolidtoASolid PTOF/Reaction.h
   "PTOF/Reaction.h"
   \brief \f$A_F + A_S \to A_S\f$ surface reaction.
*/
class SurfaceReaction_AFluidPlusASolidtoASolid {
public:
  /** Container to hold reactant surface concentrations. */
  using SurfaceConcentrations = std::unordered_map<Foam::label, double>;

  /**
     \brief Constructor.
     \param rate_constant Surface reaction rate per solid concentration.
     \param diff_coeff Diffusion coefficient.
     \param surface_concentrations Map of initial reactive surface
     concentrations over mesh faces.
  */
  SurfaceReaction_AFluidPlusASolidtoASolid(
      double rate_constant, double diff_coeff,
      SurfaceConcentrations surface_concentrations)
      : rate_constant{rate_constant}, diff_coeff{diff_coeff},
        surface_concentrations{surface_concentrations} {}

  /**
     \brief React over exposure time.
     \param state State of particle to react.
     \param state_old Previous state of particle.
     \param face Mesh face index where reaction occurs.
  */
  template <typename State>
  void operator()(State &state, State const &state_old, Foam::label face) {
    double exposure_time = state.time - state_old.time;
    double probability_first_order =
        rate(face) * std::sqrt(constants::pi * exposure_time / diff_coeff);
    double probability_second_order =
        probability_first_order / (1. + probability_first_order / 2.);
    state.mass *= std::exp(-probability_second_order);
  }

  /**
     \param face Mesh face index.
     \return Surface reaction rate at face.
  */
  double rate(Foam::label face) const {
    if (surface_concentrations.count(face))
      return rate_constant * surface_concentrations.at(face);
    return 0.;
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    return output << io::line() << "Surface reaction\n"
                  << io::line() << "A_Fluid + A_Solid -> A_Solid\n"
                  << io::line();
  }

  double rate_constant; /**< Surface reaction rate per solid concentration. */
  double diff_coeff;    /**< Diffusion coefficient. */
  SurfaceConcentrations
      surface_concentrations; /**< Map of concentrations over reactive faces. */
};

/**
   \class SurfaceReaction_AFluidPlusASolidtoNothing PTOF/Reaction.h
   "PTOF/Reaction.h"
   \brief \f$A_F + A_S \to \emptyset\f$ surface reaction.
*/
template <typename Locator, typename ParallelOption>
class SurfaceReaction_AFluidPlusASolidtoNothing {
public:
  /** Container to hold reactant surface concentrations. */
  using SurfaceConcentrations = std::unordered_map<Foam::label, double>;

  /**
     \brief Constructor.
     \param rate_constant Surface reaction rate per solid concentration
     for fluid phase.
     \param rate_constant_ratio_solid_to_fluid Ratio of reactions for solid and
     fluid phase concentrations.
     \param diff_coeff Diffusion coefficient.
     \param surface_concentrations Map of initial reactive surface
     concentrations over mesh faces.
     \param locator Object to locate points in mesh.
  */
  SurfaceReaction_AFluidPlusASolidtoNothing(
      double rate_constant, double rate_constant_ratio_solid_to_fluid,
      double diff_coeff, SurfaceConcentrations surface_concentrations,
      Locator &&locator, ParallelOption)
      : rate_constant{rate_constant},
        rate_constant_ratio_solid_to_fluid{rate_constant_ratio_solid_to_fluid},
        solid_rate_constant{rate_constant_ratio_solid_to_fluid * rate_constant},
        diff_coeff{diff_coeff}, surface_concentrations{surface_concentrations},
        _locator{std::forward<Locator>(locator)},
        _face_to_tags_times(par::get_num_threads(ParallelOption{})) {}

  /**
     \brief React over exposure time.
     \param state State of particle to react.
     \param state_old Previous state of particle.
     \param face Mesh face index where reaction occurs.
  */
  template <typename State>
  void operator()(State &state, State const &state_old, Foam::label face) {
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Serial>) {
      update_mass_and_surface_concentration(state, state_old, face);
    } else {
      _face_to_tags_times[par::get_thread_num(ParallelOption{})][face].insert(
          {state.tag, state.time});
    }
  }

  /**
     \param face Mesh face index.
     \return Surface reaction rate for fluid phase at face.
  */
  double rate(Foam::label face) const {
    if (surface_concentrations.count(face))
      return rate_constant * surface_concentrations.at(face);
    return 0.;
  }

  /**
     \param face Mesh face index.
     \return Surface reaction rate for solid phase at face.
  */
  double solid_rate(Foam::label face) const {
    if (surface_concentrations.count(face))
      return solid_rate_constant * surface_concentrations.at(face);
    return 0.;
  }

  /**
     \param face Mesh face index.
     \return Surface reaction rate for fluid and solid phases at face.
  */
  std::vector<double> rates(Foam::label face) const {
    if (surface_concentrations.count(face)) {
      double fluid_rate = solid_rate_constant * surface_concentrations.at(face);
      return {fluid_rate, rate_constant_ratio_solid_to_fluid * fluid_rate};
    }
    return {0., 0.};
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    return output << io::line() << "Surface reaction\n"
                  << io::line() << "A_Fluid + A_Solid -> Nothing\n"
                  << io::line();
  }

  /**
     \brief Update masses and surfaces concentrations according to stored
     particle-boundary collisions.
  **/
  template <typename CTRW>
  void update(CTRW &ctrw, typename CTRW::State::Time time_new,
              typename CTRW::State::Time time_old) {
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Parallel>) {
      auto merged_face_to_tags_times = merge_and_clear_hits();
#ifdef _OPENMP
#pragma omp parallel for
#endif /** _OPENMP */
      for (std::size_t ii = 0; ii < merged_face_to_tags_times.size(); ++ii) {
        // Note: Must advance iterator manually each time due to current OMP
        // limitations on some platforms
        auto it = merged_face_to_tags_times.begin();
        std::advance(it, ii);
        auto const &[face, tags_times] = *it;
        for (auto const &tag_time : tags_times) {
          // Due to a limitation of the current language standard, lambdas
          // cannot capture elements of structured bindings
          auto const &face_to_capture = face;
          ctrw.transform(
              [this, face_to_capture](typename CTRW::State state,
                                      typename CTRW::State const &state_old) {
                update_mass_and_surface_concentration(state, state_old,
                                                      face_to_capture);
              },
              tag_time.tag);
        }
      }
    }
  }

  double rate_constant; /**< Surface reaction rate per solid
                           concentration for fluid phase. */
  double rate_constant_ratio_solid_to_fluid; /** Ratio of solid to fluid rate
                                                constants. */
  double solid_rate_constant;                /**< Surface reaction rate per
                                                solid concentration for solid phase. */
  double diff_coeff;                         /**< Diffusion coefficient. */
  SurfaceConcentrations
      surface_concentrations; /**< Map of concentrations over reactive faces. */

private:
  Locator _locator; /**< To locate positions in mesh.*/

  struct TagTime {
    std::size_t tag;
    double time;
  }; /**< Particle tag and time to record face hits. */

  /**
     \struct TagTime_SetComparator
     \brief When used as a comparator function with std::set, this orders by
     times but considers elements equal if they have equal tags.
  */
  struct TagTime_SetComparator {
    bool operator()(TagTime const &v1, TagTime const &v2) const {
      return (v1.time < v2.time) & (v1.tag != v2.tag);
    }
  };

  using TagTimeContainer = std::set<TagTime, TagTime_SetComparator>;
  using FaceTagTimeContainer =
      std::unordered_map<Foam::label, TagTimeContainer>;

  std::vector<FaceTagTimeContainer>
      _face_to_tags_times; /**< Associative containers of particle tags and
                              times of hitting faces, one per thread. */

  /**
     \brief Merge collision data across threads.
     \note Merge maps, but keep only one face hit per particle for
     parallelization reasons. This is fine in the limit of small time step.
  */
  auto merge_and_clear_hits() {
    FaceTagTimeContainer merged_face_to_tags_times;
    for (auto &container : _face_to_tags_times) {
      for (auto &face_tags_times : container) {
        // Check if face is already in thread-merged container
        auto it = merged_face_to_tags_times.find(face_tags_times.first);
        if (it == std::end(merged_face_to_tags_times)) {
          // If it isn't, insert it, copying in all hits
          merged_face_to_tags_times.insert(it, face_tags_times);
        } else {
          // If it is, insert additional hits
          it->second.merge(face_tags_times.second);
        }
      }
      container.clear();
    }

    return merged_face_to_tags_times;
  }

  template <typename State>
  void update_mass_and_surface_concentration(State state,
                                             State const &state_old,
                                             Foam::label face) {
    double exposure_time = state.time - state_old.time;
    double probability_first_order =
        rate(face) * std::sqrt(constants::pi * exposure_time / diff_coeff);
    double probability_second_order =
        probability_first_order / (1. + probability_first_order / 2.);
    double delta_mass =
        state_old.mass * (std::exp(-probability_second_order) - 1.);
    double mass_solid =
        surface_concentrations[face] * face_area(face, _locator.mesh());
    if (mass_solid < rate_constant_ratio_solid_to_fluid * delta_mass) {
      state.mass -= mass_solid / rate_constant_ratio_solid_to_fluid;
      surface_concentrations[face] = 0.;
    } else {
      state.mass += delta_mass;
      surface_concentrations[face] -= rate_constant_ratio_solid_to_fluid *
                                      delta_mass /
                                      face_area(face, _locator.mesh());
    }
  }
};
template <typename Locator, typename ParallelOption>
SurfaceReaction_AFluidPlusASolidtoNothing(
    double, double, double, std::unordered_map<Foam::label, double>, Locator &&,
    ParallelOption)
    -> SurfaceReaction_AFluidPlusASolidtoNothing<Locator, ParallelOption>;

/**
   \class SurfaceReaction_DoNothing PTOF/Reaction.h "PTOF/Reaction.h"
   \brief Surface reaction that does nothing.
*/
class SurfaceReaction_DoNothing {
public:
  /**
     \param state State of particle to react (unused).
     \param state_old Previous state of particle (unused).
     \param face Mesh face index where reaction occurs (unused).
     \brief React over exposure time (do nothing).
  */
  template <typename State>
  void operator()(State &state, State const &state_old,
                  Foam::label face) const {}

  /**
     \param face Mesh face index (unused).
     \return Surface reaction rate at face (always 0.).
  */
  double rate(Foam::label face) const { return 0.; }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    return output << io::line() << "Surface reaction\n"
                  << io::line() << "None\n"
                  << io::line();
  }
};

/**
   \class BulkReaction_DoNothing PTOF/Reaction.h "PTOF/Reaction.h"
   \brief Bulk reaction that does nothing.
*/
class BulkReaction_DoNothing {
public:
  /**
     \brief React over exposure time (do nothing).
     \param state State of particle to react (unused).
     \param exposure_time Exposure time over which to react (unused).
  */
  template <typename State>
  void operator()(State &state, double exposure_time) const {}

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    return output << io::line() << "Bulk reaction\n"
                  << io::line() << "None\n"
                  << io::line();
  }
};
} // namespace ptof

#endif /* PTOF_REACTION_H */
