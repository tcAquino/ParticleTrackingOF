/**
 * @file   SurfaceReaction.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Thu Mar 10 00:00:00 2022

 * @brief Objects and utilities for surface and bulk reactions.
*/

#ifndef PTOF_SURFACEREACTION_H
#define PTOF_SURFACEREACTION_H

#include "CTRW/TimeGenerator.h"
#include "General/Constant.h"
#include "General/IO.h"
#include "General/Parallel.h"
#include "PTOF/Useful.h"
#include "Stochastic/Random.h"
#include <cmath>
#include <fieldTypes.H>
#include <ostream>
#include <random>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ptof {
/**
 * @brief Append homogeneous surface concentrations over patch faces to map.
 *
 * @param face_ids Mesh face indices where to assign concentration.
 *
 * @param surface_concentration Homogeneous surface concentration value.
 *
 * @param mesh Mesh object.
 *
 * @param surface_concentrations Map of surface concentrations over given face
 *                               indices, to append result.
 */
template <typename Container, typename Mesh>
void uniform_solid_reactant(
    Container const &face_ids, double surface_concentration, Mesh const &mesh,
    std::unordered_map<Foam::label, double> &surface_concentrations) {
  for (auto face : face_ids) {
    surface_concentrations[face] = surface_concentration;
  }
}

/**
 * @param patch_names Mesh patch names where to assign concentration.
 *
 * @param surface_concentration_values Surface concentration value in each
 *                                     patch.
 *
 * @param mesh Mesh object.
 *
 * @return Map of homogeneous surface concentrations over patch face indices.
 */
template <typename Mesh>
auto uniform_solid_reactant(
    std::vector<std::string> const &patch_names,
    std::vector<double> const &surface_concentration_values, Mesh const &mesh) {
  std::unordered_map<Foam::label, double> surface_concentrations;
  for (std::size_t ii = 0; ii < patch_names.size(); ++ii) {
    uniform_solid_reactant(patch_face_ids({patch_names[ii]}, mesh),
                           surface_concentration_values[ii], mesh,
                           surface_concentrations);
  }
  return surface_concentrations;
}

/**
 * @brief Append surface concentrations over face ids to map.
 *
 * @param face_ids Mesh face indices where to assign concentration.
 *
 * @param surface_concentration_values Surface concentration values.
 *
 * @param surface_concentrations Map of surface concentrations over given face
 *                               indices, to append result.
 */
template <typename Container>
void solid_reactant(
    Container const &face_ids,
    std::vector<double> const &surface_concentration_values,
    std::unordered_map<Foam::label, double> &surface_concentrations) {
  for (std::size_t ii = 0; ii < surface_concentration_values.size(); ++ii) {
    surface_concentrations[face_ids[ii]] = surface_concentration_values[ii];
  }
}

/**
 * @param face_ids Mesh face indices where to assign concentration.
 *
 * @param surface_concentration_values Surface concentration values.
 *
 * @return Map of homogeneous surface concentrations over patch face indices.
 */
template <typename Container>
auto solid_reactant(Container const &face_ids,
                    std::vector<double> const &surface_concentration_values) {
  std::unordered_map<Foam::label, double> surface_concentrations;
  solid_reactant(face_ids, surface_concentration_values,
                 surface_concentrations);
  return surface_concentrations;
}

/** @brief \f$A_F + A_S \to A_S\f$ surface reaction. */
class SurfaceReaction_AFluidPlusASolidtoASolid {
public:
  /** Container to hold reactant surface concentrations. */
  using SurfaceConcentrations = std::unordered_map<Foam::label, double>;

  /**
   * @brief Constructor.
   *
   * @param rate_constant Surface reaction rate per solid concentration.
   *
   * @param diff_coeff Diffusion coefficient.
   *
   * @param surface_concentrations Map of initial reactive surface
   *                               concentrations over mesh faces.
   */
  SurfaceReaction_AFluidPlusASolidtoASolid(
      double rate_constant, double diff_coeff,
      SurfaceConcentrations surface_concentrations)
      : rate_constant{rate_constant}, diff_coeff{diff_coeff},
        surface_concentrations{surface_concentrations} {}

  /**
   * @brief React over exposure time.
   *
   * @param state State of particle to react.
   *
   * @param state_old Previous state of particle.
   *
   * @param face Mesh face index where reaction occurs.
   */
  template <typename State>
  void operator()(State &state, State const &state_old, Foam::label face) {
    double exposure_time = state.time - state_old.time;
    double probability_first_order =
        rate(face) * std::sqrt(cnst::pi * exposure_time / diff_coeff);
    double probability_second_order =
        probability_first_order / (1. + probability_first_order / 2.);
    state.mass *= std::exp(-probability_second_order);
  }

  /**
   * @param face Mesh face index.
   *
   * @return Surface reaction rate at face.
   */
  double rate(Foam::label face) const {
    if (surface_concentrations.count(face)) {
      return rate_constant * surface_concentrations.at(face);
    }
    return 0.;
  }

  /**
   * @brief Output generic information about object.
   *
   * @param output Output stream.
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
 * @brief \f$A_F + A_S \leftrightarrow A_S\f$ surface reaction with distributed
 *        desorption (back-reaction) time.
 */
template <typename ParallelOption, typename TimeGenerator_Desorption,
          typename Engine = std::mt19937>
class SurfaceReaction_LinearSorption {
public:
  /** Container to hold reactant surface concentrations. */
  using SurfaceConcentrations = std::unordered_map<Foam::label, double>;

  /**
   * @brief Constructor.
   *
   * @param rate_constant Surface reaction rate per solid concentration.
   *
   * @param diff_coeff Diffusion coefficient.
   *
   * @param surface_concentrations Map of initial reactive surface
   *                               concentrations over mesh faces.
   *
   * @param desorption_time Time generator object for desorption times.
   */
  SurfaceReaction_LinearSorption(double rate_constant, double diff_coeff,
                                 SurfaceConcentrations surface_concentrations,
                                 TimeGenerator_Desorption &&desorption_time)
      : rate_constant{rate_constant}, diff_coeff{diff_coeff},
        surface_concentrations{surface_concentrations},
        _desorption_time{
            std::forward<TimeGenerator_Desorption>(desorption_time)} {}

  /**
   * @brief Constructor.
   *
   * @param rate_constant Surface reaction rate per solid concentration.
   *
   * @param diff_coeff Diffusion coefficient.
   *
   * @param surface_concentrations Map of initial reactive surface
   *                               concentrations over mesh faces.
   *
   * @param desorption_rate Exponential desorption time rate.
   */
  SurfaceReaction_LinearSorption(double rate_constant, double diff_coeff,
                                 SurfaceConcentrations surface_concentrations,
                                 double desorption_rate)
      : SurfaceReaction_LinearSorption{
            rate_constant, diff_coeff, surface_concentrations,
            ctrw::TimeGenerator_Dist{
                std::exponential_distribution<double>{desorption_rate},
                ParallelOption{}}} {}

  SurfaceReaction_LinearSorption(double rate_constant, double diff_coeff,
                                 SurfaceConcentrations surface_concentrations,
                                 TimeGenerator_Desorption &&desorption_time,
                                 ParallelOption)
      : SurfaceReaction_LinearSorption{
            rate_constant, diff_coeff, surface_concentrations,
            std::forward<TimeGenerator_Desorption>(desorption_time)} {}

  SurfaceReaction_LinearSorption(double rate_constant, double diff_coeff,
                                 SurfaceConcentrations surface_concentrations,
                                 double desorption_rate, ParallelOption)
      : SurfaceReaction_LinearSorption{rate_constant, diff_coeff,
                                       surface_concentrations,
                                       desorption_rate} {}

  /**
   * @brief React over exposure time.
   *
   * @param state State of particle to react.
   *
   * @param state_old Previous state of particle.
   *
   * @param face Mesh face index where reaction occurs.
   */
  template <typename State>
  void operator()(State &state, State const &state_old, Foam::label face) {
    double exposure_time = state.time - state_old.time;
    double probability_first_order =
        rate(face) * std::sqrt(cnst::pi * exposure_time / diff_coeff);
    double probability_second_order =
        probability_first_order / (1. + probability_first_order / 2.);
    if (_uniform_dist(_rng) < probability_second_order) {
      state.info.adsorbed = true;
    }
  }

  /** @brief Desorb after distributed time. */
  template <typename State> void desorb(State &state) {
    state.time += _desorption_time(state);
    state.info.adsorbed = false;
  }

  /**
   * @param face Mesh face index.
   *
   * @return Surface reaction rate at face.
   */
  double rate(Foam::label face) const {
    if (surface_concentrations.count(face)) {
      return rate_constant * surface_concentrations.at(face);
    }
    return 0.;
  }

  /**
   * @brief Output generic information about object.
   *
   * @param output Output stream.
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

private:
  TimeGenerator_Desorption _desorption_time;
  std::uniform_real_distribution<double> _uniform_dist{0., 1.};
  stochastic::RNGThreaded<ParallelOption> _rng{};
};
template <typename TimeGenerator_Desorption>
SurfaceReaction_LinearSorption(double, double,
                               std::unordered_map<Foam::label, double>,
                               TimeGenerator_Desorption &&)
    -> SurfaceReaction_LinearSorption<par::ParallelOptions::Serial,
                                      TimeGenerator_Desorption, std::mt19937>;
template <typename ParallelOption, typename TimeGenerator_Desorption>
SurfaceReaction_LinearSorption(double, double,
                               std::unordered_map<Foam::label, double>,
                               TimeGenerator_Desorption &&, ParallelOption)
    -> SurfaceReaction_LinearSorption<ParallelOption, TimeGenerator_Desorption,
                                      std::mt19937>;

/** @brief \f$A_F + A_S \to \emptyset\f$ surface reaction. */
template <typename Locator, typename ParallelOption>
class SurfaceReaction_AFluidPlusASolidtoNothing {
public:
  /** Container to hold reactant surface concentrations. */
  using SurfaceConcentrations = std::unordered_map<Foam::label, double>;

  /**
   * @brief Constructor.
   *
   * @param rate_constant Surface reaction rate per solid concentration
   *                      for fluid phase.
   *
   * @param rate_constant_ratio_solid_to_fluid Ratio of reactions for solid and
   *                                           fluid phase concentrations.
   *
   * @param diff_coeff Diffusion coefficient.
   *
   * @param surface_concentrations Map of initial reactive surface
   *                               concentrations over mesh faces.
   *
   * @param locator Object to locate points in mesh.
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
        _face_to_tags(par::get_num_threads(ParallelOption{})) {}

  /**
   * @brief React over exposure time.
   *
   * @param state State of particle to react.
   *
   * @param state_old Previous state of particle.
   *
   * @param face Mesh face index where reaction occurs.
   */
  template <typename State>
  void operator()(State &state, State const &state_old, Foam::label face) {
    if (!surface_concentrations.count(face)) {
      return;
    }
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Serial>) {
      update_mass_and_surface_concentration(state, state_old, face);
    } else {
      _face_to_tags[par::get_thread_num(ParallelOption{})][face].insert(
          {state.tag, state.time});
    }
  }

  /**
   * @param face Mesh face index.
   *
   * @return Surface reaction rate for fluid phase at face.
   */
  double rate(Foam::label face) const {
    if (surface_concentrations.count(face)) {
      return rate_constant * surface_concentrations.at(face);
    }
    return 0.;
  }

  /**
   * @param face Mesh face index.
   *
   * @return Surface reaction rate for solid phase at face.
   */
  double solid_rate(Foam::label face) const {
    if (surface_concentrations.count(face)) {
      return solid_rate_constant * surface_concentrations.at(face);
    }
    return 0.;
  }

  /**
   * @param face Mesh face index.
   *
   * @return Surface reaction rate for fluid and solid phases at face.
   */
  std::vector<double> rates(Foam::label face) const {
    if (surface_concentrations.count(face)) {
      double fluid_rate = solid_rate_constant * surface_concentrations.at(face);
      return {fluid_rate, rate_constant_ratio_solid_to_fluid * fluid_rate};
    }
    return {0., 0.};
  }

  /**
   * @brief Output generic information about object.
   *
   * @param output Output stream.
   */
  inline static std::ostream &info(std::ostream &output) {
    return output << io::line() << "Surface reaction\n"
                  << io::line() << "A_Fluid + A_Solid -> Nothing\n"
                  << io::line();
  }

  /**
   * @brief Update masses and surface concentrations according to stored
   *        particle-boundary collisions for parallel simulations.
   *
   * @details For serial simulations, this does nothing as updates are made
   *          on-the-fly.
   *
   * @note
   *
   * - Requires particle tags to correspond to their position in the CTRW
   *   container.
   *
   * - Requires limiting collisions to one collision per particle.
   **/
  template <typename CTRW>
  void update(CTRW &ctrw, typename CTRW::State::Time time_new,
              typename CTRW::State::Time time_old) {
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Serial>) {
      return;
    }

    auto merged_face_to_tags = merge_and_clear_hits();
#ifdef _OPENMP
#pragma omp parallel for
#endif /** _OPENMP */
    for (std::size_t ii = 0; ii < merged_face_to_tags.size(); ++ii) {
      // Note: Must advance iterator manually each time due to current OMP
      // limitations on some platforms
      auto it = merged_face_to_tags.begin();
      std::advance(it, ii);
      auto const &[face, tags] = *it;
      for (auto const &tag_time : tags) {
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
  Locator _locator; /**< To locate positions in mesh. */

  using TagContainer = std::vector<std::size_t>;
  using FaceToTagContainer = std::unordered_map<Foam::label, TagContainer>;

  std::vector<FaceToTagContainer>
      _face_to_tags; /**< Associative containers of particle tags and
                              times of hitting faces, one per thread. */

  /**
   * @brief Merge collision data across threads and clear separate data.
   *
   * @note Keep only one collision per particle for parallelization reasons.
   *       This converges to the correct behavior in in the limit of small time
   *       step.
   */
  auto merge_and_clear_hits() {
    std::set<std::size_t> tags_collided;
    FaceToTagContainer merged_face_to_tags;

    for (auto &container : _face_to_tags) {
      for (auto &[face, tags] : container) {
        // Find collisions at face, or insert it with no collisions if no
        // collisions have been recorded yet (note that it may remain empty)
        auto [it_new_face, inserted_new_face] =
            merged_face_to_tags.insert({face, {}});
        for (auto tag : tags) {
          // Check if this is the first collision identified for this particle,
          // and if so record it
          auto [it_new_tag, inserted_new_tag] = tags_collided.insert(tag);
          if (inserted_new_tag) {
            it_new_face->second.push_back(tag);
          }
        }
      }
      container.clear();
    }

    return merged_face_to_tags;
  }

  /**
   * @brief Update particle masses and surface concentration at face.
   *
   * @param state State of particle to react.
   *
   * @param state_old Previous state of particle.
   *
   * @param face Mesh face index where reaction occurs.
   */
  template <typename State>
  void update_mass_and_surface_concentration(State state,
                                             State const &state_old,
                                             Foam::label face) {
    double exposure_time = state.time - state_old.time;
    double probability_first_order =
        rate(face) * std::sqrt(cnst::pi * exposure_time / diff_coeff);
    double probability_second_order =
        probability_first_order / (1. + probability_first_order / 2.);
    double delta_mass =
        state_old.mass * (std::exp(-probability_second_order) - 1.);
    double mass_solid =
        surface_concentrations[face] * face_area(face, _locator.mesh());
    if (mass_solid <= rate_constant_ratio_solid_to_fluid * delta_mass) {
      state.mass -= mass_solid / rate_constant_ratio_solid_to_fluid;
      surface_concentrations.erase(face);
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
    double, double, double,
    typename SurfaceReaction_AFluidPlusASolidtoNothing<
        Locator, ParallelOption>::SurfaceConcentrations,
    Locator &&, ParallelOption)
    -> SurfaceReaction_AFluidPlusASolidtoNothing<Locator, ParallelOption>;

/**
   "PTOF/SurfaceReaction.h" \brief Surface reaction that does nothing.
*/
class SurfaceReaction_DoNothing {
public:
  /**
   * @param state State of particle to react (unused).
   *
   * @param state_old Previous state of particle (unused).
   *
   * @param face Mesh face index where reaction occurs (unused).
   *
   * @brief React over exposure time (do nothing).
   */
  template <typename State>
  void operator()(State &state, State const &state_old,
                  Foam::label face) const {}

  /**
   * @param face Mesh face index (unused).
   *
   * @return Surface reaction rate at face (always 0.).
   */
  double rate(Foam::label face) const { return 0.; }

  /**
   * @brief Output generic information about object.
   *
   * @param output Output stream.
   */
  inline static std::ostream &info(std::ostream &output) {
    return output << io::line() << "Surface reaction\n"
                  << io::line() << "None\n"
                  << io::line();
  }
};
} // namespace ptof

#endif /* PTOF_SURFACEREACTION_H */
