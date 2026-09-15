/**
 * @file   BoundaryInfo.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Wed Mar  9 00:00:00 2022
 *
 * @brief Objects to handle storing information upon hitting a boundary.
 */

#ifndef PTOF_BOUNDARYINFO_H
#define PTOF_BOUNDARYINFO_H

#include "General/Meta.h"
#include "General/Parallel.h"
#include "General/Useful.h"
#include "PTOF/Boundary.h"
#include "PTOF/BoundaryConditionList.h"
#include "PTOF/Meta.h"
#include <pointIndexHit.H>

namespace ptof {
/**
 * @brief Base BoundaryInfo object to use polymorphism to store in container.
 */
template <typename State_t, typename Intersection_t = Foam::pointIndexHit>
struct BoundaryInfo_Base {
  using State = State_t;
  using Intersection = Intersection_t;

  virtual ~BoundaryInfo_Base(){};

  /**
   * @brief Generic operation to be applied for absorbing boundaries.
   *
   * @param state Current particle state.
   *
   * @param state_old Previous particle state.
   *
   * @param intersection Information about intersection with boundary.
   *
   * @note The boundary condition type is selected at compile time through the
   *       type of <tt>meta::Selector<BoundaryConditionList::Type, type></tt>.
   */
  virtual void
  operator()(State &state, State const &state_old,
             Intersection const &intersection,
             meta::Selector<BoundaryConditionList::Type,
                            BoundaryConditionList::Type::absorbing>) {
    generic(state, state_old, intersection);
  }

  /**
   * @brief Generic operation to be applied for custom boundaries.
   *
   * @param state Current particle state.
   *
   * @param state_old Previous particle state.
   *
   * @param intersection Information about intersection with boundary.
   *
   * @note The boundary condition type is selected at compile time through the
   *       type of <tt>meta::Selector<BoundaryConditionList::Type, type></tt>.
   */
  virtual void operator()(State &state, State const &state_old,
                          Intersection const &intersection,
                          meta::Selector<BoundaryConditionList::Type,
                                         BoundaryConditionList::Type::custom>) {
    generic(state, state_old, intersection);
  }

  /**
   * @brief Generic operation to be applied for empty boundaries.
   *
   * @param state Current particle state.
   *
   * @param state_old Previous particle state.
   *
   * @param intersection Information about intersection with boundary.
   *
   * @note The boundary condition type is selected at compile time through the
   *       type of <tt>meta::Selector<BoundaryConditionList::Type, type></tt>.
   */
  virtual void operator()(State &state, State const &state_old,
                          Intersection const &intersection,
                          meta::Selector<BoundaryConditionList::Type,
                                         BoundaryConditionList::Type::empty>) {
    generic(state, state_old, intersection);
  }

  /**
   * @brief Generic operation to be applied for inlet boundaries.
   *
   * @param state Current particle state.
   *
   * @param state_old Previous particle state.
   *
   * @param intersection Information about intersection with boundary.
   *
   * @note The boundary condition type is selected at compile time through the
   *       type of <tt>meta::Selector<BoundaryConditionList::Type, type></tt>.
   */
  virtual void operator()(State &state, State const &state_old,
                          Intersection const &intersection,
                          meta::Selector<BoundaryConditionList::Type,
                                         BoundaryConditionList::Type::inlet>) {
    generic(state, state_old, intersection);
  }

  /**
   * @brief Generic operation to be applied for periodic boundaries.
   *
   * @param state Current particle state.
   *
   * @param state_old Previous particle state.
   *
   * @param intersection Information about intersection with boundary.
   *
   * @note The boundary condition type is selected at compile time through the
   *       type of <tt>meta::Selector<BoundaryConditionList::Type, type></tt>.
   */
  virtual void
  operator()(State &state, State const &state_old,
             Intersection const &intersection,
             meta::Selector<BoundaryConditionList::Type,
                            BoundaryConditionList::Type::periodic>) {
    generic(state, state_old, intersection);
  }

  /**
   * @brief Generic operation to be applied for reacting boundaries.
   *
   * @param state Current particle state.
   *
   * @param state_old Previous particle state.
   *
   * @param intersection Information about intersection with boundary.
   *
   * @note The boundary condition type is selected at compile time through the
   *       type of <tt>meta::Selector<BoundaryConditionList::Type, type></tt>.
   */
  virtual void
  operator()(State &state, State const &state_old,
             Intersection const &intersection,
             meta::Selector<BoundaryConditionList::Type,
                            BoundaryConditionList::Type::reacting>) {
    generic(state, state_old, intersection);
  }

  /**
   * @brief Generic operation to be applied for reflecting boundaries.
   *
   * @param state Current particle state.
   *
   * @param state_old Previous particle state.
   *
   * @param intersection Information about intersection with boundary.
   *
   * @note The boundary condition type is selected at compile time through the
   *       type of <tt>meta::Selector<BoundaryConditionList::Type, type></tt>.
   */
  virtual void
  operator()(State &state, State const &state_old,
             Intersection const &intersection,
             meta::Selector<BoundaryConditionList::Type,
                            BoundaryConditionList::Type::reflecting>) {
    generic(state, state_old, intersection);
  }

  /**
   * @brief Generic operation to be applied for all boundary types.
   *
   * @param state Current particle state.
   *
   * @param state_old Previous particle state.
   *
   * @param intersection Information about intersection with boundary.
   *
   * @note The boundary condition type is selected at compile time through the
   *       type of <tt>meta::Selector<BoundaryConditionList::Type, type></tt>.
   */
  virtual void generic(State &state, State const &state_old,
                       Intersection const &intersection) = 0;
};

/** @brief Do not store any info. */
template <typename State, typename Intersection_t = Foam::pointIndexHit>
struct BoundaryInfo_Nothing final
    : public BoundaryInfo_Base<State, Intersection_t> {
  using Intersection =
      typename BoundaryInfo_Base<State, Intersection_t>::Intersection;

  void generic(State &state, State const &state_old,
               Intersection const &intersection) override {}
};

/** @brief Store info about latest boundary face. */
template <typename State, typename Intersection_t = Foam::pointIndexHit>
struct BoundaryInfo_boundary_face final
    : public BoundaryInfo_Base<State, Intersection_t> {
  using Intersection =
      typename BoundaryInfo_Base<State, Intersection_t>::Intersection;

  void generic(State &state, State const &state_old,
               Intersection const &intersection) override {
    state.info.boundary_face = intersection.index();
  }
};

/** @brief Store info about latest contact point. */
template <typename State, typename Intersection_t = Foam::pointIndexHit>
struct BoundaryInfo_contact_point final
    : public BoundaryInfo_Base<State, Intersection_t> {
  using Intersection =
      typename BoundaryInfo_Base<State, Intersection_t>::Intersection;

  void generic(State &state, State const &state_old,
               Intersection const &intersection) override {
    state.info.contact_point = State::make_position(intersection.point());
  }
};

/** @brief Store info about number of reinjections. */
template <typename State, typename Intersection_t = Foam::pointIndexHit>
struct BoundaryInfo_reinjections final
    : public BoundaryInfo_Base<State, Intersection_t> {
  using Intersection =
      typename BoundaryInfo_Base<State, Intersection_t>::Intersection;

  void generic(State &state, State const &state_old,
               Intersection const &intersection) override {}

  void
  operator()(State &state, State const &state_old,
             Intersection const &intersection,
             meta::Selector<BoundaryConditionList::Type,
                            BoundaryConditionList::Type::custom>) override {
    state.info.reinjections++;
  }
};

/** @brief Store info about face and number of reinjections. */
template <typename State, typename Intersection_t = Foam::pointIndexHit>
struct BoundaryInfo_boundary_face_reinjections final
    : public BoundaryInfo_Base<State, Intersection_t> {
  using Intersection =
      typename BoundaryInfo_Base<State, Intersection_t>::Intersection;

  void generic(State &state, State const &state_old,
               Intersection const &intersection) override {
    state.info.boundary_face = intersection.index();
  }

  void
  operator()(State &state, State const &state_old,
             Intersection const &intersection,
             meta::Selector<BoundaryConditionList::Type,
                            BoundaryConditionList::Type::custom>) override {
    generic(state, state_old, intersection);
    state.info.reinjections++;
  }
};

/**
 * @brief Store each of face, contact point, and number of reinjections, if
 *        present in State::Info.
 */
template <typename State, typename Intersection_t = Foam::pointIndexHit>
struct BoundaryInfo_IfPresent_boundary_face_contact_point_reinjections final
    : public BoundaryInfo_Base<State, Intersection_t> {
  using Intersection =
      typename BoundaryInfo_Base<State, Intersection_t>::Intersection;

  void generic(State &state, State const &state_old,
               Intersection const &intersection) override {
    if constexpr (meta::has_boundary_face_v<typename State::Info>) {
      state.info.boundary_face = intersection.index();
    }
    if constexpr (meta::has_contact_point_v<typename State::Info>) {
      state.info.contact_point = State::make_position(intersection.point());
    }
  }

  void
  operator()(State &state, State const &state_old,
             Intersection const &intersection,
             meta::Selector<BoundaryConditionList::Type,
                            BoundaryConditionList::Type::custom>) override {
    generic(state, state_old, intersection);
    if constexpr (meta::has_reinjections_v<typename State::Info>) {
      state.info.reinjections++;
    }
  }
};

/** @brief Record info about net mass consumed at reactive boundaries. */
template <typename State,
          typename ParallelOption = typename par::ParallelOptions::Serial,
          typename Intersection_t = Foam::pointIndexHit>
struct BoundaryInfo_Record_mass_reacted_face final
    : public BoundaryInfo_Base<State, Intersection_t> {
  using Intersection =
      typename BoundaryInfo_Base<State, Intersection_t>::Intersection;
  void generic(State &state, State const &state_old,
               Intersection const &intersection) override {}

  void
  operator()(State &state, State const &state_old,
             Intersection const &intersection,
             meta::Selector<BoundaryConditionList::Type,
                            BoundaryConditionList::Type::reacting>) override {
    double mass_change = state.mass - state_old.mass;
    auto [it, inserted] = _masses[par::get_thread_num(ParallelOption{})].insert(
        {intersection.index(), mass_change});
    if (!inserted) {
      it->second += mass_change;
    }
  }

  auto const &masses() {
    for (std::size_t ii = 1; ii < _masses.size(); ++ii) {
      for (auto const &face_to_mass : _masses[ii]) {
        auto [it, inserted] = _masses[0].insert(face_to_mass);
        if (!inserted) {
          it->second += face_to_mass.second;
        }
      }
      // Clear all except first container to make future merges more efficient
      _masses[ii].clear();
    }
    // Return first container with merged data
    return _masses[0];
  }

private:
  using FaceToMass = std::unordered_map<Foam::label, double>;
  std::vector<FaceToMass> _masses;
};

/**
 * @brief Record info about net mass consumed at reactive boundaries, accounting
 *        for periodicity.
 */
template <typename State,
          typename ParallelOption = typename par::ParallelOptions::Serial,
          typename Intersection_t = Foam::pointIndexHit>
struct BoundaryInfo_Record_mass_reacted_face_periodic final
    : public BoundaryInfo_Base<State, Intersection_t> {
  using Intersection =
      typename BoundaryInfo_Base<State, Intersection_t>::Intersection;
  void generic(State &state, State const &state_old,
               Intersection const &intersection) override {}

  void
  operator()(State &state, State const &state_old,
             Intersection const &intersection,
             meta::Selector<BoundaryConditionList::Type,
                            BoundaryConditionList::Type::reacting>) override {
    double mass_change = state.mass - state_old.mass;
    auto [it, inserted] = _masses[par::get_thread_num(ParallelOption{})].insert(
        {{intersection.index(), state.periodicity}, mass_change});
    if (!inserted) {
      it->second += mass_change;
    }
  }

  auto const &masses() {
    for (std::size_t ii = 1; ii < _masses.size(); ++ii) {
      for (auto const &face_periodicity_to_mass : _masses[ii]) {
        auto [it, inserted] = _masses[0].insert(face_periodicity_to_mass);
        if (!inserted) {
          it->second += face_periodicity_to_mass.second;
        }
      }
      // Clear all except first container to make future merges more efficient
      _masses[ii].clear();
    }
    // Return first container with merged data
    return _masses[0];
  }

private:
  using Periodicity = typename State::Periodicity;
  using FacePeriodicityToMass =
      std::unordered_map<std::pair<Foam::label, Periodicity>, double,
                         useful::Hash_pair<Foam::label, Periodicity>>;
  std::vector<FacePeriodicityToMass> _masses{std::vector<FacePeriodicityToMass>(
      par::get_num_threads(ParallelOption{}))};
};
} // namespace ptof

#endif /* PTOF_BOUNDARYINFO_H */
