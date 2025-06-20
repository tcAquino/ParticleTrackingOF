/**
   \file PTOF/BoundaryInfo.h
   \author Tomás Aquino
   \date 09/03/2022
   \brief Objects to handle storing information upon hitting a boundary.
*/

#ifndef PTOF_BOUNDARYINFO_H
#define PTOF_BOUNDARYINFO_H

#include "General/Meta.h"
#include "General/Useful.h"
#include "PTOF/BoundaryConditionList.h"
#include "PTOF/Info.h"
#include "PTOF/Meta.h"

namespace ptof {
/**
   \struct BoundaryInfo PTOF/BoundaryInfo.h "PTOF/BoundaryInfo.h"
   \brief Do not store any info.
*/
struct BoundaryInfo_Nothing {
  /**
     \brief Store nothing.
     \note The boundary condition type is selected at compile time through the
     type of meta::Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param intersection Information about intersection with boundary.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  meta::Selector<BoundaryConditionList::Type, type>) const {}
};

/**
   \struct BoundaryInfo_face PTOF/BoundaryInfo.h "PTOF/BoundaryInfo.h"
   \brief Store info about latest boundary face.
*/
struct BoundaryInfo_face {
  /**
     \brief Store boundary type.
     \note The boundary condition type is selected at compile time through the
     type of meta::Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param intersection Information about intersection with boundary.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  meta::Selector<BoundaryConditionList::Type, type>) const {
    state.info.boundary_face = intersection.index();
  }
};

/**
   \struct BoundaryInfo_face PTOF/BoundaryInfo.h "PTOF/BoundaryInfo.h"
   \brief Store info about latest contact point.
*/
struct BoundaryInfo_contact_point {
  /**
     \brief Store boundary type.
     \note The boundary condition type is selected at compile time through the
     type of meta::Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param intersection Information about intersection with boundary.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  meta::Selector<BoundaryConditionList::Type, type>) const {
    state.info.contact_point = State::make_position(intersection.point());
  }
};

/**
   \struct BoundaryInfo_face PTOF/BoundaryInfo.h "PTOF/BoundaryInfo.h"
   \brief Store info about number of reinjections.
*/
struct BoundaryInfo_reinjections {
  /**
     \brief Store nothing (unless overloads exist for boundary type).
     \note The boundary condition type is selected at compile time through the
     type of meta::Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param intersection Information about intersection with boundary.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  meta::Selector<BoundaryConditionList::Type, type>) const {}

  /**
     \brief Store info about number of reinjections for custom boundaries.
     \note The boundary condition type is selected at compile time through the
     type of Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param intersection Information about intersection with boundary.
  */
  template <typename State, typename Intersection>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  meta::Selector<BoundaryConditionList::Type,
                                 BoundaryConditionList::Type::custom>) const {
    state.info.reinjections++;
  }
};

/**
   \struct BoundaryInfo_face PTOF/BoundaryInfo.h "PTOF/BoundaryInfo.h"
   \brief Store info about face number of reinjections.
*/
struct BoundaryInfo_face_reinjections {
  /**
     \brief Store nothing (unless overloads exist for boundary type).
     \note The boundary condition type is selected at compile time through the
     type of meta::Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param intersection Information about intersection with boundary.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  meta::Selector<BoundaryConditionList::Type, type>) const {
    generic(state, state_old, intersection);
  }

  /**
     \brief Store info about number of reinjections for custom boundaries.
     \note The boundary condition type is selected at compile time through the
     type of Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param intersection Information about intersection with boundary.
  */
  template <typename State, typename Intersection>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  meta::Selector<BoundaryConditionList::Type,
                                 BoundaryConditionList::Type::custom>) const {
    generic(state, state_old, intersection);
  }

  template <typename State, typename Intersection>
  void generic(State &state, State const &state_old,
               Intersection const &intersection) const {
    state.info.boundary_face = intersection.index();
  }
};

/**
   \struct BoundaryInfo_IfPresent_type_face_contact_point_reinjections
   PTOF/BoundaryInfo.h "PTOF/BoundaryInfo.h"
   \brief Store each of face, contact point, and number of reinjections, if
   present in State::Info.
*/
struct BoundaryInfo_IfPresent_face_contact_point_reinjections {
  /**
     \brief Store each of face and contact point, if present in State::Info.
     \note The boundary condition type is selected at compile time through the
     type of meta::Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param intersection Information about intersection with boundary.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  meta::Selector<BoundaryConditionList::Type, type>) const {
    generic(state, state_old, intersection);
  }

  /**
     \brief Store each of face, contact point, and number of reinjections, if
     present in State::Info.
     \note The boundary condition type is selected at compile time through the
     type of meta::Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param intersection Information about intersection with boundary.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  meta::Selector<BoundaryConditionList::Type,
                                 BoundaryConditionList::Type::custom>) const {
    generic(state, state_old, intersection);
    if constexpr (meta::has_type_v<typename State::contact_point>) {
      state.info.reinjections++;
    }
  }

  template <typename State, typename Intersection>
  void generic(State &state, State const &state_old,
               Intersection const &intersection) const {
    if constexpr (meta::has_boundary_face_v<typename State::Info>) {
      state.info.boundary_face = intersection.index();
    }
    if constexpr (meta::has_contact_point_v<typename State::Info>) {
      state.info.contact_point = State::make_position(intersection.point());
    }
  }
};
} // namespace ptof

#endif /* PTOF_BOUNDARYINFO_H */
