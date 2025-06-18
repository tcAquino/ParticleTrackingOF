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
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  BoundaryConditionList const &implemented,
                  meta::Selector<BoundaryConditionList::Type, type>) const {}
};

/**
   \struct BoundaryInfo_type PTOF/BoundaryInfo.h "PTOF/BoundaryInfo.h"
   \brief Store info about latest boundary type.
*/
struct BoundaryInfo_type {
  /**
     \brief Store boundary type.
     \note The boundary condition type is selected at compile time through the
     type of meta::Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  BoundaryConditionList const &implemented,
                  meta::Selector<BoundaryConditionList::Type, type>) const {
    state.info.type = implemented.name(type);
  }
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
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  BoundaryConditionList const &implemented,
                  meta::Selector<BoundaryConditionList::Type, type>) const {
    state.info.face = intersection.index();
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
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  BoundaryConditionList const &implemented,
                  meta::Selector<BoundaryConditionList::Type, type>) const {
    state.info.contact_point = intersection.point();
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
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  BoundaryConditionList const &implemented,
                  meta::Selector<BoundaryConditionList::Type, type>) const {}

  /**
     \brief Store info about number of reinjections for custom boundaries.
     \note The boundary condition type is selected at compile time through the
     type of Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
  */
  template <typename State, typename Intersection>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  BoundaryConditionList const &implemented,
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
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
  */
  template <typename State, typename Intersection,
            BoundaryConditionList::Type type>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  BoundaryConditionList const &implemented,
                  meta::Selector<BoundaryConditionList::Type, type>) const {
    state.info.face = intersection.index();
  }

  /**
     \brief Store info about number of reinjections for custom boundaries.
     \note The boundary condition type is selected at compile time through the
     type of Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
  */
  template <typename State, typename Intersection>
  void operator()(State &state, State const &state_old,
                  Intersection const &intersection,
                  BoundaryConditionList const &implemented,
                  meta::Selector<BoundaryConditionList::Type,
                                 BoundaryConditionList::Type::custom>) const {
    state.info.face = intersection.index();
    state.info.reinjections++;
  }
};
} // namespace ptof

#endif /* PTOF_BOUNDARYINFO_H */
