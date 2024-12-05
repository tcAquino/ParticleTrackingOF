/**
   \file PTOF/Store.h
   \author Tomás Aquino
   \date 09/03/2022
   \brief Objects to handle storing information upon hitting a boundary.
*/

#ifndef PTOF_STORE_H
#define PTOF_STORE_H

#include "General/Meta.h"
#include "General/Useful.h"
#include "PTOF/BoundaryConditionList.h"
#include "PTOF/Info.h"

namespace ptof {
/**
   \struct Store_Type PTOF/Store.h "PTOF/Store.h"
   \brief Store info about latest boundary type.
*/
struct Store_Type {
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
    store_info_type(state, implemented.name(type));
  }
};

/**
   \struct Store_Nothing PTOF/Store.h "PTOF/Store.h"
   \brief Do not store any info.
*/
struct Store_Nothing {
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
   \struct Store_Absorbed PTOF/Store.h "PTOF/Store.h"
   \brief Store info about absorption.
*/
struct Store_Absorbed {
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
     \brief Store info about absorption for absorbing boundaries.
     \note The boundary condition type is selected at compile time through the
     type of meta::Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
  */
  template <typename State, typename Intersection>
  void
  operator()(State &state, State const &state_old,
             Intersection const &intersection,
             BoundaryConditionList const &implemented,
             meta::Selector<BoundaryConditionList::Type,
                            BoundaryConditionList::Type::absorbing>) const {
    store_info_absorbed(state);
  }
};

/**
   \struct Store_Absorbed_Reinjections PTOF/Store.h "PTOF/Store.h"
   \brief Store info about absorption and number of reinjections.
*/
struct Store_Absorbed_Reinjections {
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
     \brief Store info about absorption for absorbing boundaries.
     \note The boundary condition type is selected at compile time through the
     type of Selector<BoundaryConditionList::Type, type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
  */
  template <typename State, typename Intersection>
  void
  operator()(State &state, State const &state_old,
             Intersection const &intersection,
             BoundaryConditionList const &implemented,
             meta::Selector<BoundaryConditionList::Type,
                            BoundaryConditionList::Type::absorbing>) const {
    store_info_absorbed(state);
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
    store_info_reinjections(state);
  }
};
} // namespace ptof

#endif /* PTOF_STORE_H */
