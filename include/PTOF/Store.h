/**
 \file PTOF/Store.h
 \author Tomás Aquino
 \date 09/03/2022
*/

#ifndef PTOF_STORE_H
#define PTOF_STORE_H

#include "General/Meta.h"
#include "General/Useful.h"
#include "PTOF/Boundary.h"

namespace ptof
{
  /** \struct Store_Type PTOF/Store.h "PTOF/Store.h"
   * \brief  Store info about last boundary type. */
  struct Store_Type
  {
    /**
     \brief Store info.
     \note The boundary condition type is selected at compile time through the type \c meta::Selector<BoundaryCondition::Type,type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
     */
    template
    <typename State,
    typename Intersection,
    BoundaryCondition::Type type>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     BoundaryCondition const& implemented,
     meta::Selector<
     BoundaryCondition::Type,
     type>) const
    {
      store_info_type(state, implemented.name(type));
    }
  };

  /** \struct Store_Nothing PTOF/Store.h "PTOF/Store.h"
   * \brief Do not store any info. */
  struct Store_Nothing
  {
    /**
     \brief Store info.
     \note The boundary condition type is selected at compile time through the type \c meta::Selector<BoundaryCondition::Type,type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
    */
    template
    <typename State,
    typename Intersection,
    BoundaryCondition::Type type>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     BoundaryCondition const& implemented,
     meta::Selector<
      BoundaryCondition::Type,
      type>) const
    {}
  };

  /** \struct Store_Absorbed PTOF/Store.h "PTOF/Store.h"
   *  \brief Store info about absorption. */
  struct Store_Absorbed
  {
    /**
     \brief Store info.
     \note The boundary condition type is selected at compile time through the type \c meta::Selector<BoundaryCondition::Type,type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
    */
    template
    <typename State,
    typename Intersection,
    BoundaryCondition::Type type>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     BoundaryCondition const& implemented,
     meta::Selector<
      BoundaryCondition::Type,
      type>) const
    {}
    
    /**
     \brief Store info
     \note The boundary condition type is selected at compile time through the type \c meta::Selector<BoundaryCondition::Type,type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
    */
    template
    <typename State,
    typename Intersection>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     BoundaryCondition const& implemented,
     meta::Selector<
      BoundaryCondition::Type,
      BoundaryCondition::Type::absorbing>) const
    {
      store_info_absorbed(state);
    }
  };

  /** \struct Store_Absorbed_Reinjections PTOF/Store.h "PTOF/Store.h"
    \brief Store info about number absorption and number of reinjections. */
  struct Store_Absorbed_Reinjections
  {
    /**
     \brief Store info.
     \note The boundary condition type is selected at compile time through the type \c meta::Selector<BoundaryCondition::Type,type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
    */
    template
    <typename State,
    typename Intersection,
    BoundaryCondition::Type type>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     BoundaryCondition const& implemented,
     meta::Selector<
      BoundaryCondition::Type,
      type>) const
    {}
    
    /**
     \brief Store info.
     \note The boundary condition type is selected at compile time through the type \c meta::Selector<BoundaryCondition::Type,type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
    */
    template
    <typename State,
    typename Intersection>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     BoundaryCondition const& implemented,
     meta::Selector<
      BoundaryCondition::Type,
      BoundaryCondition::Type::absorbing>) const
    {
      store_info_absorbed(state);
    }

    /**
     \brief Store info.
     \note The boundary condition type is selected at compile time through the type \c meta::Selector<BoundaryCondition::Type,type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
    */
    template
    <typename State,
    typename Intersection>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     BoundaryCondition const& implemented,
     meta::Selector<
      BoundaryCondition::Type,
      BoundaryCondition::Type::custom>) const
    {
      store_info_reinjections(state);
    }
  };

  /** \struct Store_Info_Contact PTOF/Store.h "PTOF/Store.h"
   * \brief Store info about last boundary contact. */
  struct Store_Info_Contact
  {
    /**
     \brief Store info.
     \note The boundary condition type is selected at compile time through the type \c meta::Selector<BoundaryCondition::Type,type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
    */
    template
    <typename State,
    typename Intersection,
    BoundaryCondition::Type type>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     BoundaryCondition const& implemented,
     meta::Selector<
      BoundaryCondition::Type,
      type>) const
    {}
    
    /**
     \brief Store info.
     \note The boundary condition type is selected at compile time through the type \c meta::Selector<BoundaryCondition::Type,type>.
     \param state Current particle state.
     \param state_old Previous particle state.
     \param intersection Information about intersection with boundary.
     \param implemented Information about implemented boundary conditions.
    */
    template
    <typename State,
    typename Intersection>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     BoundaryCondition const& implemented,
     meta::Selector<
      BoundaryCondition::Type,
      BoundaryCondition::Type::info>) const
    {
      store_info_type(state, "info");
      store_info_contact(state, intersection.rawPoint());
      store_info_time(state, state.time);
    }
  };
}

#endif /* PTOF_STORE_H */
