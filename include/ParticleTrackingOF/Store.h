//
//  Store.h
//
//  Created by Tomás Aquino on 09/03/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Store_OF_h
#define Store_OF_h

#include "general/useful.h"
#include "ParticleTrackingOF/Boundary.h"

namespace ptof
{
  // Store info about last boundary type
  struct Store_Type
  {
    template
    <typename State,
    typename Intersection,
    ImplementedBoundaryConditions::Type type>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     ImplementedBoundaryConditions const& implemented,
     useful::Selector<
     ImplementedBoundaryConditions::Type,
     type>) const
    {
      store_info_type(state, implemented.name(type));
    }
  };

  // Do not store any info
  struct Store_Nothing
  {
    template
    <typename State,
    typename Intersection,
    ImplementedBoundaryConditions::Type type>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     ImplementedBoundaryConditions const& implemented,
     useful::Selector<
      ImplementedBoundaryConditions::Type,
      type>) const
    {}
  };

  // Store info about absorption
  struct Store_Absorbed
  {
    template
    <typename State,
    typename Intersection,
    ImplementedBoundaryConditions::Type type>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     ImplementedBoundaryConditions const& implemented,
     useful::Selector<
      ImplementedBoundaryConditions::Type,
      type>) const
    {}
    
    template
    <typename State,
    typename Intersection>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     ImplementedBoundaryConditions const& implemented,
     useful::Selector<
      ImplementedBoundaryConditions::Type,
      ImplementedBoundaryConditions::Type::absorbing>) const
    {
      store_info_absorbed(state);
    }
  };

  // Store info about number absorption
  // and number of reinjections
  struct Store_Absorbed_Reinjections
  {
    template
    <typename State,
    typename Intersection,
    ImplementedBoundaryConditions::Type type>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     ImplementedBoundaryConditions const& implemented,
     useful::Selector<
      ImplementedBoundaryConditions::Type,
      type>) const
    {}
    
    template
    <typename State,
    typename Intersection>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     ImplementedBoundaryConditions const& implemented,
     useful::Selector<
      ImplementedBoundaryConditions::Type,
      ImplementedBoundaryConditions::Type::absorbing>) const
    {
      store_info_absorbed(state);
    }

    template
    <typename State,
    typename Intersection>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     ImplementedBoundaryConditions const& implemented,
     useful::Selector<
      ImplementedBoundaryConditions::Type,
      ImplementedBoundaryConditions::Type::custom>) const
    {
      store_info_reinjections(state);
    }
  };

  // Store info about last boundary contact
  struct Store_Info_Contact
  {
    template
    <typename State,
    typename Intersection,
    ImplementedBoundaryConditions::Type type>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     ImplementedBoundaryConditions const& implemented,
     useful::Selector<
      ImplementedBoundaryConditions::Type,
      type>) const
    {}
    
    template
    <typename State,
    typename Intersection>
    void operator()
    (State& state,
     State const& state_old,
     Intersection const& intersection,
     ImplementedBoundaryConditions const& implemented,
     useful::Selector<
      ImplementedBoundaryConditions::Type,
      ImplementedBoundaryConditions::Type::info>) const
    {
      store_info_type(state, "info");
      store_info_contact(state, intersection.rawPoint());
      store_info_time(state, state.time);
    }
  };
}

#endif /* Store_OF_h */
