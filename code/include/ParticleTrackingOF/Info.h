//
//  BoundaryInfo.h
//
//  Created by Tomás Aquino on 09/03/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef BoundaryInfo_h
#define BoundaryInfo_h

#include <cstddef>
#include <string>
#include <fieldTypes.H>
#include "general/useful.h"

namespace ptof
{
  // Store information about absorption
  template <typename State>
  void store_info_absorbed
  (State& state)
  {
    state.info.absorbed = 1;
  }

  // Store information about boundary condition type
  template <typename State, typename BC>
  void store_info_type
  (State& state, BC const& type)
  {
    state.info.type = type;
  }

  // Store information about contact point
  template <typename State>
  void store_info_contact
  (State& state, Foam::point const& contact_point)
  {
    state.info.set_position(contact_point);
  }

  // Store information about cell face
  template <typename State>
  void store_info_face
  (State& state, Foam::label face)
  {
    state.info.face = face;
  }

  // Store information about time
  template <typename State>
  void store_info_time
  (State& state, Foam::scalar time)
  {
    state.info.time = time;
  }

  // Store information about number of reinjections
  template <typename State>
  void store_info_reinjections
  (State& state)
  {
    ++state.info.reinjections;
  }

  // Information about type
  struct Info_Type
  {
    std::string type{};
  };

  // Information about absorption
  struct Info_Absorbed
  {
    bool absorbed{ 0 };
  };

  // Information about number of reinjections
  struct Info_Absorbed_Reinjections
  {
    bool absorbed{ 0 };
    std::size_t reinjections{ 0 };
  };
}

#endif /* BoundaryInfo_h */
