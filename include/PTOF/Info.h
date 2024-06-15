/**
\file PTOF/Info.h
\author Tomás Aquino
\date 03/09/2022
*/

#ifndef PTOF_INFO_H
#define PTOF_INFO_H

#include "General/Useful.h"
#include <cstddef>
#include <fieldTypes.H>
#include <point.H>
#include <string>

namespace ptof {
/** \brief Store information about absorption. */
template <typename State> void store_info_absorbed(State &state) {
  state.info.absorbed = 1;
}

/** \brief Store information about boundary condition type. */
template <typename State, typename BC>
void store_info_type(State &state, BC const &type) {
  state.info.type = type;
}

/** \brief Store information about contact point. */
template <typename State>
void store_info_contact(State &state, Foam::point const &contact_point) {
  state.info.set_position(contact_point);
}

/** \brief Store information about cell face.*/
template <typename State> void store_info_face(State &state, Foam::label face) {
  state.info.face = face;
}

/** \brief Store information about time.*/
template <typename State>
void store_info_time(State &state, Foam::scalar time) {
  state.info.time = time;
}

/** \brief Store information about number of reinjections.*/
template <typename State> void store_info_reinjections(State &state) {
  ++state.info.reinjections;
}

/** \struct Info_type PTOF/Info.h "PTOF/Info.h"
 * \brief Information about type.*/
struct Info_Type {
  std::string type{};
};

/** \struct Info_Absorbed PTOF/Info.h "PTOF/Info.h"
 * \brief Information about absorption.*/
struct Info_Absorbed {
  bool absorbed{0};
};

/** \struct Info_Absorbed_Reinjections PTOF/Info.h "PTOF/Info.h"
 * \brief Information about number of reinjections.*/
struct Info_Absorbed_Reinjections {
  bool absorbed{0};
  std::size_t reinjections{0};
};
} // namespace ptof

#endif /* PTOF_INFO_H */
