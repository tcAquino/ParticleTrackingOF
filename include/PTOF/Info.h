/**
 * @file   Info.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Sat Sep  3 00:00:00 2022
 *
 * @brief Objects and utilities to store information in particle states.
 */

#ifndef PTOF_INFO_H
#define PTOF_INFO_H

#include <cstddef>
#include <fieldTypes.H>
#include <point.H>

namespace ptof {
/** @brief Information about absorption. */
struct Info_absorbed {
  bool absorbed;
};

/** @brief Information about absorption and adsorption. */
struct Info_absorbed_adsorbed {
  bool absorbed;
  bool adsorbed;
};

/**
 * @brief Information about absorption and boundary face where absorption
 *        happened.
 */
struct Info_absorbed_boundary_face {
  bool absorbed;
  Foam::label boundary_face;
};

/**
 * @brief Information about absorption and boundary face where absorption
 *        happened.
 */
struct Info_absorbed_adsorbed_boundary_face {
  bool absorbed;
  bool adsorbed;
  Foam::label boundary_face;
};

/** @brief Information about absorption and number of reinjections. */
struct Info_absorbed_reinjections {
  bool absorbed;
  std::size_t reinjections;
};

/** @brief Information about absorption and number of reinjections. */
struct Info_absorbed_boundary_face_reinjections {
  bool absorbed;
  Foam::label boundary_face;
  std::size_t reinjections;
};

/** @brief Information about absorption and number of reinjections. */
struct Info_absorbed_adsorbed_boundary_face_reinjections {
  bool absorbed;
  bool adsorbed;
  Foam::label boundary_face;
  std::size_t reinjections;
};

/**
 * @brief Information about absorption, adsorption, and number of reinjections.
 */
struct Info_absorbed_adsorbed_reinjections {
  bool absorbed;
  bool adsorbed;
  std::size_t reinjections;
};

/** @brief Information about absorption and contact point. */
struct Info_absorbed_contact_point {
  bool absorbed;
  Foam::point contact_point;
};

/** @brief Information about absorption, adsorption, and contact point. */
struct Info_absorbed_adsorbed_contact_point {
  bool absorbed;
  bool adsorbed;
  Foam::point contact_point;
};
} // namespace ptof

#endif /* PTOF_INFO_H */
