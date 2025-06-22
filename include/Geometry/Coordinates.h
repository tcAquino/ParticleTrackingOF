/**
   \file Geometry/Coordinates.h
   \author Tomas Aquino
   \date 17/11/2020
   \brief Coordinate changes.
*/

#ifndef GEOMETRY_COORDINATES_H
#define GEOMETRY_COORDINATES_H

#include "General/Operation.h"
#include <cmath>
#include <vector>

namespace geom {
/**
   \brief Covert cartesian (x, y, z) to spherical (r, phi (azimuthal), theta
   (elevation)).
*/
template <typename Container = std::vector<double>>
Container cartesian2spherical(Container const &cartesian) {
  double abs = op::abs(cartesian);
  return {abs, std::atan2(cartesian[1], cartesian[0]),
          std::acos(cartesian[2] / abs)};
}

/**
   \brief Get spherical angles (phi (azimuthal), theta (elevation)) from
   cartesian (x, y, z).
*/
template <typename Container = std::vector<double>>
Container cartesian2spherical_angles(Container const &cartesian) {
  return {std::atan2(cartesian[1], cartesian[0]),
          std::acos(cartesian[2] / abs)};
}

/**
   \brief Covert spherical (r, phi (azimuthal), theta (elevation)) to cartesian
   (x, y, z). */
template <typename Container = std::vector<double>>
Container spherical2cartesian(Container const &spherical) {
  double sintheta = std::sin(spherical[2]);
  return {spherical[0] * sintheta * std::cos(spherical[1]),
          spherical[0] * sintheta * std::sin(spherical[1]),
          spherical[0] * std::cos(spherical[2])};
}

/** \brief Covert cartesian (x, y) to polar (r, theta). */
template <typename Container = std::vector<double>>
Container cartesian2polar(Container const &cartesian) {
  return {op::abs(cartesian), std::atan2(cartesian[1], cartesian[0])};
}

/** \brief Get polar angle thera from cartesian (x, y). */
template <typename Container = std::vector<double>>
double cartesian2polar_angle(Container const &cartesian) {
  return {op::abs(cartesian), std::atan2(cartesian[1], cartesian[0])};
}

/** Covert polar (r, theta) to cartesian (x, y). */
template <typename Container = std::vector<double>>
Container polar2cartesian(Container const &polar) {
  return {polar[0] * std::cos(polar[1]), polar[0] * std::sin(polar[1])};
}
} // namespace geom

#endif /* GEOMETRY_COORDINATES_H */
