/**
 \file Geometry/Coordinates.h
 \author Tomás Aquino
 \date 17/11/2020
*/

#ifndef Coordinates_h
#define Coordinates_h

#include <cmath>
#include <vector>
#include "general/Operations.h"

namespace geometry
{
  /** \brief Covert cartesian (x, y, z) to spherical (r, phi (azimuthal), theta (elevation)). */
  template <typename Container = std::vector<double>>
  Container cartesian2spherical(Container const& cartesian)
  {
    double abs = operation::abs(cartesian);
    return {
      abs,
      std::atan2(cartesian[1], cartesian[0]),
      std::acos(cartesian[2]/abs) };
  }

  /** \brief Covert spherical (r, phi (azimuthal), theta (elevation)) to cartesian (x, y, z). */
  template <typename Container = std::vector<double>>
  Container spherical2cartesian(std::vector<double> const& spherical)
  {
    double sintheta = std::sin(spherical[2]);
    return {
      spherical[0]*sintheta*std::cos(spherical[1]),
      spherical[0]*sintheta*std::sin(spherical[1]),
      spherical[0]*std::cos(spherical[2]) };
  }

  /** \brief Covert cartesian (x, y) to polar (r, theta). */
  template <typename Container = std::vector<double>>
  Container cartesian2polar(std::vector<double> const& cartesian)
  {
    return {
      operation::abs(cartesian),
      std::atan2(cartesian[1], cartesian[0]) };
  }

  /** Covert polar (r, theta) to cartesian (x, y). */
  template <typename Container = std::vector<double>>
  Container polar2cartesian(std::vector<double> const& polar)
  {
    return {
      polar[0]*std::cos(polar[1]),
      polar[0]*std::sin(polar[1]) };
  }
}

#endif /* Coordinates_h */
