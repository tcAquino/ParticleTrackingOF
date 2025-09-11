/**
   \file Geometry/Boundary.h
   \author Tomas Aquino
   \date 08/02/2019

   \brief Boundary-enforcing classes for particle tracking.

   \details
   A Boundary class should implement the following basic functionality:

   \code{.cpp}
   class Example_Boundary
   {
   public:

   // Check if position is out of bounds
   template <typename Position>
   bool out_of_bounds(Position const& position) const
   {
   // Return true if boundary condition is to be applied, false otherwise
   }

   // Enforce boundary condition
   template <typename State>
   bool operator()(State& state, State const& state_old) const
   {
   // Apply the boundary condition given current and previous state
   // Return true if anything was done, false otherwise
   }
   };
   \endcode
*/

#ifndef GEOMETRY_BOUNDARY_H
#define GEOMETRY_BOUNDARY_H

#include "General/Operation.h"
#include "Geometry/SymmetryPlanes.h"
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>

namespace geom {
/**
   \param position 1D position.
   \param boundaries Pair of lower and upper boundary location.
   \return Displacement to be added to \p position to place it within periodic
   boundaries.
*/
double boundary_periodic(double position,
                         std::pair<double, double> const &boundaries) {
  double box_size = boundaries.second - boundaries.first;
  double pos = position - boundaries.first;
  return -std::floor(pos / box_size) * box_size;
}

/**
   \param position 1D position .
   \param boundaries Pair of lower and upper boundary location.
   \return Pair of displacement to be added to \p position and the associated
   signed number of domain displacements out to place it within periodic
   boundaries.
*/
std::pair<double, int> boundary_periodic_with_outside_info(
    double position, std::pair<double, double> const &boundaries) {
  double box_size = boundaries.second - boundaries.first;
  double pos = position - boundaries.first;
  int outside = int(std::floor(pos / box_size));
  return {-outside * box_size, outside};
}

/**
   \param position 1D position.
   \param boundaries Pair of lower and upper boundary location.
   \return Displacement to be added to 1D \p position to place it within
   reflecting boundaries.
*/
double boundary_reflecting(double position,
                           std::pair<double, double> const &boundaries) {
  double box_size = boundaries.second - boundaries.first;
  double pos = position - boundaries.first;
  double nr_boxes_jumped = std::trunc(pos / box_size);
  if (int(nr_boxes_jumped) % 2 == 0) {
    return -pos + std::abs(pos - nr_boxes_jumped * box_size);
  } else {
    return -pos + box_size - std::abs(pos - nr_boxes_jumped * box_size);
  }
}

/**
   \param position 1D position.
   \param boundaries Pair of lower and upper boundary location.
   \return \c true if \p position is out of bounds, \c false otherwise.
*/
bool out_of_bounds_box(double position,
                       std::pair<double, double> const &boundaries) {
  return position < boundaries.first || position > boundaries.second;
}

/**
   \param position position container.
   \param boundaries Container of pairs of lower and upper boundary locations
   along each dimension.
   \return \c true if \p position is out of bounds, \c false otherwise. */
template <typename Position = std::vector<double>,
          typename Boundaries = std::vector<std::pair<double, double>>>
bool out_of_bounds_box(Position const &position, Boundaries const &boundaries) {
  for (std::size_t dd = 0; dd < boundaries.size(); ++dd) {
    if (out_of_bounds_box(position[dd], boundaries[dd])) {
      return true;
    }
  }
  return false;
}

/**
   \class Boundary_DoNothing Geometry/Boundary.h "Geometry/Boundary.h"
   \brief Boundary with no effect.
*/
class Boundary_DoNothing {
public:
  /**
     \brief Check if \p position is out of bounds.
     \return \c false (always in bounds).
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return false;
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC.
     \param state_old Previous particle state (unused).
     \return \c false (no effect).
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) const {
    return false;
  }
};

/**
   \class Boundary_Periodic_1d Geometry/Boundary.h "Geometry/Boundary.h"
   \brief Periodic boundary in 1d./
   \note State must define:
   - \c position [Scalar]
*/
class Boundary_Periodic_1d {
public:
  const std::pair<double, double>
      boundaries; /**< Lower and upper boundary location. */

  /**
     \brief Constructor.
     \param boundaries Pair of lower and upper boundary location.
  */
  Boundary_Periodic_1d(std::pair<double, double> boundaries)
      : boundaries{boundaries} {}

  /**
     \brief Check if position is out of bounds.
     \param position Position to check
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position> bool out_of_bounds(Position position) const {
    return OutOfBounds_Box(position, boundaries);
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC to.
     \param state_old Previous particle state (unused).
     \return \c true if there was an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) const {
    if (!out_of_bounds(state.position)) {
      return false;
    }
    state.position += boundary_periodic(state.position, boundaries);
    return true;
  }

  /**
     \brief Translate a position according to projection in periodic cell.
     \param position Position to translate.
     \param projection Position within periodic cell.
  */
  template <typename Position, typename Projection = double>
  void translate(Position &position, Projection projection) const {
    op::plus_inplace(position,
                     (boundaries.second - boundaries.first) * projection);
  }
};

/**
   \class Boundary_Periodic Geometry/Boundary.h "Geometry/Boundary.h"
   \brief  Periodic boundaries along each dimension.
   \note State must define:
   - \c position [Container]
*/
class Boundary_Periodic {
public:
  const std::vector<std::pair<double, double>>
      boundaries; /**< Lower and upper boundary locations along each dimension.
                   */
  const std::vector<double>
      domain_dimensions; /**< Size of domain along each dimension. */

  /**
     \brief Constructor.
     \param boundaries Container of pairs of lower and upper boundary locations
     along each dimension. */
  Boundary_Periodic(std::vector<std::pair<double, double>> boundaries)
      : boundaries{boundaries}, domain_dimensions{make_dimensions(boundaries)} {
  }

  /**
     \brief Check if position is out of bounds.
     \param position Position to check
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return out_of_bounds_box(position, boundaries);
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC to.
     \param state_old Previous particle state (unused).
     \return \c true if there was an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) const {
    if (!out_of_bounds(state.position)) {
      return false;
    }
    for (std::size_t dd = 0; dd < boundaries.size(); ++dd) {
      state.position[dd] +=
          boundary_periodic(state.position[dd], boundaries[dd]);
    }
    return true;
  }

  /**
     \brief Translate a position according to projection in periodic cell.
     \param position Position to translate.
     \param projections Projection components in periodic cell.
  */
  template <typename Position, typename Projections = std::vector<double>>
  void translate(Position &position, Projections const &projections) const {
    auto increment = op::times(domain_dimensions, projections);
    for (std::size_t dd = 0; dd < increment.size(); ++dd) {
      position[dd] += increment[dd];
    }
  }

private:
  /**
     \param boundaries Container of pairs of lower and upper boundary locations
     along each dimension.
     \return Domain dimensions.
  */
  std::vector<double> make_dimensions(
      std::vector<std::pair<double, double>> const &boundaries) const {
    std::vector<double> dimensions;
    for (auto const &val : boundaries) {
      dimensions.push_back(val.second - val.first);
    }

    return dimensions;
  }
};

/**
   \class Boundary_Periodic_Dim Geometry/Boundary.h "Geometry/Boundary.h"
   \brief Periodic boundary along dimension dd.
   \note State must define:
   - \c position [Container]
*/
template <std::size_t dd> class Boundary_Periodic_Dim {
public:
  const std::pair<double, double>
      boundaries; /**< Lower and upper boundary location along dimension \c dd.
                   */
  static const std::size_t dim =
      dd; /**< Cartesian dimension along which BC holds. */

  /**
     \brief Constructor.
     \param half_width Half-size of boundary, centered at zero.
  */
  Boundary_Periodic_Dim(double half_width)
      : boundaries{-half_width, half_width} {}

  /**
     \brief Constructor.
     \param boundaries Pair of lower and upper boundary location.
  */
  Boundary_Periodic_Dim(std::pair<double, double> boundaries)
      : boundaries{boundaries} {}

  /**
     \brief Check if position is out of bounds.
     \param position Position to check
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return out_of_bounds_box(position[dd], boundaries);
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC to.
     \param state_old Previous particle state (unused).
     \return \c true if there was an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) const {
    if (!out_of_bounds(state.position)) {
      return false;
    }
    state.position[dd] += boundary_periodic(state.position[dd], boundaries);
    return true;
  }
};

/**
   \class Boundary_Periodic_WithOutsideInfo Geometry/Boundary.h
   "Geometry/Boundary.h"
   \brief Periodic boundaries along each dimension with information about where
   position would be outside domain.
   \note State must define:
   - \c position [Container]
   - \c periodicity [Container of signed integers] (counting how many domains
   have been traveled along each dimension)
*/
class Boundary_Periodic_WithOutsideInfo {
public:
  const std::vector<std::pair<double, double>>
      boundaries; /**< Lower and upper boundaries along each dimension. */
  const std::vector<double>
      domain_dimensions; /**< Size of domain along each dimension. */

  /**
     \brief Constructor.
     \param boundaries Container of pairs of lower and upper boundary locations
     along each dimension. */
  Boundary_Periodic_WithOutsideInfo(
      std::vector<std::pair<double, double>> boundaries)
      : boundaries{boundaries}, domain_dimensions{make_dimensions(boundaries)} {
  }

  /**
     \brief Check if position is out of bounds.
     \param position Position to check.
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return out_of_bounds_box(position, boundaries);
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC to.
     \param state_old Previous particle state (unused).
     \return \c true if there was an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) const {
    if (!out_of_bounds(state.position)) {
      return false;
    }
    for (std::size_t dd = 0; dd < boundaries.size(); ++dd) {
      auto change_outside = boundary_periodic_with_outside_info(
          state.position[dd], boundaries[dd]);
      state.position[dd] += change_outside.first;
      state.periodicity[dd] += change_outside.second;
    }
    return true;
  }

  /**
     \brief Translate a position according to projection in periodic cell.
     \param position Position to translate.
     \param projections Projection components in periodic cell.
  */
  template <typename Position, typename Projections = std::vector<double>>
  void translate(Position &position, Projections const &projections) const {
    auto increment = op::times(domain_dimensions, projections);
    for (std::size_t dd = 0; dd < increment.size(); ++dd) {
      position[dd] += increment[dd];
    }
  }

private:
  /**
     \param boundaries Container of pairs of lower and upper boundary locations
     along each dimension.
     \return Domain dimensions.
  */
  std::vector<double> make_dimensions(
      std::vector<std::pair<double, double>> const &boundaries) const {
    std::vector<double> dimensions;
    for (auto const &val : boundaries) {
      dimensions.push_back(val.second - val.first);
    }
    return dimensions;
  }
};

/**
   \class Boundary_Reflecting_1d Geometry/Boundary.h "Geometry/Boundary.h"
   \brief Reflecting boundary in one dimension.
   \note State must define:
   - \c position [Scalar]
*/
class Boundary_Reflecting_1d {
public:
  const std::pair<double, double>
      boundaries; /**< Lower and upper boundaries. */

  /**
     \brief Constructor.
     \param boundaries Pairs of lower and upper boundary location.
  */
  Boundary_Reflecting_1d(std::pair<double, double> boundaries)
      : boundaries{boundaries} {}

  /**
     \brief Check if position is out of bounds.
     \param position Position to check
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return out_of_bounds_box(position, boundaries);
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC to.
     \param state_old Previous particle state (unused).
     \return \c true if there was an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) const {
    if (!out_of_bounds(state.position)) {
      return false;
    }
    state.position += boundary_reflecting(state.position, boundaries);
    return true;
  }
};

/**
   \class Boundary_Reflecting Geometry/Boundary.h "Geometry/Boundary.h"
   \brief Reflecting boundaries along along each dimension.
   \note State must define:
   - \c position [Container]
*/
class Boundary_Reflecting {
public:
  const std::vector<std::pair<double, double>>
      boundaries; /**< Lower and upper boundaries along each dimension. */

  /**
     \brief Constructor.
     \param boundaries Container of pairs of lower and upper boundary locations
     along each dimension.
  */
  Boundary_Reflecting(std::vector<std::pair<double, double>> boundaries)
      : boundaries{boundaries} {}

  /**
     \brief Check if position is out of bounds.
     \param position Position to check
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return out_of_bounds_box(position, boundaries);
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC to.
     \param state_old Previous particle state (unused).
     \return \c true if there was an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) const {
    if (!out_of_bounds(state.position)) {
      return false;
    }
    for (std::size_t dd = 0; dd < boundaries.size(); ++dd) {
      state.position[dd] +=
          boundary_reflecting(state.position[dd], boundaries[dd]);
    }
    return true;
  }
};

/**
   \class Boundary_Reflecting_Dim Geometry/Boundary.h "Geometry/Boundary.h"
   \brief Reflecting boundary along dimension dd.
   \note State must define:
   - \c position [Container]
*/
template <std::size_t dd> class Boundary_Reflecting_Dim {
public:
  const std::pair<double, double>
      boundaries; /**< Lower and upper boundaries along dimension dd. */
  static const std::size_t dim = dd; /**< Dimension for outside visibility.*/

  /** Constructor.
   \param half_width Half-size of boundary, centered at zero. */
  Boundary_Reflecting_Dim(double half_width)
      : boundaries{-half_width, half_width} {}

  /** Constructor.
  \param boundaries Pairs of lower and upper boundary location. */
  Boundary_Reflecting_Dim(std::pair<double, double> boundaries)
      : boundaries{boundaries} {}

  /**
     \brief Check if position is out of bounds.
     \param position Position to check
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return out_of_bounds_box(position[dd], boundaries);
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC to.
     \param state_old Previous particle state (unused).
     \return \c true if there was an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) const {
    if (!out_of_bounds(state.position)) {
      return false;
    }
    state.position[dd] += boundary_reflecting(state.position[dd], boundaries);
    return true;
  }
};

/**
   \class Boundary_RadialReflecting_2d Geometry/Boundary.h "Geometry/Boundary.h"
   \brief Reflecting boundaries on the inside of a circle.
   \note State must define:
   - \c position [Container]
*/
class Boundary_RadialReflecting_2d {
public:
  const double radius;

  /**
     \brief Constructor.
     \param radius Domain radius.
  */
  Boundary_RadialReflecting_2d(double radius) : radius{radius} {}

  /**
     \brief Check if position is out of bounds.
     \param position Position to check
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return op::abs_sq(position) > _radius_sq;
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC to.
     \param state_old Previous particle state.
     \return \c true if there was an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old) const {
    if (!out_of_bounds(state.position)) {
      return false;
    }

    auto position_old = state_old.position;
    std::vector<double> jump = op::minus(state.position, position_old);
    std::vector<double> tangent_at_contact(2);

    // Iterate reflection procedure until position is inside boundary
    while (1) {
      // Fraction alpha of jump up to the boundary
      double pos_dot_jump = op::dot(position_old, jump);
      double norm_sq_jump = op::abs_sq(jump);
      double pos_sq_old = op::abs_sq(position_old);
      double aux = pos_dot_jump / norm_sq_jump;
      double alpha =
          std::sqrt(aux * aux + (_radius_sq - pos_sq_old) / norm_sq_jump) - aux;

      // Avoid numerical issues
      if (alpha <= 0. || std::isnan(alpha)) {
        op::times_scalar_inplace(radius / op::abs(state.position),
                                 state.position);
        break;
      }

      // Place old position at contact
      for (std::size_t dd = 0; dd < 2; ++dd) {
        position_old[dd] += alpha * jump[dd];
      }

      // Tangent to boundary at contact
      tangent_at_contact[0] = position_old[1] / radius;
      tangent_at_contact[1] = -position_old[0] / radius;

      // Angle between jump and tangent to boundary at contact
      double tangent_dot_jump = op::dot(tangent_at_contact, jump);
      double cos_contact_angle = tangent_dot_jump / std::sqrt(norm_sq_jump);
      double contact_angle = std::acos(cos_contact_angle);

      // Outgoing jump at collision, rotate outer jump and flip the direction
      op::rotate(jump, -2. * contact_angle);
      op::times_scalar_inplace(1. - alpha, jump);

      // Place particle at reflected position
      op::plus(position_old, jump, state.position);

      // If inside boundary, done
      if (!out_of_bounds(state.position)) {
        break;
      }
    }

    return true;
  }

private:
  const double _radius_sq{radius * radius}; /**< Domain radius squared. */
};

/**
   \class Boundary_Open_RadialReflecting_3d "Geometry/Boundary.h"
   \brief Reflecting boundaries on the inside of an infinite cylinder in 3d.
   \tparam dd_open Cylinder longitudinal axis, must be 0 or 2.
   \note State must define:
   - \c position [Container]
*/
template <std::size_t dd_open = 0> class Boundary_Open_RadialReflecting_3d {
  static_assert(dd_open == 0 || dd_open == 2,
                "Cylinder axis dimension must be 0 or 2");

public:
  const double radius; /**< Domain radius along reflecting dimensions. */

  /**
     \brief Constructor.
     \param radius Domain radius along reflecting dimensions.
  */
  Boundary_Open_RadialReflecting_3d(double radius) : radius{radius} {}

  /**
     \brief Check if position is out of bounds.
     \param position Position to check
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    double radial_pos_sq =
        std::inner_product(position.cbegin() + _begin_transverse,
                           position.cbegin() + _begin_transverse + 2,
                           position.cbegin() + _begin_transverse, 0.);
    return radial_pos_sq > _radius_sq;
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC to.
     \param state_old Previous particle state.
     \return \c true if there was an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old) const {
    if (!out_of_bounds(state.position)) {
      return false;
    }

    std::vector<double> position_radial_old(2);
    std::vector<double> jump_radial(2);
    std::vector<double> tangent_at_contact(2);

    for (std::size_t dd = 0; dd < 2; ++dd) {
      jump_radial[dd] = state.position[dd + _begin_transverse] -
                        state_old.position[dd + _begin_transverse];
      position_radial_old[dd] = state_old.position[dd + _begin_transverse];
    }

    // Iterate reflection procedure until position is inside boundary
    while (1) {
      // Fraction alpha of jump up to the boundary
      double pos_dot_jump = std::inner_product(position_radial_old.cbegin(),
                                               position_radial_old.cend(),
                                               jump_radial.cbegin(), 0.);
      double norm_sq_jump = std::inner_product(
          jump_radial.cbegin(), jump_radial.cend(), jump_radial.cbegin(), 0.);
      double radial_pos_sq_old = op::abs_sq(position_radial_old);
      double aux = pos_dot_jump / norm_sq_jump;
      double alpha = std::sqrt(aux * aux + (_radius_sq - radial_pos_sq_old) /
                                               norm_sq_jump) -
                     aux;

      // Avoid numerical issues
      if (alpha <= 0. || std::isnan(alpha)) {
        double radial_pos_sq =
            std::inner_product(state.position.cbegin() + _begin_transverse,
                               state.position.cbegin() + _begin_transverse + 2,
                               state.position.cbegin() + _begin_transverse, 0.);
        for (std::size_t dd = 0; dd < 2; ++dd) {
          state.position[dd + _begin_transverse] *=
              radius / std::sqrt(radial_pos_sq);
        }
        break;
      }

      // Place old position at contact
      for (std::size_t dd = 0; dd < 2; ++dd) {
        position_radial_old[dd] += alpha * jump_radial[dd];
      }

      // Tangent to boundary at contact
      tangent_at_contact[0] = position_radial_old[1] / radius;
      tangent_at_contact[1] = -position_radial_old[0] / radius;

      // Angle between jump and tangent to boundary at contact
      double tangent_dot_jump = std::inner_product(tangent_at_contact.cbegin(),
                                                   tangent_at_contact.cend(),
                                                   jump_radial.cbegin(), 0.);
      double cos_contact_angle = tangent_dot_jump / std::sqrt(norm_sq_jump);
      double contact_angle = std::acos(cos_contact_angle);

      // Outgoing jump at collision, rotate outer jump and flip the direction
      op::rotate(jump_radial, -2. * contact_angle);
      op::times_scalar_inplace(1. - alpha, jump_radial);

      // Place particle at reflected position
      for (std::size_t dd = 0; dd < 2; ++dd) {
        state.position[dd + _begin_transverse] =
            position_radial_old[dd] + jump_radial[dd];
      }

      // If inside, done
      if (!out_of_bounds(state.position)) {
        break;
      }
    }

    return true;
  }

private:
  const double _radius_sq{
      radius *
      radius}; /**< Domain radius along reflecting dimensions squared. */
  const std::size_t _begin_transverse{
      dd_open == 0 ? 1 : 0}; /** First reflecting dimension. */
};

/**
   \class Boundary_Periodic_SymmetryPlanes Geometry/Boundary.h
   "Geometry/Boundary.h"
   \brief Periodic boundaries along each symmetry plane.
   \note State must define:
   - \c position [Container]
*/
template <typename SymmetryPlanes> class Boundary_Periodic_SymmetryPlanes {
public:
  const SymmetryPlanes
      symmetry_planes; /**> Object describing domain symmetry planes (see
                          Geometry/SymmetryPlanes.h). */
  const double scale;  /**> Scale factor. */
  const std::vector<double> origin; /**> Coordinate origin. */

  /**
     \brief Constructor.
     \param symmetry_planes Object describing domain symmetry planes (see
     "Geometry/SymmetryPlanes.h").
     \param scale Scale factor.
     \param origin Coordinate origin.
  */
  Boundary_Periodic_SymmetryPlanes(SymmetryPlanes symmetry_planes,
                                   double scale = 1.,
                                   std::vector<double> origin = {})
      : symmetry_planes{symmetry_planes}, scale{scale},
        origin{origin.size() == 0.
                   ? std::vector<double>(symmetry_planes.dim, 0.)
                   : origin} {}

  /**
     \brief Constructor.
     \brief Use default-constructed domain symmetry planes.
     \param scale Scale factor.
     \param origin Coordinate origin.
  */
  Boundary_Periodic_SymmetryPlanes(double scale = 1.,
                                   std::vector<double> origin = {})
      : symmetry_planes{}, scale{scale}, origin{
                                             origin.size() == 0.
                                                 ? std::vector<double>(
                                                       symmetry_planes.dim, 0.)
                                                 : origin} {}

  /**
     \brief Check if position is out of bounds.
     \param position Position to check.
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd) {
      if (std::floor(project(position, symmetry_planes, dd, scale, origin)) !=
          0.) {
        return 1;
      }
    }
    return 0;
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC to.
     \param state_old Previous particle state (unused).
     \return \c true if there was an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) const {
    auto projections =
        place_in_unit_cell(state.position, symmetry_planes, scale, origin);

    for (auto const &val : projections) {
      if (val != 0) {
        return 1;
      }
    }
    return 0;
  }

  /**
     \brief Translate a position according to projection in periodic cell.
     \param position Position to translate.
     \param projections Projection components in periodic cell.
  */
  template <typename Position, typename Projections = std::vector<double>>
  void translate(Position &position, Projections const &projections) const {
    geom::translate(position, symmetry_planes, projections, scale);
  }
};

/**
   \class Boundary_Periodic_SymmetryPlanes_WithOutsideInfo Geometry/Boundary.h
   "Geometry/Boundary.h"
   \brief Periodic boundaries along each symmetry plane with information about
   where position would be outside domain.
   \note State must define:
   - \c position [Container]
   - \c periodicity [Container of signed integers] (counting how many domains
   have been traveled along each dimension)
*/
template <typename SymmetryPlanes>
class Boundary_Periodic_SymmetryPlanes_WithOutsideInfo {
public:
  const SymmetryPlanes
      symmetry_planes; /**> Object describing domain symmetry planes (see
                          Geometry/SymmetryPlanes.h). */
  const double scale;  /**> Scale factor. */
  const std::vector<double> origin; /**> Coordinate origin. */

  /**
     \brief Constructor.
     \param symmetry_planes Object describing domain symmetry planes.
     \param scale Scale factor.
     \param origin Coordinate origin.
  */
  Boundary_Periodic_SymmetryPlanes_WithOutsideInfo(
      SymmetryPlanes symmetry_planes, double scale = 1.,
      std::vector<double> origin = {})
      : symmetry_planes{symmetry_planes}, scale{scale},
        origin{origin.size() == 0.
                   ? std::vector<double>(symmetry_planes.dim, 0.)
                   : origin} {}

  /**
     \brief Constructor.
     \brief Use default-constructed domain symmetry planes.
     \param scale Scale factor.
     \param origin Coordinate origin.
  */
  Boundary_Periodic_SymmetryPlanes_WithOutsideInfo(
      double scale = 1., std::vector<double> origin = {})
      : symmetry_planes{}, scale{scale}, origin{
                                             origin.size() == 0.
                                                 ? std::vector<double>(
                                                       symmetry_planes.dim, 0.)
                                                 : origin} {}

  /**
     \brief Check if position is out of bounds.
     \param position Position to check
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd) {
      if (std::floor(project(position, symmetry_planes, dd, scale, origin)) !=
          0.) {
        return 1;
      }
    }
    return 0;
  }

  /**
     \brief Enforce boundary condition.
     \param state Particle state to apply BC to.
     \param state_old Previous particle state (unused).
     \return \c true if there was an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) const {
    auto projections =
        place_in_unit_cell(state.position, symmetry_planes, scale, origin);
    op::plus_inplace(state.periodicity, projections);
    for (auto const &val : projections) {
      if (val != 0) {
        return 1;
      }
    }
    return 0;
  }

  /**
     \brief Translate a position according to projection in periodic cell.
     \param position Position to translate.
     \param projections Projection components in periodic cell.
  */
  template <typename Position, typename Projections = std::vector<double>>
  void translate(Position &position, Projections const &projections) const {
    geom::translate(position, symmetry_planes, projections, scale);
  }
};
} // namespace geom

#endif /* GEOMETRY_BOUNDARY_H */
