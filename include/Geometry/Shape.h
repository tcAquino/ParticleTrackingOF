/**
 \file Geometry/Shape.h
 \author Tomás Aquino
 \date 07/04/2020
*/

#ifndef Shape_h
#define Shape_h

#include "General/Operations.h"
#include <cmath>
#include <vector>

/** \namespace geom Miscallaneous geometry-related objects and utilities. */
namespace geom {
/** \struct Parallelepiped  Geometry/Shape.h "Geometry/Shape.h"
 \brief Rectangle in any dimension, aligned with Cartesian axes. */
template <typename Position_t = std::vector<double>> struct Parallelepiped {
  using Position = Position_t; /**> Type of coordinate container. */

  Position corner;     /**> Lower left corner. */
  Position dimensions; /**> Sizes of rectangle along each dimension. */

  /** Constructor.
  \note Object with unassigned properties. */
  Parallelepiped() {}

  /** Constructor.
   \param corner Lower left corner.
   \param dimensions Side sizes.
  */
  Parallelepiped(Position corner, Position dimensions)
      : corner{corner}, dimensions{dimensions} {}

  /** Spatial dimension of rectangle.
  \note: Zero if side sizes have not been assigned. */
  std::size_t dim() const { return dimensions.size(); }

  /** \return \c true if \p position is inside rectangle, \c false otherwise. */
  template <typename Position> bool inside(Position const &position) const {
    auto position_helper = operation::minus(position, get_center());
    auto half_dim = get_half_dimensions();
    for (std::size_t dd = 0; dd < dim(); ++dd)
      if (std::abs(position_helper[dd]) > half_dim[dd])
        return false;
    return true;
  }

  /** \return Half the size along each dimension. */
  Position get_half_dimensions() const {
    return operation::div_scalar(dimensions, 2.);
  }

  /** \return Position of rectangle center. */
  Position get_center() const {
    return operation::plus(corner, get_half_dimensions());
  };
};

/** \struct Sphere  Geometry/Shape.h "Geometry/Shape.h"
\brief Sphere in any dimension. */
template <typename Position_t = std::vector<double>> struct Sphere {
  using Position = Position_t; /**< Type of center position. */
  Position center;             /**< Sphere center position. */
  double radius;               /**< Sphere radius. */

  /** Constructor.
  \note Object with unassigned properties. */
  Sphere() {}

  /** Constructor.
   \param center Sphere center.
   \param radius Sphere radius. */
  Sphere(Position center, double radius) : center{center}, radius{radius} {}

  /** Spatial dimension of sphere.
  \note: Zero if center has not been assigned. */
  std::size_t dim() const { return center.size(); }

  /** \return Position of sphere center. */
  Position get_center() const { return center; };

  /** \return \c true if \p position is inside spherer, \c false otherwise. */
  template <typename Position> bool inside(Position const &position) const {
    auto position_helper = operation::minus(position, center);
    if (operation::abs_sq(position_helper) > radius * radius)
      return false;
    return true;
  }
};

/** \return \c true if any \p ctrw particle is outside \p domain, \c false
 * otherwise. */
template <typename CTRW, typename Domain>
bool out_of_bounds(CTRW const &ctrw, Domain const &domain) {
  for (auto const &part : ctrw.particles())
    if (domain.out_of_bounds(part.state_new().position))
      return 1;
  return 0;
}

/** Domain with arbitrary shape, with rectangular and/or spherical inclusions in
 * any spatial dimension. */
template <typename Shape> struct Domain {
  using DomainShape = Shape; /**< Bounding domain shape. */

  /** Constructor. */
  Domain() {}

  /** Constructor.
   \param box Domain shape. */
  Domain(DomainShape box) : box{box} {}

  /** Constructor.
  \param box Domain shape.
  \param parallelepipeds Rectangular inclusions.
  \param spheres Spherical inclusions. */
  Domain(DomainShape box, std::vector<Parallelepiped<>> parallelepipeds,
         std::vector<Sphere<>> spheres)
      : box{box}, parallelepipeds{parallelepipeds}, spheres{spheres} {}

  Shape box;                                     /**< Bounding domain shape. */
  std::vector<Parallelepiped<>> parallelepipeds; /**< Rectangular inclusions. */
  std::vector<Sphere<>> spheres;                 /**< Spherical inclusions. */

  /** \return Domain side sizes. */
  std::vector<double> dimensions() const {
    if constexpr (std::is_same_v<Shape, Sphere<>>)
      return std::vector<double>(box.dim(), 2. * box.radius);
    else
      return box.dimensions;
  }

  /** \return \c true if \c position is outside domain, \c false otherwise .*/
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    if (!box.inside(position))
      return 1;
    for (auto const &shape : parallelepipeds)
      if (shape.inside(position))
        return 1;
    for (auto const &shape : spheres)
      if (shape.inside(position))
        return 1;

    return 0;
  }
};

/** \brief Reflect off sphere exterior (in any dimension).
 \param position_new Position inside sphere.
 \param position_old Position outside sphere.
 \param sphere Sphere to reflect off.
 \param reflected Resulting reflected position, computed here.
 \param contact_point Resulting contact point with sphere, computed here.
 \note - It is safe to pass the same containers for inputs and outputs.
 - Output vectors must be passed with correct size.
 */
template <typename Position, typename Sphere>
void reflectOffSphere_outsideToInside(Position const &position_new,
                                      Position const &position_old,
                                      Sphere const &sphere, Position &reflected,
                                      Position &contact_point) {
  // Vector from outside position to inside position
  auto insideminusoutside = operation::minus(position_new, position_old);

  // Vector from outside position to sphere center
  auto outsideminuscenter = operation::minus(position_old, sphere.center);

  // Find contact point with sphere
  double coeff_a = operation::abs_sq(insideminusoutside);
  double coeff_b = 2. * operation::dot(insideminusoutside, outsideminuscenter);
  double coeff_c = operation::dot(outsideminuscenter, outsideminuscenter) -
                   sphere.radius * sphere.radius;
  double fraction_to_intersection =
      -(coeff_b + std::sqrt(coeff_b * coeff_b - 4. * coeff_a * coeff_c)) /
      (2. * coeff_a);
  operation::linearOp(fraction_to_intersection, insideminusoutside,
                      position_old, contact_point);

  // Vector from contact point to inside position
  auto leftover = operation::times_scalar(1. - fraction_to_intersection,
                                          insideminusoutside);

  // Compute unit normal to sphere surface at contact
  auto normal = operation::minus(contact_point, sphere.center);
  operation::div_scalar_InPlace(normal, sphere.radius);

  // Compute reflected position
  operation::linearOp(-2. * operation::dot(leftover, normal), normal,
                      position_new, reflected);
}

/** \brief Reflect off sphere exterior (in any dimension).
\param position_new Position inside sphere.
\param position_old Position outside sphere.
\param sphere Sphere to reflect off.
\param reflected Resulting reflected position, computed here.
\note - It is safe to pass the same containers for inputs and outputs.
 - Output vectors must be passed with correct size. */
template <typename Position, typename Sphere>
void reflectOffSphere_outsideToInside(Position const &position_new,
                                      Position const &position_old,
                                      Sphere const &sphere,
                                      Position &reflected) {
  auto contact_point(position_new.size());
  reflectOffSphere_outsideToInside(position_new, position_old, sphere,
                                   reflected, contact_point);
}

/** \brief Reflect off sphere exterior (in any dimension).
 \details Trajectory is approximated as a straight line.
 \param position_new Position inside sphere.
 \param position_old Position outside sphere.
 \param velocity_old Velocity at \c position_old.
 \param sphere Sphere to reflect off.
 \param position_reflected Resulting reflected position, computed here.
 \param velocity_reflected Resulting reflected velocity, computed here.
 \param contact_point Resulting contact point with sphere, computed here.
 \param contact_velocity Resulting velocity at contact point with sphere,
 computed here. \param time_to_contact Resulting time until contact, computed
 here. \param time_step Time step. \param acceleration Acceleration vector.
 \param restitution_coeff_norm Fraction of normal velocity convserved in
 collision. \param restitution_coeff_tang Fraction of tangential velocity
 convserved in collision. \param radius Particle radius. \note - It is safe to
 pass the same containers for inputs and outputs.
 - Output vectors must be passed with correct size.
*/
template <typename Position, typename Sphere>
void reflectOffSphere_velocity_outsideToInside(
    Position const &position_new, Position const &position_old,
    Position const &velocity_old, Sphere const &sphere,
    Position &position_reflected, Position &velocity_reflected,
    Position &contact_point, Position &contact_velocity,
    double &time_to_contact, double time_step, Position const &acceleration,
    double restitution_coeff_norm = 1., double restitution_coeff_tang = 1.,
    double radius = 0.) {
  double effective_radius = sphere.radius + radius;

  // Vector from outside position to inside position
  auto insideminusoutside = operation::minus(position_new, position_old);

  // Vector from outside position to sphere center
  auto outsideminuscenter = operation::minus(position_old, sphere.center);

  // Find contact point with sphere
  double coeff_a = operation::abs_sq(insideminusoutside);
  double coeff_b = 2. * operation::dot(insideminusoutside, outsideminuscenter);
  double coeff_c = operation::dot(outsideminuscenter, outsideminuscenter) -
                   effective_radius * effective_radius;
  double fraction_to_intersection =
      -(coeff_b + std::sqrt(coeff_b * coeff_b - 4. * coeff_a * coeff_c)) /
      (2. * coeff_a);
  operation::linearOp(fraction_to_intersection, insideminusoutside,
                      position_old, contact_point);

  // Vector from outside position to contact point
  auto to_contact =
      operation::times_scalar(fraction_to_intersection, insideminusoutside);

  // Compute unit normal to sphere surface at contact
  auto normal = operation::minus(contact_point, sphere.center);
  operation::div_scalar_InPlace(normal, effective_radius);

  time_to_contact = std::sqrt(operation::abs_sq(to_contact) /
                              operation::abs_sq(velocity_old));
  operation::linearOp(time_to_contact, acceleration, contact_velocity,
                      contact_velocity);
  operation::linearOp(-2. * operation::dot(contact_velocity, normal), normal,
                      velocity_old, contact_velocity);

  auto normal_vel =
      operation::times_scalar(operation::dot(contact_velocity, normal), normal);
  auto tangential_vel = operation::minus(contact_velocity, normal_vel);
  operation::times_scalar_InPlace(restitution_coeff_norm, normal_vel);
  operation::times_scalar_InPlace(restitution_coeff_tang, tangential_vel);
  operation::plus(normal_vel, tangential_vel, contact_velocity);

  // Compute reflected position and velocity
  double time_leftover = time_step - time_to_contact;
  operation::linearOp(time_leftover, contact_velocity, contact_point,
                      position_reflected);
  operation::linearOp(time_leftover, acceleration, contact_velocity,
                      velocity_reflected);
}

/** \brief Reflect off sphere exterior (in any dimension).
 \details Trajectory is approximated as a straight line.
 \param position_new Position inside sphere.
 \param position_old Position outside sphere.
 \param velocity_old Velocity at \c position_old.
 \param sphere Sphere to reflect off.
 \param position_reflected Resulting reflected position, computed here.
 \param velocity_reflected Resulting reflected velocity, computed here.
 \param time_step Time step.
 \param acceleration Acceleration vector.
 \param restitution_coeff Fraction of velocity convserved in collision.
 \param radius Particle radius.
 \note - It is safe to pass the same containers for inputs and outputs.
 - Output vectors must be passed with correct size.
*/
template <typename Position, typename Sphere, typename State,
          typename Acceleration>
void reflectOffSphere_velocity_outsideToInside(
    Position const &position_new, Position const &position_old,
    Position const &velocity_old, Sphere const &sphere,
    Position &position_reflected, Position &velocity_reflected,
    double time_step, Position const &acceleration,
    double restitution_coeff = 1., double radius = 0.) {
  auto contact_point(position_new.size());
  auto contact_velocity(position_new.size());
  double time_to_contact;
  reflectOffSphere_velocity_outsideToInside(
      position_new, position_old, velocity_old, sphere, position_reflected,
      velocity_reflected, contact_point, contact_velocity, time_to_contact,
      time_step, acceleration, restitution_coeff, restitution_coeff, radius);
}
} // namespace geom

#endif /* Shape_h */
