/**
 \file PTOF/Measure.h
 \author Tomás Aquino
 \date 07/03/2022
*/

#ifndef PTOF_MEASURE_H
#define PTOF_MEASURE_H

#include <map>
#include <string>

namespace ptof {
/** \struct Measure PTOF/Output.h "PTOF/Output.h"
 * \brief Keep track of measurement types */
struct Measure {
  /** \enum Type
   *  \brief Implemented types. */
  enum class Type {
    position,               /**< Time, particle tags, positions, and masses. */
    position_in_regions,    /**< Time, particle tags, positions, and masses in
                               regions specified by masks. */
    position_mean,          /**< Time and mean position. */
    position_second_moment, /**< Time and position second moment. */
    position_variance,      /**< Time and position variance. */
    mass,                   /**< Time and total mass. */
    mass_in_regions, /**< Time and total mass in regions speciefied by masks. */
    velocity,        /**< Time, particle tags, and local velocities. */
    velocity_mean,   /**< Time and mean of velocity field over particles. */
    velocity_gradient, /**< Time, particle tags, and local velocity gradients.
                        */
    velocity_gradient_mean, /**< Time and mean of velocity gradient field over
                             * particles.
                             */
    scalar_field, /**< Time, particle tags, and local values of scalar field. */
    vector_field, /**< Time, particle tags, and local values of vector field. */
    tensor_field, /**< Time, particle tags and local values of tensor field. */
    scalar_field_mean, /**< Time and mean of scalar field over particles. */
    vector_field_mean, /**< Time and mean of vector field over particles. */
    tensor_field_mean, /**< Time and mean of tensor field over particles. */
    position_periodic, /**< Time, particle tags, true positions accounting for
                          periodicity, and masses. */
    position_in_regions_periodic, /**< Time, particle tags, true positions
                                     accounting for periodicity, and masses in
                                     regions specified by masks. */
    position_mean_periodic, /**< Time and true mean position accounting for
                               periodicity. */
    position_second_moment_periodic, /**< Time and true position second moment
                                        accounting for periodicity. */
    position_variance_periodic, /**< Time and true position variance accounting
                                   for periodicity. */
    first_crossing_time, /**< First crossing time, tag, and mass for specified
                            position along specified dimension. */
    absorption_time /**< Particle absorption times, tags, and masses at end of
                       dynamics. */
  };

  /** \return Type from name. */
  static auto type(std::string const &name) { return name_to_type.at(name); }

  /** \return Name from type. */
  static auto name(Type type) { return type_to_name.at(type); }

  /** \return \c true if criterion name exists, \c false otherwise. */
  static bool exists(std::string const &name) {
    return name_to_type.count(name);
  }

  /** Map of names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"position", Type::position},
      {"position_in_regions", Type::position_in_regions},
      {"position_mean", Type::position_mean},
      {"position_second_moment", Type::position_second_moment},
      {"position_variance", Type::position_variance},
      {"mass", Type::mass},
      {"mass_in_regions", Type::mass_in_regions},
      {"velocity", Type::velocity},
      {"velocity_mean", Type::velocity_mean},
      {"velocity_gradient", Type::velocity_gradient},
      {"velocity_gradient_mean", Type::velocity_gradient_mean},
      {"scalar_field", Type::scalar_field},
      {"vector_field", Type::vector_field},
      {"tensor_field", Type::tensor_field},
      {"scalar_field_mean", Type::scalar_field_mean},
      {"vector_field_mean", Type::vector_field_mean},
      {"tensor_field_mean", Type::scalar_field_mean},
      {"position_periodic", Type::position_periodic},
      {"position_in_regions_periodic", Type::position_in_regions_periodic},
      {"position_mean_periodic", Type::position_mean_periodic},
      {"position_second_moment_periodic",
       Type::position_second_moment_periodic},
      {"position_variance_periodic", Type::position_variance_periodic},
      {"first_crossing_time", Type::first_crossing_time},
      {"absorption_time", Type::absorption_time}};

  /** Map of types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::position, "position"},
      {Type::position_in_regions, "position_in_regions"},
      {Type::position_mean, "position_mean"},
      {Type::position_second_moment, "position_second_moment"},
      {Type::position_variance, "position_variance"},
      {Type::mass, "mass"},
      {Type::mass_in_regions, "mass_in_regions"},
      {Type::velocity, "velocity"},
      {Type::velocity_mean, "velocity_mean"},
      {Type::velocity_gradient, "velocity_gradient"},
      {Type::velocity_gradient_mean, "velocity_gradient_mean"},
      {Type::scalar_field, "scalar_field"},
      {Type::vector_field, "vector_field"},
      {Type::tensor_field, "tensor_field"},
      {Type::scalar_field_mean, "scalar_field_mean"},
      {Type::vector_field_mean, "vector_field_mean"},
      {Type::scalar_field_mean, "tensor_field_mean"},
      {Type::position_periodic, "position_periodic"},
      {Type::position_in_regions_periodic, "position_in_regions_periodic"},
      {Type::position_mean_periodic, "position_mean_periodic"},
      {Type::position_second_moment_periodic,
       "position_second_moment_periodic"},
      {Type::position_variance_periodic, "position_variance_periodic"},
      {Type::first_crossing_time, "first_crossing_time"},
      {Type::absorption_time, "absorption_time"}};
};

/** \struct MeasureSpacing PTOF/Output.h "PTOF/Output.h"
 *  \brief Keep track of names and types of measure spacing */
struct MeasureSpacing {
  /** \enum Type
   *  \brief Implemented types. */
  enum class Type {
    linear, /**< Fixed number of linearly-spaced measurements, starting and
               ending at specified  times */
    log,    /**< Fixed number of log-spaced measurements, starting and ending at
               specified times.      */
    step    /**< At fixed time intervals, starting at specified time.    */
  };

  /** \return Type from name. */
  static auto type(std::string const &name) { return name_to_type.at(name); }

  /** \return Name from type. */
  static auto name(Type type) { return type_to_name.at(type); }

  /** \return \c true if name exists, \c false otherwise. */
  static bool exists(std::string const &name) {
    return name_to_type.count(name);
  }

  /** Map of names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"linear", Type::linear}, {"log", Type::log}, {"step", Type::step}};

  /** Map of types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::linear, "linear"}, {Type::log, "log"}, {Type::step, "step"}};
};

/** \struct MeasureSpacingUnits PTOF/Output.h "PTOF/Output.h"
 * \brief Keep track of names and types of measure spacing units (scaling) */
struct MeasureSpacingUnits {
  /** \enum Type
   *  \brief Implemented types. */
  enum class Type {
    diffusion, /**< Units of diffusion scale */
    advection, /**< Units of advection scale */
    reaction,  /**<  Units of reaction scale*/
    arbitrary  /**< Arbitrary units*/
  };

  /** \return Type from name. */
  static auto type(std::string const &name) { return name_to_type.at(name); }

  /** \return Name from type. */
  static auto name(Type type) { return type_to_name.at(type); }

  /** \return \c true if criterion name exists, \c false otherwise. */
  static bool exists(std::string const &name) {
    return name_to_type.count(name);
  }

  /** Map of names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"diffusion", Type::diffusion},
      {"advection", Type::advection},
      {"reaction", Type::reaction},
      {"arbitrary", Type::arbitrary}};

  /** Map of types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::diffusion, "diffusion"},
      {Type::advection, "advection"},
      {Type::reaction, "reaction"},
      {Type::arbitrary, "arbitrary"}};
};
} // namespace ptof

#endif /* PTOF_MEASURE_H */
