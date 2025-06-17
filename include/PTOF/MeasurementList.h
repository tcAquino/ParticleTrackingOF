/**
   \file PTOF/MeasurementList.h
   \author Tomás Aquino
   \date 07/03/2022
   \brief Measurement types.
*/

#ifndef PTOF_MEASUREMENTLIST_H
#define PTOF_MEASUREMENTLIST_H

#include <map>
#include <string>

/**
   \struct MeasurementList PTOF/MeasurementList.h "PTOF/MeasurementList.h"
   \brief Names and types of measurements.
*/
struct MeasurementList {
  /**
     \enum Type
     \brief Implemented types.
  */
  enum class Type {
    position,               /**< Time, particle tags, positions, and masses. */
    position_in_regions,    /**< Time, particle tags, positions, and masses in
                               regions specified by masks. */
    position_mean,          /**< Time and mean position. */
    position_second_moment, /**< Time and position second moment. */
    position_nth_moment,    /**< Time and position nth moment. */
    position_moment,   /**< Time and position moment with specified exponents
                          along each dimension. */
    position_variance, /**< Time and position variance. */
    mass,              /**< Time and total mass. */
    mass_absorbed,     /**< Time and total absorbed mass. */
    mass_adsorbed,     /**< Time and total adsorbed mass. */
    mass_in_regions, /**< Time and total mass in regions speciefied by masks. */
    velocity,        /**< Time, particle tags, and local velocities. */
    velocity_mean,   /**< Time and mean of velocity field over particles. */
    velocity_gradient, /**< Time, particle tags, and local velocity gradients.
                        */
    velocity_gradient_mean, /**< Time and mean of velocity gradient field over
                               particles. */
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
    position_nth_moment_periodic,    /**< Time and true position nth moment
                                           accounting for periodicity. */
    position_moment_periodic,   /**< Time and position moment with specified
                                   exponents along each dimension. */
    position_variance_periodic, /**< Time and true position variance accounting
                                   for periodicity. */
    first_crossing_time, /**< First crossing time, tag, and mass for specified
                            position along specified dimension. */
    adsorbed_position, /**< Time, tag, position, and mass of adsorbed particles.
                        */
    adsorbed_position_periodic, /**< Time, tag, position, and mass of adsorbed
                                   particles. */
    absorption_time, /**< Particle absorption times, tags, and masses at end of
                       dynamics. */
    absorption_time_patch,    /**< Particle absorption times, particle tags,
                  particle masses, and absorption patch at end of dynamics. */
    absorption_time_position, /**< Particle absorption times, particle tags,
       particle masses, and particle position at end of dynamics. */
    absorption_time_patch_position,    /**< Particle absorption times, particle
                                         tags, particle masses, closest boundary
                                         patch, and particle position at end of
                                         dynamics. */
    absorption_time_position_periodic, /**< Particle absorption times, particle
       tags, particle masses, and particle position accounting for periodicity
       at end of dynamics. */
    absorption_time_patch_position_periodic /**< Particle absorption times,
                                      particle tags, particle masses, closest
                                      boundary patch, and particle position
                                      accounting for periodicity at end of
                                      dynamics. */
  };

  /**
     \brief Type from name.
     \param name Boundary condition name.
     \return Boundary condition type.
  */
  static auto type(std::string const &name) { return name_to_type.at(name); }

  /**
     \brief Name from type.
     \param type Boundary condition type.
     \return Boundary condition name.
  */
  static auto name(Type type) { return type_to_name.at(type); }

  /**
     \brief Check if name exists.
     \param name condition name.
     \return \c true if name exists, \c false otherwise.
  */
  static bool contains(std::string const &name) {
    return name_to_type.count(name);
  }

  /** Map of names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"position", Type::position},
      {"position_in_regions", Type::position_in_regions},
      {"position_mean", Type::position_mean},
      {"position_second_moment", Type::position_second_moment},
      {"position_nth_moment", Type::position_nth_moment},
      {"position_moment", Type::position_moment},
      {"position_variance", Type::position_variance},
      {"mass", Type::mass},
      {"mass_absorbed", Type::mass_absorbed},
      {"mass_adsorbed", Type::mass_adsorbed},
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
      {"position_nth_moment_periodic", Type::position_nth_moment},
      {"position_moment_periodic", Type::position_moment_periodic},
      {"position_variance_periodic", Type::position_variance_periodic},
      {"first_crossing_time", Type::first_crossing_time},
      {"adsorbed_position", Type::adsorbed_position},
      {"adsorbed_position_periodic", Type::adsorbed_position_periodic},
      {"absorption_time", Type::absorption_time},
      {"absorption_time_patch", Type::absorption_time_patch},
      {"absorption_time_position", Type::absorption_time_position},
      {"absorption_time_patch_position", Type::absorption_time_patch_position},
      {"absorption_time_position_periodic",
       Type::absorption_time_position_periodic},
      {"absorption_time_patch_position_periodic",
       Type::absorption_time_patch_position_periodic}};

  /** Map of types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::position, "position"},
      {Type::position_in_regions, "position_in_regions"},
      {Type::position_mean, "position_mean"},
      {Type::position_second_moment, "position_second_moment"},
      {Type::position_nth_moment, "position_nth_moment"},
      {Type::position_moment, "position_moment"},
      {Type::position_variance, "position_variance"},
      {Type::mass, "mass"},
      {Type::mass_absorbed, "mass_absorbed"},
      {Type::mass_adsorbed, "mass_adsorbed"},
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
      {Type::position_nth_moment_periodic, "position_nth_moment_periodic"},
      {Type::position_moment_periodic, "position_moment_periodic"},
      {Type::position_variance_periodic, "position_variance_periodic"},
      {Type::first_crossing_time, "first_crossing_time"},
      {Type::adsorbed_position, "adsorbed_position"},
      {Type::adsorbed_position_periodic, "adsorbed_position_periodic"},
      {Type::absorption_time, "absorption_time"},
      {Type::absorption_time_patch, "absorption_time_patch"},
      {Type::absorption_time_position, "absorption_time_position"},
      {Type::absorption_time_position_periodic,
       "absorption_time_position_periodic"},
      {Type::absorption_time_patch_position_periodic,
       "absorption_time_patch_position_periodic"}};
};

#endif /* PTOF_MEASUREMENTLIST_H */
