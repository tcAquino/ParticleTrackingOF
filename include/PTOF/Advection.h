/**
   \file PTOF/Advection.h
   \author Tomas Aquino
   \date 09/03/2022
   \brief Utilities for handling velocity field.
*/

#ifndef PTOF_ADVECTION_H
#define PTOF_ADVECTION_H

#include "PTOF/Field.h"
#include <IOobject.H>
#include <fieldTypes.H>
#include <utility>
#include <volFields.H>

namespace ptof {
/**
   \brief Get the velocity field U from the OpenFOAM case time associated with
   a mesh.
   \param geometry Domain geometry info and utilities.
   \return OpenFOAM velocity field data.
*/
template <typename Geometry> auto get_velocity_data(Geometry const &geometry) {
  auto const &mesh = geometry.mesh();
  return Foam::volVectorField{Foam::IOobject{"U", mesh.time().timeName(), mesh,
                                             Foam::IOobject::MUST_READ,
                                             Foam::IOobject::NO_WRITE},
                              mesh};
}

/**
   \brief Update velocity data, and interpolator if needed.
   \param velocity_field Velocity field as a function of state.
   \param geometry Domain geometry info and utilities.
   \param params_transport Transport parameters.
*/
template <typename VelocityField, typename Geometry,
          typename TransportParameters>
static void update_velocity_field(VelocityField &velocity_field,
                                  Geometry const &geometry,
                                  TransportParameters const &params_transport) {
  velocity_field.set(ptof::get_velocity_data(geometry));
  if (params_transport.velocity_rescaling_factor != 1.) {
    velocity_field.rescale(params_transport.velocity_rescaling_factor);
  } else if (params_transport.time_dependent_bcs) {
    velocity_field.recompute_interpolator();
  }
}

/**
   \brief Make a linear interpolator for a field using OpenFOAM interpolation.
   \param geometry Domain geometry info and utilities.
   \param field OpenFOAM vector field data.
   \return Vector field interpolator.
*/
template <typename Geometry, typename Field>
auto makeLinearVelocityInterpolator(Geometry const &geometry, Field &&field) {
  return ptof::VectorField_LinearInterpolation_OF{
      std::forward<Field>(field), geometry.locator,
      CheckOptions::Warn{}};
};

/**
   \brief Make a linear interpolator for a field using OpenFOAM interpolation.
   \param geometry Domain geometry info and utilities.
   \param field OpenFOAM vector field data.
   \param uninterpolated OpenFOAM vector field data to be added to \c field
   without interpolation.
   \return Vector field interpolator.
*/
template <typename Geometry, typename Field, typename Uninterpolated>
auto makeLinearVelocityInterpolator(Geometry const &geometry, Field &&field,
                                    Uninterpolated &&uninterpolated) {
  return ptof::VectorField_LinearInterpolation_OF{
      std::forward<Field>(field), geometry.locator,
      std::forward<Uninterpolated>(uninterpolated),
      CheckOptions::Warn{}};
};
} // namespace ptof

#endif /* PTOF_ADVECTION_H */
