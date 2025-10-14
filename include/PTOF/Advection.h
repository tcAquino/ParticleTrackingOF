/**
   \file PTOF/Advection.h
   \author Tomas Aquino
   \date 09/03/2022
   \brief Utilities for handling velocity field.
*/

#ifndef PTOF_ADVECTION_H
#define PTOF_ADVECTION_H

#include "General/Meta.h"
#include "PTOF/CheckOptions.h"
#include "PTOF/Field.h"
#include "PTOF/Useful.h"
#include <IOobject.H>
#include <fieldTypes.H>
#include <utility>
#include <volFields.H>

namespace ptof {
/**
   \brief Get the velocity field U from the OpenFOAM case time associated with
   a mesh, at a prescribed time.
   \return OpenFOAM velocity field data.
*/
template <typename Mesh, bool advection = true>
auto get_velocity_data(
    Mesh const &mesh, Foam::word const &timeName,
    meta::Selector<bool, advection> = meta::Selector<bool, true>{}) {
  if constexpr (advection) {
    return Foam::volVectorField{Foam::IOobject{"U", timeName, mesh,
                                               Foam::IOobject::MUST_READ,
                                               Foam::IOobject::NO_WRITE},
                                mesh};
  } else {
    return meta::Empty{};
  }
}

/**
   \brief Get the velocity field U from the OpenFOAM case time associated with
   a mesh, at a prescribed time.
   \return OpenFOAM velocity field data.
*/
template <typename Mesh, bool advection = true>
auto get_velocity_data(Mesh const &mesh, Foam::scalar time,
                       meta::Selector<bool, advection> get_velocity =
                           meta::Selector<bool, true>{}) {
  return get_velocity_data(mesh, closest_time_name(mesh, time), get_velocity);
}

/**
   \brief Get the velocity field U from the OpenFOAM case time associated with
   a mesh.
   \return OpenFOAM velocity field data.
*/
template <typename Mesh, bool advection = true>
auto get_velocity_data(Mesh const &mesh,
                       meta::Selector<bool, advection> get_velocity =
                           meta::Selector<bool, true>{}) {
  return get_velocity_data(mesh, mesh.time().timeName(), get_velocity);
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
  velocity_field.set(ptof::get_velocity_data(geometry.mesh()));
  if (params_transport.velocity_rescaling_factor != 1.) {
    velocity_field.rescale(params_transport.velocity_rescaling_factor);
  }
  if (params_transport.time_dependent_bcs) {
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
  return ptof::VectorField_Interpolation{
      std::forward<Field>(field), geometry.locator,
      InterpolationTypes::Linear{}, CheckOptions::Warn{}};
};
} // namespace ptof

#endif /* PTOF_ADVECTION_H */
