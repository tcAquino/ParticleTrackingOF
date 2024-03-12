/**
 \file PTOF/Advection.h
 \author Tomás Aquino
 \date 09/03/2022
*/

#ifndef PTOF_ADVECTION_H
#define PTOF_ADVECTION_H

#include <utility>
#include <fieldTypes.H>
#include <IOobject.H>
#include <meshSearch.H>
#include "PTOF/Field.h"

namespace ptof
{ 
  /**
   \brief Get the velocity field U from the OpenFOAM case time associated with a mesh.
   \param mesh OpenFOAM mesh.
   \return OpenFOAM velocity field data.
   */
  template <typename Mesh>
  auto get_velocity_data(Mesh const& mesh)
  {
    return Foam::volVectorField{
      Foam::IOobject{
        "U",
        mesh.time().timeName(),
        mesh,
        Foam::IOobject::MUST_READ,
        Foam::IOobject::NO_WRITE },
      mesh };
  }
  
  /**
   \brief Get the velocity field U from the OpenFOAM case time associated with a mesh and rescale it to a given average.
   \param mesh OpenFOAM mesh.
   \param average Absolute value of velocity field average.
   \return OpenFOAM velocity field data.
  */
  template <typename Mesh>
  auto get_velocity_data_rescaled(Mesh const& mesh, double average)
  {
    auto data = get_velocity_data(mesh);
    rescale_to_average(data, mesh, average);
    return data;
  }
  
  /**
   \brief Make a linear interpolator for a field using OpenFOAM interpolation.
   \param geometry Domain geometry info and utilities.
   \param field OpenFOAM vector field data
   \return Vector field interpolator
  */
  template <typename Geometry, typename Field>
  auto makeLinearInterpolator
  (Geometry const& geometry, Field&& field)
  {
    return ptof::VectorField_LinearInterpolation_OF{
      std::forward<Field>(field),
      geometry.locator,
      CheckOptions::Warn{} };
  };
  
  /**
   \brief Make a linear interpolator for a field using OpenFOAM interpolation.
   \param geometry Domain geometry info and utilities.
   \param field OpenFOAM vector field data
   \param uninterpolated OpenFOAM vector field data to be added to \c field without interpolation
   \return Vector field interpolator
  */
  template <typename Geometry, typename Field, typename Uninterpolated>
  auto makeLinearInterpolator
  (Geometry const& geometry, Field&& field, Uninterpolated&& uninterpolated)
  {
    return ptof::VectorField_LinearInterpolation_OF{
      std::forward<Field>(field),
      geometry.locator,
      std::forward<Uninterpolated>(uninterpolated),
      CheckOptions::Warn{} };
  };
}

#endif /* PTOF_ADVECTION_H */
