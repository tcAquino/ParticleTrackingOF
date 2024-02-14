/**
* \file PTOF/Advection.h
* \author Tomás Aquino
* \date 09/03/2022
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
  /** Get the velocity field U from the OpenFOAM case time associated with a mesh.*/  
  template <typename Mesh>
  auto get_velocity_data(Mesh const& mesh)
  {
    return Foam::volVectorField{
      Foam::IOobject{
        "U",
        mesh.time().timeName(),
        mesh,
        Foam::IOobject::MUST_READ,
        Foam::IOobject::NO_WRITE
      },
      mesh };
  }
  
  /** Get velocity data and rescale it to a given average. */  
  template <typename Mesh>
  auto get_velocity_data_rescaled(Mesh const& mesh, double average)
  {
    auto data = get_velocity_data(mesh);
    rescale_to_average(data, mesh, average);
    return data;
  }
  
  /** Make a linear interpolator for a field using OpenFOAM interpolation.*/  
  template <typename Geometry, typename Field>
  auto makeLinearInterpolator
  (Geometry const& geometry, Field&& field)
  {
    return ptof::VectorField_LinearInterpolation_OF{
      std::forward<Field>(field),
      geometry.locator,
      FieldOptions::Warn{} };
  };
  
  /** Make a linear interpolator for a field using OpenFOAM interpolation.*/
  template <typename Geometry, typename Field>
  auto makeLinearInterpolator
  (Geometry const& geometry, Field&& field, std::size_t thread)
  {
    return ptof::VectorField_LinearInterpolation_OF{
      std::forward<Field>(field),
      geometry.locator[thread],
      FieldOptions::Warn{} };
  };
}

#endif /* PTOF_ADVECTION_H */
