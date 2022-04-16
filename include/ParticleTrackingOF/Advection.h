//
//  Advection.h
//
//  Created by Tomás Aquino on 09/03/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Advection_OF_h
#define Advection_OF_h

#include <utility>
#include <fieldTypes.H>
#include <IOobject.H>
#include <meshSearch.H>
#include "ParticleTrackingOF/Field.h"

namespace ptof
{
  // Rescale a volumetric field (set of values associated with mesh cells)
  // to a given average value
  template <typename Field, typename Mesh>
  void rescale_to_average(Field& field, Mesh const& mesh, double average)
  {
    Foam::scalar mesh_volume
      = Foam::sum(mesh.cellVolumes());
    auto average_weighted_vector_data
      = Foam::sum(field*mesh.cellVolumes());
    auto average_data
      = Foam::mag(average_weighted_vector_data)/mesh_volume;
    field *= average/average_data;
  }
  
  // Get the velocity field U
  // from the OpenFOAM case time associated with a mesh
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
  
  // Get velocity data and rescale it to a given average
  template <typename Mesh>
  auto get_velocity_data_rescaled(Mesh const& mesh, double average)
  {
    auto data = get_velocity_data(mesh);
    rescale_to_average(data, mesh, average);
    return data;
  }
  
  // Make a linear interpolator for a field
  // using OpenFOAM interpolation
  template <typename Geometry, typename Field>
  auto makeLinearInterpolator
  (Geometry const& geometry, Field&& field)
  {
    return ptof::Field_LinearInterpolation_OF{
      std::forward<Field>(field),
      geometry.locator,
      FieldOptions::Warn{}
    };
  };
}

#endif /* Advection_OF_h */
