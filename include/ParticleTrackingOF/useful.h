//
//  useful_OF.h
//  ParticleTracking_OpenFOAM
//
//  Created by Tomás Aquino on 17/02/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef useful_OF_h
#define useful_OF_h

#include <fieldTypes.H>

namespace ptof
{
  // Make 3D point from 2D point
  auto make_point(Foam::Vector2D<Foam::scalar> const& point)
  {
    return Foam::point{ point[0], point[1], 0. };
  }

  // Make 3D point from 1D point
  auto make_point(Foam::scalar point)
  {
    return Foam::point{ point, 0., 0. };
  }

  // Make 3D point 3D point (simple copy)
  auto make_point(Foam::point const& point)
  { return point; }
  
  // Unit normal at boundary face, pointing outward
  template <typename Mesh>
  auto unit_normal_outward
  (Foam::label face, Mesh const& mesh)
  {
    return mesh.faces()[face].unitNormal(mesh.points());
  }
  
  // Unit normal at boundary face, pointing inward
  template <typename Mesh>
  auto unit_normal_inward
  (Foam::label face, Mesh const& mesh)
  {
    return -mesh.faces()[face].unitNormal(mesh.points());
  }
  
  // Area (magnitude) of a face
  template <typename Mesh>
  auto face_area
  (Foam::label face, Mesh const& mesh)
  {
    return mesh.faces()[face].mag(mesh.points());
  }
  
  // Center of a face
  template <typename Mesh>
  auto face_center
  (Foam::label face, Mesh const& mesh)
  {
    return mesh.faces()[face].centre(mesh.points());
  }
  
  // Volume of a cell
  template <typename Mesh>
  auto cell_volume
  (Foam::label cell, Mesh const& mesh)
  {
    return mesh.cellVolumes()[cell];
  }
  
  // Center of a cell
  template <typename Mesh>
  auto cell_center
  (Foam::label cell, Mesh const& mesh)
  {
    return mesh.cellCentres()[cell];
  }
}

#endif /* useful_OF_h */
