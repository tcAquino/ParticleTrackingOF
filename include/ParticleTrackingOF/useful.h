//
//  useful_OF.h
//
//  Created by Tomás Aquino on 17/02/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef useful_OF_h
#define useful_OF_h

#include <type_traits>
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
  
  // Small offset forward given current face and direction
  template <typename MeshSearch>
  Foam::vector offset_forward_face
  (Foam::point const& begin,
   Foam::label face,
   Foam::vector const& direction,
   MeshSearch const& mesh_search)
  {
    Foam::label owner_cell = mesh_search.mesh().faceOwner()[face];
    Foam::point const& center = mesh_search.mesh().
      cellCentres()[owner_cell];
    Foam::scalar typ_dim = Foam::mag(center - begin);
    
    return begin + mesh_search.tol_*typ_dim*
      direction/Foam::mag(direction);
  }
  
  // Small offset backward given current face and direction
  template <typename MeshSearch>
  Foam::vector offset_backward_face
  (Foam::point const& begin,
   Foam::label face,
   Foam::vector const& direction,
   MeshSearch const& mesh_search)
  {
    Foam::label owner_cell = mesh_search.mesh().faceOwner()[face];
    Foam::point const& center = mesh_search.mesh().
      cellCentres()[owner_cell];
    Foam::scalar typ_dim = Foam::mag(center - begin);
    
    return begin - mesh_search.tol_*typ_dim*
      direction/Foam::mag(direction);
  }
  
  // Small offset forward
  // given current cell and direction
  template <typename MeshSearch>
  Foam::vector offset_forward_cell
  (Foam::point const& begin,
   Foam::label cell,
   Foam::vector const& direction,
   MeshSearch const& mesh_search)
  {
    Foam::point const& center = mesh_search.mesh().cellCentres()[cell];
    Foam::scalar typ_dim = Foam::mag(center - begin);
    
    return begin - mesh_search.tol_*typ_dim*
      direction/Foam::mag(direction);
  }
  
  // Small offset backward
  // given current cell and direction
  template <typename MeshSearch>
  Foam::vector offset_backward_cell
  (Foam::point const& begin,
   Foam::label cell,
   Foam::vector const& direction,
   MeshSearch const& mesh_search)
  {
    Foam::point const& center = mesh_search.mesh().cellCentres()[cell];
    Foam::scalar typ_dim = Foam::mag(center - begin);
    
    return begin - mesh_search.tol_*typ_dim*
      direction/Foam::mag(direction);
  }
  
  // Check for existence of periodicity info
  // (to use with State)
  template <typename State, typename = int>
  struct Has_periodicity : std::false_type { };
  template <typename State>
  struct Has_periodicity
  <State, decltype((void) State::periodicity, 0)>
  : std::true_type { };
}

#endif /* useful_OF_h */
