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
  
  // Sometimes the face center associated with a mesh face
  // is not considered within the cell,
  // verify this
  template <typename Mesh, typename MeshSearch>
  bool face_center_is_in_cell
  (Foam::label face, Mesh const& mesh, MeshSearch const& mesh_search)
  {
    return mesh_search.findCell(face_center(face, mesh_search.mesh()))
      == mesh.owner()[face];
  }
  
  // Small offset given current face, direction, and begin point
  template <typename MeshSearch>
  Foam::vector offset_face
  (Foam::point const& begin,
   Foam::label face,
   Foam::vector const& direction,
   MeshSearch const& mesh_search)
  {
    Foam::label owner_cell = mesh_search.mesh().faceOwner()[face];
    Foam::point const& cell_center = mesh_search.mesh().
      cellCentres()[owner_cell];
    Foam::scalar typ_dim = Foam::mag(cell_center
                                     - face_center(face,
                                                   mesh_search.mesh()));
    
    return mesh_search.tol_*typ_dim*
      direction/Foam::mag(direction);
  }
  
  // Small offset along face normal (outward)
  // given current face
  // with face center as begin point
  template <typename MeshSearch>
  Foam::vector offset_face
  (Foam::label face,
   MeshSearch const& mesh_search)
  {
    return offset_face(face_center(face, mesh_search.mesh()),
                       face,
                       unit_normal_outward(face, mesh_search.mesh()),
                       mesh_search);
  }
  
  // Small offset forward from begin
  // given current face and direction
  template <typename MeshSearch>
  Foam::vector offset_forward_face
  (Foam::point const& begin,
   Foam::label face,
   Foam::vector const& direction,
   MeshSearch const& mesh_search)
  {
    return begin
      + offset_face(begin, face, direction, mesh_search);
  }
  
  // Small offset backward from begin
  // given current face and direction
  template <typename MeshSearch>
  Foam::vector offset_backward_face
  (Foam::point const& begin,
   Foam::label face,
   Foam::vector const& direction,
   MeshSearch const& mesh_search)
  {
    return begin
      - offset_face(begin, face, direction, mesh_search);
  }
  
  // Small offset along face normal (outward)
  // given current face
  // with face center as begin point
  template <typename MeshSearch>
  Foam::vector offset_outward_face
  (Foam::label face,
   MeshSearch const& mesh_search)
  {
    return face_center(face, mesh_search.mesh())
      + offset_face(face, mesh_search);
  }
  
   // Small offset along face normal (inward)
   // given current face
   // with face center as begin point
   template <typename MeshSearch>
   Foam::vector offset_inward_face
   (Foam::label face,
    MeshSearch const& mesh_search)
   {
     return face_center(face, mesh_search.mesh())
      - offset_face(face, mesh_search);
   }
  
  // Small offset
  // given current cell, direction, and begin point
  template <typename MeshSearch>
  Foam::vector offset_cell
  (Foam::point const& begin,
   Foam::label cell,
   Foam::vector const& direction,
   MeshSearch const& mesh_search)
  {
    Foam::point const& center = mesh_search.mesh().cellCentres()[cell];
    Foam::scalar typ_dim = Foam::mag(center - begin);
    
    return mesh_search.tol_*typ_dim*
      direction/Foam::mag(direction);
  }
  
  // Small offset forward from begin
  // given current cell and direction
  template <typename MeshSearch>
  Foam::vector offset_forward_cell
  (Foam::point const& begin,
   Foam::label cell,
   Foam::vector const& direction,
   MeshSearch const& mesh_search)
  {
    return begin
      + offset_cell(begin, cell, direction, mesh_search);
  }
  
  // Small offset backward from begin
  // given current cell and direction
  template <typename MeshSearch>
  Foam::vector offset_backward_cell
  (Foam::point const& begin,
   Foam::label cell,
   Foam::vector const& direction,
   MeshSearch const& mesh_search)
  {
    return begin
      - offset_cell(begin, cell, direction, mesh_search);
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
