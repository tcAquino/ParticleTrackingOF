/**
* \file PTOF/Useful.h
* \author Tomás Aquino
* \date 2022/02/17
*/

#ifndef PTOF_USEFUL_H
#define PTOF_USEFUL_H

#include <type_traits>
#include <fieldTypes.H>
#include "CTRW/StateGetter.h"

namespace ptof
{
  /** Make 3D point from 2D point. */
  auto make_point(Foam::Vector2D<Foam::scalar> const& point)
  {
    return Foam::point{ point[0], point[1], 0. };
  }

  /**  Make 3D point from 1D point. */
  auto make_point(Foam::scalar point)
  {
    return Foam::point{ point, 0., 0. };
  }

  /** Make 3D point 3D point (simple copy). */
  auto make_point(Foam::point const& point)
  { return point; }
  
  /** Unit normal at boundary face, pointing outward. */
  template <typename Mesh>
  auto unit_normal_outward
  (Foam::label face, Mesh const& mesh)
  {
    return mesh.faces()[face].unitNormal(mesh.points());
  }
  
  /** Unit normal at boundary face, pointing inward. */
  template <typename Mesh>
  auto unit_normal_inward
  (Foam::label face, Mesh const& mesh)
  {
    return -mesh.faces()[face].unitNormal(mesh.points());
  }
  
  /** Area (magnitude) of a face. */
  template <typename Mesh>
  auto face_area
  (Foam::label face, Mesh const& mesh)
  {
    return mesh.faces()[face].mag(mesh.points());
  }
  
  /** Center of a face. */
  template <typename Mesh>
  auto face_center
  (Foam::label face, Mesh const& mesh)
  {
    return mesh.faces()[face].centre(mesh.points());
  }
  
  /** Volume of a cell. */
  template <typename Mesh>
  auto cell_volume
  (Foam::label cell, Mesh const& mesh)
  {
    return mesh.cellVolumes()[cell];
  }
  
  /** Center of a cell. */
  template <typename Mesh>
  auto cell_center
  (Foam::label cell, Mesh const& mesh)
  {
    return mesh.cellCentres()[cell];
  }
  
  /** Sometimes the face center associated with a mesh face
   *  is not considered within the cell,
   *  verify this */
  template <typename Mesh, typename MeshSearch>
  bool face_center_is_in_cell
  (Foam::label face, Mesh const& mesh, MeshSearch const& mesh_search)
  {
    return mesh_search.findCell(face_center(face, mesh_search.mesh()))
      == mesh.owner()[face];
  }
  
  /** Small offset forward given current face and direction. */
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
  
  /** Small offset along face normal (outward)
   * given current face
   * with face center as begin point */
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
  
  /** Small offset forward from begin
   * given current face and direction */
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
  
  /**Small offset backward from begin
   * given current face and direction */
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
  
  /** Small offset along face normal (outward)
   * given current face
   *  with face center as begin point */
  template <typename MeshSearch>
  Foam::vector offset_outward_face
  (Foam::label face,
   MeshSearch const& mesh_search)
  {
    return face_center(face, mesh_search.mesh())
      + offset_face(face, mesh_search);
  }
  
   /**Small offset along face normal (inward)
    * given current face
    * with face center as begin point */
   template <typename MeshSearch>
   Foam::vector offset_inward_face
   (Foam::label face,
    MeshSearch const& mesh_search)
   {
     return face_center(face, mesh_search.mesh())
      - offset_face(face, mesh_search);
   }
  
  /** Small offset
   * given current cell, direction, and begin point */
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
  
  /** Small offset forward
   * given current cell and direction. */
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
  
  /** Small offset backward
   * given current cell and direction. */
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
  
  /** Check for existence of periodicity info
   * (to use with State). */
  template <typename State, typename = int>
  struct Has_periodicity : std::false_type { };
  template <typename State>
  struct Has_periodicity
  <State, decltype((void) State::periodicity, 0)>
  : std::true_type { };
  
  /** Compute number of absorbed particles. */
  template <typename Subject>
  auto nr_absorbed
  (Subject const& subject, double time)
  {
    std::size_t absorbed = 0;
    for (auto const& part : subject.particles())
      
      if (part.state_new().info.absorbed
          && part.state_old().time <= time)
        ++absorbed;
    return absorbed;
  };
  
  /** Compute total mass. */
  template <typename Subject>
  auto mass
  (Subject const& subject, double time)
  {
    double mass = 0.;
    for (auto const& part : subject.particles())
      if (!part.state_new().info.absorbed
          && part.state_old().time <= time){
        mass += part.state_new().mass;
        }
    return mass;
  };
  
  /** Compute mean position (weighted by mass). */
  template
  <typename Subject, typename GetterPosition = ctrw::Get_position>
  auto position_mean
  (Subject const& subject, double time,
   GetterPosition getter_position = {})
  {
    decltype(getter_position(subject.particles(0).state_new()))
      position_mean = Foam::zero{};
    for (auto const& part : subject.particles())
    {
      auto const& state = part.state_new();
      if (!state.info.absorbed
          && part.state_old().time <= time)
      {
        position_mean += state.mass*getter_position(state);
      }
    }
    return position_mean/mass(subject, time);
  };
  
  /** Compute second moment of position (weighted by mass). */
  template
  <typename Subject, typename GetterPosition = ctrw::Get_position>
  auto position_second_moment
  (Subject const& subject, double time,
   GetterPosition getter_position = {})
  {
    decltype(getter_position(subject.particles(0).state_new()))
      second_moment = Foam::zero{};
    for (auto const& part : subject.particles())
    {
      auto const& state = part.state_new();
      if (!state.info.absorbed
          && part.state_old().time <= time)
      {
        auto position = getter_position(state);
        for (std::size_t dd = 0; dd < second_moment.size(); ++dd)
          second_moment[dd] += position[dd]*position[dd];
        second_moment *= state.mass;
      }
    }
    return second_moment/mass(subject, time);
  };
  
  /** Compute position variance (weighted by mass). */
  template
  <typename Subject, typename GetterPosition = ctrw::Get_position>
  auto position_variance
  (Subject const& subject, double time,
   GetterPosition getter_position = {})
  {
    auto position_mean_sq = position_mean(subject, time, getter_position);
    for (std::size_t dd = 0; dd < position_mean_sq.size(); ++dd)
      position_mean_sq[dd] *= position_mean_sq[dd];
    return (position_second_moment(subject, time, getter_position)-
              position_mean_sq)/mass(subject, time);
  };
}

#endif /* PTOF_USEFUL_H */
