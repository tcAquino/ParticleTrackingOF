/**
* \file PTOF/Locator.h
* \author Tomás Aquino
* \date 09/03/2022
*/

#ifndef PTOF_LOCATOR_H
#define PTOF_LOCATOR_H

#include <boost/iterator/counting_iterator.hpp>
#include <fieldTypes.H>
#include <meshSearch.H>
#include "PTOF/Useful.h"

namespace ptof
{
  /** \struct Locator_Cell PTOF/Locator.h "PTOF/Locator.h"
   * \brief Object to locate positions or states in mesh. */
  struct Locator_Cell
  {
    // Typedefs for position, vector, and mesh cell index types
    using Point = Foam::point;
    using Point2D = Foam::Vector2D<Foam::scalar>;
    using Vector = Foam::vector;
    using Vector2D = Foam::Vector2D<Foam::scalar>;
    using Scalar = Foam::scalar;
    using Index = Foam::label;
    
    using MeshSearch = Foam::meshSearch; /**< Type of mesh searching object. */
    
    /** Construct given mesh searching object. */
    Locator_Cell(MeshSearch const& mesh_search)
    : mesh_search{ mesh_search }
    {}
    
    /** Find mesh cell for 3D position given hint. */
    auto operator()
    (Point const& position, Index hint = -1) const
    {
      return mesh_search.findCell(position, hint);
    }
    
    /** Find mesh cell for 2D position, given hint. */
    auto operator()
    (Point2D const& position, Index hint = -1) const
    {
      return mesh_search.findCell(make_point(position), hint);
    }
    
    /** Find mesh cell for 1D position, given hint. */
    auto operator()
    (Scalar position, Index hint = -1) const
    {
      return mesh_search.findCell(make_point(position), hint);
    }
    
    /** Find mesh cell for state, using state's cell hint. */
    template <typename State>
    auto operator()
    (State const& state) const
    {
      return mesh_search.findCell(make_point(state.position),
                                  state.cell);
    }
    
    private:
      MeshSearch const& mesh_search;  /**< Object to search mesh. */
  };
  
  /** Get all cell indexes of mesh. */
  template <typename Mesh>
  auto all_cell_ids(Mesh const& mesh)
  {
    return std::vector<Foam::label>{
      boost::counting_iterator<Foam::label>(0),
      boost::counting_iterator<Foam::label>(mesh.nCells()) };
  }
  
  /** Get face indexes in given patches of mesh. */
  template <typename Mesh>
  auto patches_face_ids
  (Mesh const& mesh, std::vector<std::string> const& patch_names)
  {
    std::vector<Foam::label> face_ids;
    for (auto const& name : patch_names)
    {
      Foam::label patch_id = mesh.boundary().findPatchID(name);
      auto const& patch = mesh.boundaryMesh()[patch_id];
      face_ids.reserve(face_ids.size()+patch.size());
      Foam::label start = patch.start();
      for (std::size_t idx = 0; idx < std::size_t(patch.size()); ++idx)
        face_ids.push_back(start+idx);
    }
    return face_ids;
  }
}


#endif /* PTOF_LOCATOR_H */
