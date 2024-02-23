/**
 \file PTOF/Locator.h
 \author Tomás Aquino
 \date 09/03/2022
*/

#ifndef PTOF_LOCATOR_H
#define PTOF_LOCATOR_H

#include <set>
#include <vector>
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
    using Point = Foam::point;                      /**> 3D point. */
    using Point2D = Foam::Vector2D<Foam::scalar>;   /**> 2D point. */
    using Vector = Foam::vector;                    /**> 3D vector. */
    using Vector2D = Foam::Vector2D<Foam::scalar>;  /**> 2D vector. */
    using Scalar = Foam::scalar;                    /**> Scalar (also 1D point). */
    using Index = Foam::label;                      /**> Cell index.*/
    using MeshSearch = Foam::meshSearch;            /**< Type of mesh searching object. */
    
    /** Constructor.
     \param mesh_search Mesh searching object. */
    Locator_Cell(MeshSearch const& mesh_search)
    : _mesh_search{ mesh_search }
    {}
    
    /** \return Mesh search object. */
    MeshSearch const& mesh_search() const
    { return _mesh_search; }
    
    /**
    \param position Position to locate.
    \param hint Mesh cell index hint.
    \return Mesh cell index. */
    auto operator()
    (Point const& position, Index hint = -1) const
    {
      if (hint != -1)
      {
        auto const& mesh = mesh_search().mesh();
        if (mesh.pointInCell(position, hint))
          return hint;

        for (auto const& cell_index : mesh.cellCells()[hint])
          if (mesh.pointInCell(position, cell_index))
            return cell_index;
      }

      return mesh_search().findCell(position, hint);
    }
    
    /**
    \param position Position to locate.
    \param hint Mesh cell index hint.
    \return Mesh cell index. */
    auto operator()
    (Point2D const& position, Index hint = -1) const
    {
      return this->operator()(make_point(position), hint);
    }
    
    /**
    \param position Position to locate.
    \param hint Mesh cell index hint.
    \return Mesh cell index. */
    auto operator()
    (Scalar position, Index hint = -1) const
    {
      return this->operator()(make_point(position), hint);
    }
    
    /**
    \param state State with position to locate.
    \return Mesh cell for state, using state's cell hint.
    \note State must implement:
    - position
    - cell [cell index, used as hint] */
    template <typename State>
    auto operator()
    (State const& state) const
    {
      return this->operator()(make_point(state.position),
                              state.cell);
    }
    
    /**
     \param position Position to locate.
     \return Nearest mesh cell index.
     \details Does not use hint because holes are not handled correctly. */
    auto nearest_cell
    (Point const& position) const
    {
      return mesh_search().findNearestCell(position);
    }
    
    /**
    \param position Position to locate.
    \return Nearest mesh cell index.
    \details Does not use hint because holes are not handled correctly. */
    auto nearest_cell
    (Point2D const& position) const
    {
      return nearest_cell(make_point(position));
    }
    
    /**
    \param position Position to locate.
    \return Nearest mesh cell index.
    \details Does not use hint because holes are not handled correctly. */
    auto nearest_cell
    (Scalar position) const
    {
      return nearest_cell(make_point(position));
    }
    
    /**
     \param state State with position to locate.
     \return Nearest mesh cell index.
     \details Does not use hint because holes are not handled correctly.
     \note State must implement:
     - position */
    template <typename State>
    auto nearest_cell(State const& state) const
    {
      return nearest_cell(make_point(state.position));
    }
    
  private:
    MeshSearch const& _mesh_search;  /**< Mesh searching tools. */
  };
  
  /**
   \param mesh OpenFOAM mesh.
   \return All mesh cell indices. */
  template <typename Mesh>
  auto all_cell_ids(Mesh const& mesh)
  {
    return std::vector<Foam::label>{
      boost::counting_iterator<Foam::label>(0),
      boost::counting_iterator<Foam::label>(mesh.nCells()) };
  }
  
  /**
   \param mesh OpenFOAM mesh.
   \param patch_names Names of patches.
   \return Mesh face indices of faces given patches. */
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
  
  /**
   \param boundaries Container of pairs of lower and upper boundary locations along each dimension.
   \param mesh OpenFOAM mesh.
   \param locator Object to locate positions in mesh.
   \return Mesh cell indices within boundaries. */
  template <typename Mesh, typename Locator>
  auto cell_ids_region_cartesian
  (std::vector<std::pair<double, double>> boundaries,
   Mesh const& mesh,
   Locator const& locator)
  {
    // Identify degenerate dimensions
    // (if some minimum and maximum boundary values are equal,
    // the Cartesian region has lower dimension than the space)
    std::vector<std::size_t> degenerate_dimensions;
    std::vector<std::size_t> non_degenerate_dimensions;
    for (std::size_t dd = 0; dd < boundaries.size(); ++dd)
      if (boundaries[dd].first == boundaries[dd].second)
        degenerate_dimensions.push_back(dd);
      else
        non_degenerate_dimensions.push_back(dd);
    
    std::set<Foam::label> cell_ids;
    for (Foam::label cc = 0; cc < mesh.nCells(); ++cc)
    {
      auto center = cell_center(cc, mesh);
      
      bool cell_is_within_non_degenerate_boundaries = 1;
      for (auto dd : non_degenerate_dimensions)
      {
        if (center[dd] < boundaries[dd].first ||
            center[dd] > boundaries[dd].second)
        {
          cell_is_within_non_degenerate_boundaries = 0;
          break;
        }
      }
      if (!cell_is_within_non_degenerate_boundaries)
        continue;
      if (degenerate_dimensions.size() == 0)
        cell_ids.insert(cc);
      
      // Consider a position equal to the cell center of the candidate cell
      // but with components along the degenerate dimension
      // equal to the prescribed value. If it is in the
      // mesh, include it in the region
      for (auto dd : degenerate_dimensions)
        center[dd] = degenerate_dimensions[dd];
      auto cell_id = locator.mesh_search().findCell(center);
      if (cell_id >= 0)
        cell_ids.insert(cell_id);
    }

    std::vector<Foam::label> cells(cell_ids.begin(), cell_ids.end());
    return cells;
  }
}


#endif /* PTOF_LOCATOR_H */
