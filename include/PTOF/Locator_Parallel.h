/**
 \file PTOF/Locator_Parallel.h
 \author Tomás Aquino
 \date 15/02/2024
*/

#ifndef PTOF_LOCATOR_PARALLEL_H
#define PTOF_LOCATOR_PARALLEL_H

#include <set>
#include <vector>
#include <omp.h>
#include <boost/iterator/counting_iterator.hpp>
#include <fieldTypes.H>
#include <meshSearch.H>
#include "PTOF/Locator.h"
#include "PTOF/Useful.h"

namespace ptof
{
  /** \struct Locator_Cell_Parallel PTOF/Locator_Parallel.h "PTOF/Locator_Parallel.h"
   * \brief Object to locate positions or states in mesh. */
  struct Locator_Cell_Parallel
  {
    using Point = Foam::point;                      /**> 3D point. */
    using Point2D = Foam::Vector2D<Foam::scalar>;   /**> 2D point. */
    using Vector = Foam::vector;                    /**> 3D vector. */
    using Vector2D = Foam::Vector2D<Foam::scalar>;  /**> 2D vector. */
    using Scalar = Foam::scalar;                    /**> Scalar (also 1D point). */
    using Index = Foam::label;                      /**> Cell index.*/
    using MeshSearch = Foam::meshSearch;            /**< Type of mesh searching object. */
    
    /** Constructor.
    \param mesh_searches Vector of unique pointers to mesh searching objects. */
    Locator_Cell_Parallel
    (std::vector<std::unique_ptr<MeshSearch>> const& mesh_searches)
    : _mesh_searches{ mesh_searches }
    {}
    
    /**
     \return Mesh searching tools for the current thread. */
    auto const& mesh_search() const
    { return *_mesh_searches[omp_get_thread_num()]; }
    
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
     \brief Mesh cell index for. */
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
    std::vector<std::unique_ptr<MeshSearch>> const& _mesh_searches;    /**< Mesh searching tools. */
  };
}


#endif /* PTOF_LOCATOR_PARALLEL_H */
