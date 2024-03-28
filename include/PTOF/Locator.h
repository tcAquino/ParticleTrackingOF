/**
 \file PTOF/Locator.h
 \author Tomás Aquino
 \date 09/03/2022
*/

#ifndef PTOF_LOCATOR_H
#define PTOF_LOCATOR_H

#include <fieldTypes.H>
#include <meshSearch.H>
#include "PTOF/Useful.h"

namespace ptof
{
  /** \class Locator_Cell PTOF/Locator.h "PTOF/Locator.h"
   * \brief Object to locate positions or states in mesh. */
  template <typename SearchOption = SearchOptions::FirstNeighborPrecheck>
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
        if constexpr (SearchOption::neighbor_check_level >= 0)
          if (mesh_search().mesh().pointInCell(position, hint))
            return hint;

        if constexpr (SearchOption::neighbor_check_level >= 1)
          for (auto const& cell_index : mesh_search().mesh().cellCells()[hint])
            if (mesh_search().mesh().pointInCell(position, cell_index))
              return cell_index;
        
        if constexpr (SearchOption::neighbor_check_level >= 2)
          for (auto const& cell_index_1 : mesh_search().mesh().cellCells()[hint])
            for (auto const& cell_index : mesh_search().mesh().cellCells()[cell_index_1])
              if (mesh_search().mesh().pointInCell(position, cell_index))
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
}


#endif /* PTOF_LOCATOR_H */
