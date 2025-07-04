/**
   \file PTOF/Locator.h
   \author Tomas Aquino
   \date 09/03/2022
   \brief Objects to locate positions in mesh.
*/

#ifndef PTOF_LOCATOR_H
#define PTOF_LOCATOR_H

#include "PTOF/Useful.h"
#include "PTOF/SearchOptions.h"
#include <Vector2D.H>
#include <fieldTypes.H>
#include <point.H>
#include <set>

namespace ptof {
/**
   \class Locator_Cell PTOF/Locator.h "PTOF/Locator.h"
   \brief Object to locate positions or states in mesh.
*/
template <typename Geometry,
          typename SearchOption = SearchOptions::FirstNeighborPrecheck>
struct Locator_Cell {
  using Point = Foam::point;                    /**> 3D point. */
  using Point2D = Foam::Vector2D<Foam::scalar>; /**> 2D point. */
  using Scalar = Foam::scalar;                  /**> Scalar (also 1D point). */
  using Index = Foam::label;                    /**> Cell index.*/
  using Mesh = typename Geometry::Mesh;         /**< Mesh. */
  using MeshSearch =
      typename Geometry::MeshSearch; /**< Mesh searching tools. */

  /**
     \brief Constructor.
     \param geometry Domain geometry info and utilities.
  */
  Locator_Cell(Geometry const &geometry) : _geometry{geometry} {}

  Locator_Cell(Geometry &&) = delete;

  /** \return Mesh search object. */
  MeshSearch const &mesh_search() const { return _geometry.mesh_search(); }

  /** \return Mesh object. */
  Mesh const &mesh() const { return _geometry.mesh(); }

  /**
     \param position Position to locate.
     \param hint Mesh cell index hint.
     \return Mesh cell index.
  */
  auto operator()(Point const &position, Index hint = -1) const {
    if (!outside(hint)) {
      if constexpr (SearchOption::neighbor_check_level >= 0) {
        if (mesh().pointInCell(position, hint)) {
          return hint;
        }
      }
      if constexpr (SearchOption::neighbor_check_level >= 1) {
        for (auto cell_index : mesh().cellCells()[hint]) {
          if (mesh().pointInCell(position, cell_index)) {
            return cell_index;
          }
        }
      }
      if constexpr (SearchOption::neighbor_check_level >= 2) {
        std::set<Index> second_neighbors;
        for (auto cell_index_1 : mesh().cellCells()[hint]) {
          for (auto cell_index : mesh().cellCells()[cell_index_1]) {
            second_neighbors.insert(cell_index);
          }
        }
        for (auto cell_index : second_neighbors) {
          if (mesh().pointInCell(position, cell_index)) {
            return cell_index;
          }
        }
      }
    }

    auto cell_id = mesh_search().findCell(position, hint);
    // Sometimes findCell will not find a mesh cell with hint, but will without
    return outside(cell_id) && !outside(hint) ? mesh_search().findCell(position)
                                              : cell_id;
  }

  /**
     \param position Position to locate.
     \param hint Mesh cell index hint.
     \return Mesh cell index.
  */
  auto operator()(Point2D const &position, Index hint = -1) const {
    return this->operator()(make_point(position), hint);
  }

  /**
     \param position Position to locate.
     \param hint Mesh cell index hint.
     \return Mesh cell index.
  */
  auto operator()(Scalar position, Index hint = -1) const {
    return this->operator()(make_point(position), hint);
  }

  /**
     \param state State with position to locate.
     \return Mesh cell for state, using state's cell hint.
     \note State must implement:
     - \c position
     - \c cell [cell index, used as hint]
  */
  template <typename State> auto operator()(State const &state) const {
    return this->operator()(make_point(state.position), state.cell);
  }

  /**
     \param position Position to locate.
     \return Nearest mesh cell index.
     \note Does not use hint because holes are not handled correctly.
  */
  auto nearest_cell(Point const &position) const {
    return mesh_search().findNearestCell(position);
  }

  /**
     \param position Position to locate.
     \return Nearest mesh cell index.
     \note Does not use hint because holes are not handled correctly.
  */
  auto nearest_cell(Point2D const &position) const {
    return nearest_cell(make_point(position));
  }

  /**
     \param position Position to locate.
     \return Nearest mesh cell index.
     \details Does not use hint because holes are not handled correctly.
  */
  auto nearest_cell(Scalar position) const {
    return nearest_cell(make_point(position));
  }

  /**
     \param state State with position to locate.
     \return Nearest mesh cell index.
     \details Does not use hint because holes are not handled correctly.
     \note State must implement:
     - \c position
  */
  template <typename State> auto nearest_cell(State const &state) const {
    return nearest_cell(make_point(state.position));
  }

private:
  Geometry const &_geometry; /**< Domain geometry info and utilities. */
};
} // namespace ptof

#endif /* PTOF_LOCATOR_H */
