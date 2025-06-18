/**
   \file PTOF/Boundary_Cases.h
   \author Tomás Aquino
   \date 17/02/2022
   \brief Enforce different types of boundary condition depending on patch.
*/

#ifndef PTOF_BOUNDARY_CASES_H
#define PTOF_BOUNDARY_CASES_H

#include "General/IO.h"
#include "General/Meta.h"
#include "Geometry/Boundary.h"
#include "PTOF/Boundary.h"
#include "PTOF/BoundaryConditionList.h"
#include "PTOF/Reaction.h"
#include "PTOF/Useful.h"
#include <algorithm>
#include <cstddef>
#include <fieldTypes.H>
#include <iomanip>
#include <ios>
#include <map>
#include <ostream>
#include <point.H>
#include <stdexcept>
#include <string>
#include <utility>

namespace ptof {
/**
   \class Boundary_Cases PTOF/Boundary_Cases.h "PTOF/Boundary_Cases.h"
   \brief Boundary object to handle implemented boundary types.
*/
template <typename Locator, typename BoundaryInfo, typename VelocityField,
          typename Boundary_Periodic, typename Boundary_Custom,
          typename SurfaceReaction>
class Boundary_Cases {
public:
  /** \brief Container type to hold patch names and associated bc type names. */
  using BCs = std::map<std::string, std::string>;

  /**
     \brief Constructor
     \param boundary_conditions Patch names and associated boundary condition
     types.
     \param locator Object to locate positions in mesh.
     \param store_info Object to handle storing of information upon hitting a
     boundary.
     \param velocity_field Velocity field as a function of state.
     \param boundary_periodic Periodic boundary condition enforcer.
     \param boundary_custom Custom boundary condition enforcer.
     \param surface_reaction Surface Reaction.
  */
  Boundary_Cases(
      BCs boundary_conditions, Locator &&locator, BoundaryInfo &&store_info,
      VelocityField &&velocity_field,
      Boundary_Periodic &&boundary_periodic = Boundary_DoNothing{},
      Boundary_Custom &&boundary_custom = Boundary_DoNothing{},
      SurfaceReaction &&surface_reaction = SurfaceReaction_DoNothing{})
      : _boundary_conditions{boundary_conditions},
        _locator{std::forward<Locator>(locator)},
        _store_info{std::forward<BoundaryInfo>(store_info)},
        _velocity_field{std::forward<VelocityField>(velocity_field)},
        _boundary_periodic{std::forward<Boundary_Periodic>(boundary_periodic)},
        _boundary_custom{std::forward<Boundary_Custom>(boundary_custom)},
        _patch_names{_locator.mesh().boundaryMesh().names()},
        surface_reaction{std::forward<SurfaceReaction>(surface_reaction)} {
    add_unspecified_patches(_boundary_conditions, _patch_names);
    verify_boundary_conditions(_boundary_conditions, _locator.mesh(),
                               _boundary_condition_types);
  }

  /**
     \brief Check if position is out of bounds.
     \param position Position to check.
     \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    Foam::label cell = _locator(make_point(position));
    return outside(cell);
  }

  /**
     \brief Enforce boundary conditions if necessary by choosing appropriate
     types.
     \param state Current particle state to apply BCs if needed (possibly out of
     bounds)
     \param state_old Previous particle state (should be in bounds).
     \return \c true if some boundary had an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) {
    // Note: The code flow is chosen to avoid trying to locate particles which
    // are outside as much as possible, because it is expensive

    // Find first intersection with boundary patch.
    auto intersection = _locator.mesh_search().intersection(
        make_point(state_old.position), make_point(state.position));

    // When the start point is on a cell face,
    // sometimes the intersection with it is not found, and sometimes it is
    if (!intersection.hit()) {
      state.cell = _locator(state);
      if (outside(state.cell)) {
        // If the final state is outside, avoid particles leaving the domain by
        // checking for intersections with a small backwards offset
        intersection = _locator.mesh_search().intersection(
            offset_backward_cell(make_point(state_old.position), state_old.cell,
                                 make_point(state.position) -
                                     make_point(state_old.position),
                                 _locator),
            make_point(state.position));
        // If this didn't work, place particle at nearest cell center to avoid
        // breaking
        if (!intersection.hit()) {
          state.cell = _locator.nearest_cell(state.position);
          state.set_position(cell_center(state.cell, _locator.mesh()));
          return true;
        }
      }
    } else if (((make_point(state.position) - intersection.point()) &
                unit_normal_outward(intersection.index(), _locator.mesh())) <
               0.) {
      // If the displacement is against the local face normal (into the domain
      // at the boundary), ignore intersection and find
      // the next intersection with a boundary patch
      auto old_intersection = intersection;
      intersection =
          next_intersection(intersection, make_point(state.position));

      // Ignore new intersection if offset went beyond final point or if there
      // is not enough precision to find a different intersection
      if (next_intersection_is_beyond_final_point(
              intersection, old_intersection.point(), state) ||
          intersection.index() == old_intersection.index()) {
        intersection.setMiss();
      }

      // If there was no new hit and the particle is still outside, there is a
      // precision issue in finding intersections for the mesh. Place particle
      // at nearest cell center.
      if (!intersection.hit() && outside(_locator(state))) {
        state.cell = _locator.nearest_cell(state.position);
        state.set_position(cell_center(state.cell, _locator.mesh()));
        return true;
      }
    }

    bool had_effect = 0;

    // While some patch is intersected
    while (intersection.hit()) {
      auto patch_id = ptof::patch_id(intersection.index(), _locator.mesh());
      auto patch = patch_name(patch_id);
      auto type_name = _boundary_conditions.at(patch);

      // Apply the approriate boundary and store associated info
      switch (_boundary_condition_types.type(type_name)) {
      case BoundaryConditionList::Type::reflecting: {
        _store_info(state, state_old, intersection, _boundary_condition_types,
                    meta::Selector<BoundaryConditionList::Type,
                                   BoundaryConditionList::Type::reflecting>{});
        boundary_reflecting(state, intersection.point(),
                            reflection_normal(intersection.index()));
        had_effect = 1;
        break;
      }
      case BoundaryConditionList::Type::reacting_reflecting: {
        _store_info(
            state, state_old, intersection, _boundary_condition_types,
            meta::Selector<BoundaryConditionList::Type,
                           BoundaryConditionList::Type::reacting_reflecting>{});
        surface_reaction(state, state_old, intersection.index());
        boundary_reflecting(state, intersection.point(),
                            reflection_normal(intersection.index()));
        had_effect = 1;
        break;
      }
      case BoundaryConditionList::Type::periodic: {
        _store_info(state, state_old, intersection, _boundary_condition_types,
                    meta::Selector<BoundaryConditionList::Type,
                                   BoundaryConditionList::Type::periodic>{});
        had_effect += _boundary_periodic(state, intersection);
        break;
      }
      case BoundaryConditionList::Type::absorbing: {
        _store_info(state, state_old, intersection, _boundary_condition_types,
                    meta::Selector<BoundaryConditionList::Type,
                                   BoundaryConditionList::Type::absorbing>{});
        state.info.absorbed = true;
        state.set_position(intersection.point());
        had_effect = 1;
        break;
      }
      case BoundaryConditionList::Type::inlet: {
        if (face_flux_outward(intersection.index(), _velocity_field, _locator) >
            0.) {
          _store_info(state, state_old, intersection, _boundary_condition_types,
                      meta::Selector<BoundaryConditionList::Type,
                                     BoundaryConditionList::Type::absorbing>{});
          state.set_position(intersection.point());
          state.cell = _locator(state);
        } else {
          _store_info(
              state, state_old, intersection, _boundary_condition_types,
              meta::Selector<BoundaryConditionList::Type,
                             BoundaryConditionList::Type::reflecting>{});
          boundary_reflecting(state, intersection.point(),
                              reflection_normal(intersection.index()));
        }
        had_effect = 1;
        break;
      }
      case BoundaryConditionList::Type::custom: {
        _store_info(state, state_old, intersection, _boundary_condition_types,
                    meta::Selector<BoundaryConditionList::Type,
                                   BoundaryConditionList::Type::custom>{});
        had_effect += _boundary_custom(state, state_old, intersection);
        break;
      }
      case BoundaryConditionList::Type::empty: {
      }
      default:
        throw std::runtime_error{"Boundary condition type " + type_name +
                                 " not supported"};
      }

      // Stop if right at the last intersection point
      // This avoids numerical issues for some edge cases
      if (make_point(state.position) == intersection.point()) {
        state.cell = _locator(state);
        // For some meshes, some intersections can be considered outside the
        // mesh. If so, use the owner cell center instead.
        if (outside(state.cell)) {
          auto owner_cell = _locator.mesh().faceOwner()[intersection.index()];
          state.set_position(cell_center(owner_cell, _locator.mesh()),
                             owner_cell);
          return true;
        }
        break;
      }

      // Find the next intersection with a boundary patch
      auto old_intersection = intersection;
      intersection =
          next_intersection(intersection, make_point(state.position));

      // Ignore new intersection if offset went beyond final point or if there
      // is not enough precision to find a different intersection
      if (next_intersection_is_beyond_final_point(
              intersection, old_intersection.point(), state) ||
          intersection.index() == old_intersection.index()) {
        intersection.setMiss();
      }

      if (!intersection.hit()) {
        state.cell = _locator(state);
        // Sometimes there is not enough precision to find the next
        // intersection. Avoid particles leaving the domain by placing the final
        // state at the last intersection.
        if (outside(state.cell)) {
          state.set_position(old_intersection.point());
          state.cell = _locator(state);
          // For some meshes, some intersections can be considered outside the
          // mesh. If so, use the owner cell center instead.
          if (outside(state.cell)) {
            auto owner_cell =
                _locator.mesh().faceOwner()[old_intersection.index()];
            state.set_position(cell_center(owner_cell, _locator.mesh()),
                               owner_cell);
            return true;
          }
        }
        break;
      }
    }

    return had_effect;
  }

  /**
     \brief Output information about current object.
     \param output Output stream.
  */
  std::ostream &info_runtime(std::ostream &output) const {
    io::StreamScopeFormat guard{output};
    output << io::line() << "Boundary conditions\n" << io::line();
    if (_boundary_conditions.empty()) {
      output << "None\n" << io::line();
      return output;
    }
    std::size_t width_patch =
        std::max_element(_boundary_conditions.begin(),
                         _boundary_conditions.end(),
                         [](auto const &aa, auto const &bb) {
                           return aa.first.length() < bb.first.length();
                         })
            ->first.length() +
        2;
    std::size_t width_bc =
        std::max_element(_boundary_conditions.begin(),
                         _boundary_conditions.end(),
                         [](auto const &aa, auto const &bb) {
                           return aa.second.length() < bb.second.length();
                         })
            ->second.length() +
        2;
    width_bc = std::max(width_bc, _boundary_custom.name().length());
    output << std::left;
    for (auto const &bc : _boundary_conditions) {
      output << std::setw(width_patch) << bc.first;
      if (bc.second == "custom")
        output << std::setw(width_bc) << _boundary_custom.name();
      else
        output << std::setw(width_bc) << bc.second;
      output << "\n";
    }
    output << io::line();
    return output;
  }

private:
  BCs _boundary_conditions; /**< Patch names and associated bc type names. */
  Locator _locator;         /**< Object to locate positions in mesh. */
  BoundaryInfo _store_info; /**< Object to handle boundary-related info storing
                               in states. */
  VelocityField _velocity_field;
  Boundary_Periodic
      _boundary_periodic; /**< Boundary object to handle 'periodic' bc type. */
  Boundary_Custom
      _boundary_custom; /**< Boundary object to handle 'custom' bc type. */
  Foam::wordList _patch_names;
  ; /**< Names of patches in mesh. */

public:
  SurfaceReaction
      surface_reaction; /**< Surface reaction for reactive boundaries*/

private:
  const BoundaryConditionList
      _boundary_condition_types{}; /**< Boundary condition names and types. */

  /**
     \param patch Index of patch in mesh.
     \return Name of patch.
  */
  auto patch_name(std::size_t patch) const { return _patch_names[patch]; }

  /**
     \param face Face index.
     \return Inward normal to face (used for reflection).
  */
  auto reflection_normal(Foam::label face) const {
    return unit_normal_inward(face, _locator.mesh());
  }

  /**
     \brief Find next intersection with a boundary.
     \details Find the next intersection, if any, along straight line from
     previous intersection to endpoint.
     \param intersection Previous intersection.
     \param end Endpoint.
     \return Next intersection.
  */
  template <typename Intersection>
  auto next_intersection(Intersection const &intersection,
                         Foam::point const &end) const {
    return _locator.mesh_search().intersection(
        offset_forward_face_keep_inside(intersection.point(),
                                        intersection.index(),
                                        end - intersection.point(), _locator),
        end);
  }

  /**
     \brief Check if next intersection is beyond final position.
     \details Find the next intersection, if any, along straight line from
     previous intersection to endpoint.
     \param next_intersection Next intersection.
     \param old_intersection_point Position of previous intersection.
     \param state Final state before new intersection.
     \return \c true if beyond final position, \c false otherwise.
  */
  template <typename Intersection, typename State>
  auto next_intersection_is_beyond_final_point(
      Intersection const &next_intersection,
      Foam::point const &old_intersection_point, State const &state) {
    return next_intersection.hit() &&
           ((next_intersection.point() - old_intersection_point) &
            (make_point(state.position) - next_intersection.point())) < 0.;
  }
};
template <typename Locator, typename BoundaryInfo, typename VelocityField,
          typename Boundary_Periodic, typename Boundary_Custom,
          typename Surface_Reaction>
Boundary_Cases(typename Boundary_Cases<Locator, BoundaryInfo, VelocityField,
                                       Boundary_Periodic, Boundary_Custom,
                                       Surface_Reaction>::BCs,
               Locator &&, BoundaryInfo &&, VelocityField &&,
               Boundary_Periodic &&, Boundary_Custom &&, Surface_Reaction &&)
    -> Boundary_Cases<Locator, BoundaryInfo, VelocityField, Boundary_Periodic,
                      Boundary_Custom, Surface_Reaction>;
template <typename Locator, typename BoundaryInfo, typename VelocityField,
          typename Boundary_Periodic, typename Boundary_Custom>
Boundary_Cases(typename Boundary_Cases<Locator, BoundaryInfo, VelocityField,
                                       Boundary_Periodic, Boundary_Custom,
                                       SurfaceReaction_DoNothing>::BCs,
               Locator &&, BoundaryInfo &&, VelocityField &&,
               Boundary_Periodic &&, Boundary_Custom &&)
    -> Boundary_Cases<Locator, BoundaryInfo, VelocityField, Boundary_Periodic,
                      Boundary_Custom, SurfaceReaction_DoNothing>;
template <typename Locator, typename BoundaryInfo, typename VelocityField,
          typename Boundary_Periodic>
Boundary_Cases(typename Boundary_Cases<Locator, BoundaryInfo, VelocityField,
                                       Boundary_Periodic, Boundary_DoNothing,
                                       SurfaceReaction_DoNothing>::BCs,
               Locator &&, BoundaryInfo &&, VelocityField &&,
               Boundary_Periodic &&)
    -> Boundary_Cases<Locator, BoundaryInfo, VelocityField, Boundary_Periodic,
                      Boundary_DoNothing, SurfaceReaction_DoNothing>;
template <typename Locator, typename BoundaryInfo, typename VelocityField>
Boundary_Cases(typename Boundary_Cases<Locator, BoundaryInfo, VelocityField,
                                       Boundary_DoNothing, Boundary_DoNothing,
                                       SurfaceReaction_DoNothing>::BCs,
               Locator &&, BoundaryInfo &&, VelocityField &&)
    -> Boundary_Cases<Locator, BoundaryInfo, VelocityField, Boundary_DoNothing,
                      Boundary_DoNothing, SurfaceReaction_DoNothing>;
} // namespace ptof

#endif /* PTOF_BOUNDARYCONDITION_CASES_H */
