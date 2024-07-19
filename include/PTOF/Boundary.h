/**
\file PTOF/Boundary.h
\author Tomás Aquino
\date 17/02/2022
*/

#ifndef PTOF_BOUNDARY_H
#define PTOF_BOUNDARY_H

#include "General/Meta.h"
#include "General/Operations.h"
#include "General/Useful.h"
#include "PTOF/Directories.h"
#include "PTOF/Reaction.h"
#include "PTOF/Useful.h"
#include <algorithm>
#include <cstddef>
#include <fieldTypes.H>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <point.H>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ptof {
/** \struct BoundaryCondition PTOF/Boundary.h "PTOF/Boundary.h"
 \brief Names and types of boundary conditions. */
struct BoundaryCondition {
  /** \enum Type
   *  \brief Implemented types. */
  enum class Type {
    reflecting, /**< Reflecting                                         */
    reacting_reflecting, /**< Reacting and reflecting                   */
    periodic,  /**< Periodic                                            */
    absorbing, /**< Absorbing                                           */
    info,      /**< Information upon hitting                            */
    custom,    /**< Enforced by custom Boundary object                  */
    empty      /**< No effect (default for unspecified patches in mesh) */
  };

  /** Constructor.
   \brief Default constructor. */
  BoundaryCondition() {}

  /** \brief Type from name.
   \param name Boundary condition name.
   \return Boundary condition type.
   */
  static auto type(std::string const &name) { return name_to_type.at(name); }

  /** \brief Name from type.
   \param type Boundary condition type.
   \return Boundary condition name.
  */
  static auto name(Type type) { return type_to_name.at(type); }

  /** \brief Check if name exists.
   \param name condition name.
   \return \c true if name exists, \c false otherwise.
   */
  static bool contains(std::string const &name) {
    return name_to_type.count(name);
  }

  /** \brief Map names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"reflecting", Type::reflecting},
      {"reacting_reflecting", Type::reacting_reflecting},
      {"periodic", Type::periodic},
      {"absorbing", Type::absorbing},
      {"custom", Type::custom},
      {"empty", Type::empty},
  };

  /** \brief Map types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::reflecting, "reflecting"},
      {Type::reacting_reflecting, "reacting_reflecting"},
      {Type::periodic, "periodic"},
      {Type::absorbing, "absorbing"},
      {Type::custom, "custom"},
      {Type::empty, "empty"}};
};

/**
 \brief Apply reflecting boundary condition.
 \param state Particle state to update.
 \param contact_point Contact point at boundary.
 \param unit_normal Unit normal at the contact point, pointing inwards. */
template <typename State>
void boundary_reflecting(State &state, Foam::point const &contact_point,
                         Foam::vector const &unit_normal) {
  auto leftover_jump = make_point(state.position) - contact_point;

  // When a particle starts at a reflecting wall,
  // it must not be reflected if it jumps towards the interior
  // and not through the boundary
  // This dot product can only be positive in that case
  auto leftover_dot_normal = leftover_jump & unit_normal;
  state.set_position(contact_point + leftover_jump +
                     (leftover_dot_normal < 0.
                          ? -2. * leftover_dot_normal * unit_normal
                          : Foam::zero{}));
}

/**
\brief Apply absorbing boundary condition.
\param state Particle state to update.
\param contact_point Contact point at boundary. */
template <typename State>
void boundary_absorbing(State &state, Foam::point const &contact_point) {
  state.set_position(contact_point);
}

/**
 \brief Get the boundary conditions associated with a set of patches from file.
 \param filename Name of file holding a patch name followed by a boundary
 condition type in each line.
 \param delims Possible delimiters between patch
 names and boundary condition types.
*/
inline auto get_boundary_conditions(std::string const &filename,
                                    std::string const &delims = "\t,| ") {
  std::map<std::string, std::string> boundary_conditions;
  std::string comment_sequence = "#";

  auto file = useful::open_read(filename);
  std::string line;
  while (std::getline(file, line)) {
    if (line.find(comment_sequence) == 0)
      continue;
    useful::remove_carriage_return_in_place(line);
    useful::clear_escape_in_place(line, comment_sequence);
    auto split_line = useful::split(line, "\t,| ");
    if (split_line.size() != 2)
      throw useful::parse_error(filename, line);
    boundary_conditions[split_line[0]] = split_line[1];
  }
  file.close();

  return boundary_conditions;
}

/**
 \brief Get boundary conditions for dynamics type.
 \param directories Current case directory information.
 \param delims Possible delimiters between patch names and boundary condition
 types.
 */
template <Dynamics::Type dynamics>
auto get_boundary_conditions(Directories const &directories,
                             std::string const &delims = "\t,| ") {
  if constexpr (dynamics == Dynamics::Type::transport)
    return get_boundary_conditions(directories.dir_boundaryconditions +
                                       "/boundary_conditions_" +
                                       Dynamics::name(dynamics) + ".dat",
                                   delims);
  if constexpr (dynamics == Dynamics::Type::firstpassage)
    return get_boundary_conditions(directories.dir_boundaryconditions +
                                       "/boundary_conditions_" +
                                       Dynamics::name(dynamics) + ".dat",
                                   delims);
}

/** \brief Throw if any patch name does not exist in mesh or associated boundary
 condition type is not implemented.
 \param boundary_conditions Container of
 pairs of patch names and boundary condition types.
 \param mesh Mesh object.
 \param implemented BoundaryCondition object holding information about
 implemented boundary conditions. */
template <typename BCs, typename Mesh, typename Implemented>
void verify_boundary_conditions(BCs const &boundary_conditions,
                                Mesh const &mesh,
                                Implemented const &implemented) {
  for (auto const &bc : boundary_conditions)
    mesh.boundaryMesh().findPatchID(bc.first, false);

  for (auto const &bc : boundary_conditions)
    if (!implemented.contains(bc.second))
      throw std::runtime_error{"Boundary condition type " + bc.second +
                               " not supported"};
}

/**
 \brief Add default type 'empty' to boundary conditions for unspecified patches.
 \param boundary_conditions Associative container of pairs of patch names and
 boundary condition types, to add to.
 \param patch_names Container with all
 patch names.
*/
template <typename BCs, typename Container>
void add_unspecified_patches(BCs &boundary_conditions,
                             Container const &patch_names) {
  for (auto const &name : patch_names)
    if (!boundary_conditions.count(name))
      boundary_conditions[name] = "empty";
}

/**
 \brief Add default type 'empty' to boundary conditions for unspecified patches
 in mesh.
 \param boundary_conditions Associative container of pairs of patch
 names and boundary condition types, to add to.
 \param mesh Mesh object.
*/
template <typename BCs, typename Mesh>
void add_unspecified_patches_in_mesh(BCs &boundary_conditions,
                                     Mesh const &mesh) {
  add_unspecified_patches(boundary_conditions, mesh.boundaryMesh().names());
}

/** \struct Boundary_DoNothing PTOF/Boundary.h "PTOF/Boundary.h"
 * \brief Boundary object with no effect. */
struct Boundary_DoNothing {
  /**
   \brief Check if position is out of bounds.
   \param position Position to check.
   \return \c false (always in bounds).
   */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return false;
  }

  /**
   \brief Apply boundary condition (do nothing).
   \param state Current particle state.
   \param state_old Previous particle state.
   \param dummy Dummy argument (unused).
   \return \c false (no effect).
  */
  template <typename State, typename Dummy>
  bool operator()(State const &state, State const &state_old,
                  Dummy const &dummy) const {
    return false;
  }

  /**
   \brief Apply boundary condition (do nothing).
   \param state Current particle state.
   \param dummy Dummy argument (unused).
   \return \c false (no effect).
  */
  template <typename State, typename Dummy>
  bool operator()(State const &state, Dummy const &dummy) const {
    return false;
  }

  /**
   \brief Apply boundary condition (do nothing).
   \param state Current particle state.
   \return \c false (no effect).
  */
  template <typename State> bool operator()(State const &state) const {
    return false;
  }

  /** \brief Name of boundary condition type. */
  std::string name() const { return "empty"; }
};

/** \class Boundary_Cases PTOF/Boundary.h "PTOF/Boundary.h"
 *  \brief Boundary object to handle implemented boundary types. */
template <typename Locator, typename Store_Info, typename Boundary_Periodic,
          typename Boundary_Custom, typename SurfaceReaction>
class Boundary_Cases {
public:
  /** \brief Container type to hold patch names and associated bc type names. */
  using BCs = std::map<std::string, std::string>;

  /** Constructor
   \param boundary_conditions Patch names and associated  boundary condition
   types.
   \param locator Object to locate positions in mesh.
   \param store_info Object to handle storing of information upon boundary
   hitting. \param boundary_periodic Periodic boundary condition enforcer.
   \param boundary_custom Custom boundary condition enforcer.
   \param surface_reaction
   Surface Reaction.
  */
  Boundary_Cases(
      BCs boundary_conditions, Locator &&locator, Store_Info &&store_info,
      Boundary_Periodic &&boundary_periodic = Boundary_DoNothing{},
      Boundary_Custom &&boundary_custom = Boundary_DoNothing{},
      SurfaceReaction &&surface_reaction = SurfaceReaction_DoNothing{})
      : _boundary_conditions{boundary_conditions},
        _locator{std::forward<Locator>(locator)},
        _store_info{std::forward<Store_Info>(store_info)},
        _boundary_periodic{std::forward<Boundary_Periodic>(boundary_periodic)},
        _boundary_custom{std::forward<Boundary_Custom>(boundary_custom)},
        surface_reaction{std::forward<SurfaceReaction>(surface_reaction)} {
    add_unspecified_patches(_boundary_conditions, patch_names());
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
    Foam::label cell = _locator.mesh_search().findCell(make_point(position));
    return outside(cell);
  }

  /**
   \brief Enforce boundary conditions if necessary by choosing appropriate
   types.
   \param state Current particle state to apply BCs if needed (possibly
   out of bounds).
   \param state_old Previous particle state (should be in bounds).
   \return \c true if some boundary had an effect, \c false otherwise.
  */
  template <typename State>
  bool operator()(State &state, State const &state_old = {}) {
    // Find first intersection with boundary patch.
    auto intersection = _locator.mesh_search().intersection(
        make_point(state_old.position), make_point(state.position));

    // When the start point is on a cell face,
    // sometimes the intersection with it is not found, and sometimes it is
    if (!intersection.hit()) {
      // Note: The code flow is chosen to avoid trying to locate particles which
      // are outside as much as possible, because it is expensive
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
        if (!intersection.hit())
          state.set_position(cell_center(_locator.nearest_cell(state.position),
                                         _locator.mesh()));
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
          intersection.index() == old_intersection.index())
        intersection.setMiss();

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
      /** Find patch in mesh and associated boundary type name. */
      auto patch_id =
          _locator.mesh().boundaryMesh().whichPatch(intersection.index());
      auto patch = patch_name(patch_id);
      auto type_name = _boundary_conditions.at(patch);

      // Apply the approriate boundary and store associated info
      switch (_boundary_condition_types.type(type_name)) {
      case BoundaryCondition::Type::reflecting: {
        _store_info(state, state_old, intersection, _boundary_condition_types,
                    meta::Selector<BoundaryCondition::Type,
                                   BoundaryCondition::Type::reflecting>{});
        boundary_reflecting(state, intersection.point(),
                            reflection_normal(intersection.index()));
        had_effect = 1;
        break;
      }
      case BoundaryCondition::Type::reacting_reflecting: {
        _store_info(state, state_old, intersection, _boundary_condition_types,
                    meta::Selector<BoundaryCondition::Type,
                                   BoundaryCondition::Type::reflecting>{});
        surface_reaction(state, state_old, intersection.index());
        boundary_reflecting(state, intersection.point(),
                            reflection_normal(intersection.index()));
        had_effect = 1;
        break;
      }
      case BoundaryCondition::Type::periodic: {
        _store_info(state, state_old, intersection, _boundary_condition_types,
                    meta::Selector<BoundaryCondition::Type,
                                   BoundaryCondition::Type::periodic>{});
        had_effect += _boundary_periodic(state, intersection);
        break;
      }
      case BoundaryCondition::Type::absorbing: {
        _store_info(state, state_old, intersection, _boundary_condition_types,
                    meta::Selector<BoundaryCondition::Type,
                                   BoundaryCondition::Type::absorbing>{});
        boundary_absorbing(state, intersection.point());
        had_effect = 1;
        break;
      }
      case BoundaryCondition::Type::custom: {
        _store_info(state, state_old, intersection, _boundary_condition_types,
                    meta::Selector<BoundaryCondition::Type,
                                   BoundaryCondition::Type::custom>{});
        had_effect += _boundary_custom(state, state_old, intersection);
        break;
      }
      case BoundaryCondition::Type::empty: {
      }
      default:
        throw std::runtime_error{"Boundary condition type " + type_name +
                                 " not supported"};
      }

      // Stop if right at the last intersection point
      // This avoids numerical issues for some edge cases
      if (make_point(state.position) == intersection.point())
        break;

      // Find the next intersection with a boundary patch
      auto old_intersection = intersection;
      intersection = next_intersection(intersection, make_point(state.position));
      // Ignore new intersection if offset went beyond final point or if there
      // is not enough precision to find a different intersection
      if (next_intersection_is_beyond_final_point(
              intersection, old_intersection.point(), state) ||
          intersection.index() == old_intersection.index())
        intersection.setMiss();

      // Sometimes there is not enough precision to find the next intersection.
      // Avoid particles leaving the domain by placing the final state at the
      // last intersection
      if (!intersection.hit()) {
        state.cell = _locator(state);
        if (outside(state.cell))
          state.set_position(old_intersection.point());
        break;
      }
    }

    return had_effect;
  }

  /** \brief Output information about current object. */
  template <typename OStream> void info_runtime(OStream &output) const {
    output
        << "--------------------------------------------------------------\n"
           "Boundary conditions\n"
           "--------------------------------------------------------------\n";
    if (_boundary_conditions.empty()) {
      output
          << "None\n"
          << "--------------------------------------------------------------\n";
      return;
    }
    int width_patch =
        int(std::max_element(_boundary_conditions.begin(),
                             _boundary_conditions.end(),
                             [](auto const &aa, auto const &bb) {
                               return aa.first.length() < bb.first.length();
                             })
                ->first.length()) +
        1;
    int width_bc =
        int(std::max_element(_boundary_conditions.begin(),
                             _boundary_conditions.end(),
                             [](auto const &aa, auto const &bb) {
                               return aa.second.length() < bb.second.length();
                             })
                ->second.length()) +
        1;
    width_bc = std::max(width_bc, int(_boundary_custom.name().length()));
    for (auto const &bc : _boundary_conditions) {
      output << std::left << std::setw(width_patch) << bc.first;
      if (bc.second == "custom")
        output << std::left << std::setw(width_bc) << _boundary_custom.name();
      else
        output << std::left << std::setw(width_bc) << bc.second;
      output << "\n";
    }
    output
        << "--------------------------------------------------------------\n";
  }

private:
  BCs _boundary_conditions; /**< Patch names and associated bc type names. */
  Locator _locator;         /**< Object to locate positions in mesh. */
  Store_Info
      _store_info; /**< Object to handle boundary info storing in states. */
  Boundary_Periodic
      _boundary_periodic; /**< Boundary object to handle 'periodic' bc type. */
  Boundary_Custom
      _boundary_custom; /**< Boundary object to handle 'custom' bc type. */

public:
  SurfaceReaction
      surface_reaction; /**< Surface reaction for reactive boundaries*/

private:
  const BoundaryCondition
      _boundary_condition_types{}; /**< Boundary condition names and types. */

  /** \return Names of patches in mesh. */
  auto patch_names() const { return _locator.mesh().boundaryMesh().names(); }

  /**
   \param patch Index of patch in mesh.
   \return Name of patch. */
  auto patch_name(std::size_t patch) const { return patch_names()[patch]; }

  /**
   \param face Face index.
   \return Inward normal to face (used for reflection). */
  auto reflection_normal(Foam::label face) const {
    return unit_normal_inward(face, _locator.mesh());
  }

  /** \brief Find next intersection with a boundary.
   \details Find the next intersection, if any, along straight line from
   previous intersection to endpoint.
   \param intersection previous_intersection.
   \param end Endpoint.
   \param cell Cell index for previous intersection.
   \return Next intersection. */
  template <typename Intersection>
  auto next_intersection(Intersection const &intersection,
                         Foam::point const &end) const {
    return _locator.mesh_search().intersection(
        offset_forward_face_keep_inside(intersection.point(),
                                        intersection.index(),
                                        end - intersection.point(), _locator),
        end);
  }

  template <typename Intersection, typename State>
  auto next_intersection_is_beyond_final_point(
      Intersection const &next_intersection,
      Foam::point const &old_intersection_point, State const &state) {
    return next_intersection.hit() &&
           ((next_intersection.point() - old_intersection_point) &
            (make_point(state.position) - next_intersection.point())) < 0.;
  }
};
template <typename Locator, typename Store_Info, typename Boundary_Periodic,
          typename Boundary_Custom, typename Surface_Reaction>
Boundary_Cases(typename Boundary_Cases<Locator, Store_Info, Boundary_Periodic,
                                       Boundary_Custom, Surface_Reaction>::BCs,
               Locator &&, Store_Info &&, Boundary_Periodic &&,
               Boundary_Custom &&, Surface_Reaction &&)
    -> Boundary_Cases<Locator, Store_Info, Boundary_Periodic, Boundary_Custom,
                      Surface_Reaction>;
template <typename Locator, typename Store_Info, typename Boundary_Periodic,
          typename Boundary_Custom>
Boundary_Cases(
    typename Boundary_Cases<Locator, Store_Info, Boundary_Periodic,
                            Boundary_Custom, SurfaceReaction_DoNothing>::BCs,
    Locator &&, Store_Info &&, Boundary_Periodic &&, Boundary_Custom &&)
    -> Boundary_Cases<Locator, Store_Info, Boundary_Periodic, Boundary_Custom,
                      SurfaceReaction_DoNothing>;
template <typename Locator, typename Store_Info, typename Boundary_Periodic>
Boundary_Cases(
    typename Boundary_Cases<Locator, Store_Info, Boundary_Periodic,
                            Boundary_DoNothing, SurfaceReaction_DoNothing>::BCs,
    Locator &&, Store_Info &&, Boundary_Periodic &&)
    -> Boundary_Cases<Locator, Store_Info, Boundary_Periodic,
                      Boundary_DoNothing, SurfaceReaction_DoNothing>;
template <typename Locator, typename Store_Info>
Boundary_Cases(
    typename Boundary_Cases<Locator, Store_Info, Boundary_DoNothing,
                            Boundary_DoNothing, SurfaceReaction_DoNothing>::BCs,
    Locator &&, Store_Info &&)
    -> Boundary_Cases<Locator, Store_Info, Boundary_DoNothing,
                      Boundary_DoNothing, SurfaceReaction_DoNothing>;

/** \brief Find periodic image of intersection point.
 \param state_outside Out-of-bounds state.
 \param state_image Periodic image of state, in bounds.
 \param intersection Boundary intersection info.
 \param boundary Periodic boundary enforcer.
 \param locator Object to locate positions in mesh.
 \return Periodic image of intersection. */
template <typename State, typename Intersection, typename Boundary,
          typename Locator>
auto periodic_intersection(State state_outside, State const &state_image,
                           Intersection const &intersection,
                           Boundary const &boundary, Locator const &locator) {
  // Compute displacement from intersection to outside position.
  auto displacement = make_point(state_outside.position) - intersection.point();

  // Slight offset from intersection towards the outside point.
  // The goal is to find the periodic image of a point close to the
  // intersection, so we can find its periodic image within the domain and then
  // use it to find the periodic image of the intersection itself.
  auto offset_pos = offset_forward_face(
      intersection.point(), intersection.index(), displacement, locator);
  auto offset = offset_pos - intersection.point();

  // Ensure offset position is closer to intersection than outside position.
  // Otherwise, use displacement.
  if (Foam::magSqr(offset) < Foam::magSqr(displacement))
    state_outside.set_position(offset_pos);
  else {
    offset = make_point(state_outside.position) - intersection.point();
    state_outside.set_position(intersection.point() + offset);
  }

  // Place the offset point at its periodic image
  boundary(state_outside);

  // Find the periodic image of the intersection by looking from outside the
  // domain towards the offset point

  // For some meshes, numerical error can lead to the intersection point not
  // being found. Avoid breaking by repeatedly doubling the size of the interval
  // until the intersection is found or too many tries have been made.
  std::size_t max_nr_doublings = 20;
  std::size_t iter = 0;
  do {
    offset *= 2;
    auto new_intersection = locator.mesh_search().intersection(
        make_point(state_outside.position) - offset,
        make_point(state_outside.position));
    if (new_intersection.hit())
      return new_intersection;
  } while (iter++ < max_nr_doublings);

  return Intersection{};
}

/** \class Boundary_Periodic_OF PTOF/Boundary.h "PTOF/Boundary.h"
 *  \brief Interface for periodic boundary enforcer. */
template <typename Boundary, typename Locator> struct Boundary_Periodic_OF {
  /** Constructor
   \param boundary Periodic boundary enforcer.
   \param locator Object to locate positions in mesh.
  */
  Boundary_Periodic_OF(Boundary &&boundary, Locator &&locator)
      : boundary_periodic{std::forward<Boundary>(boundary)},
        locator{std::forward<Locator>(locator)} {}

  /**
   \brief Check if position is out of bounds.
   \param position Position to check
   \return \c true if out of bounds, \c false otherwise.
  */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return boundary_periodic.out_of_bounds(position);
  }

  /** \brief Enforce boundary condition and place intersection at its periodic
   image.
   \param state Particle state to apply BC to.
   \param intersection Boundary intersection info, to be replaced by periodic
   image.
   \return \c true
   */
  template <typename State, typename Intersection>
  bool operator()(State &state, Intersection &intersection) const {
    auto state_outside = state;
    boundary_periodic(state);
    intersection = periodic_intersection(state_outside, state, intersection,
                                         boundary_periodic, locator);
    if (!intersection.hit())
      throw std::runtime_error{
        "Could not find image of periodic boundary intersection"};

    return true;
  }

  Boundary boundary_periodic; /**< Periodic boundary type. */
  Locator locator;            /**< Object to locate positions in mesh. */
};
template <typename Boundary, typename Locator>
Boundary_Periodic_OF(Boundary &&, Locator &&)
    -> Boundary_Periodic_OF<Boundary, Locator>;

/** \class Boundary_Reinject PTOF/Boundary.h "PTOF/Boundary.h"
 \brief Boundary object that reinjects particles according to prescribed initial
 condition. */
template <typename InitialCondition> struct Boundary_Reinject {
  /** Constructor.
   \param initial_condition Initial condition object. */
  Boundary_Reinject(InitialCondition &&initial_condition)
      : initial_condition{std::forward<InitialCondition>(initial_condition)} {}

  /** \brief Enforce boundary condition.
  \param state Particle state to apply BC to.
  \param state_old Old particle state (unused).
  \param intersection Boundary intersection info (unused).
  \return \c true
  */
  template <typename State, typename Intersection>
  bool operator()(State &state, State const &state_old,
                  Intersection const &intersection) const {
    state.set_position(initial_condition.make_position());
    return true;
  }

  /** \brief Name of boundary condition type. */
  std::string name() const { return "reinject"; }

  InitialCondition initial_condition; /**< Initial condition according to which
                                         to reinject. */
};
template <typename InitialCondition>
Boundary_Reinject(InitialCondition &&) -> Boundary_Reinject<InitialCondition>;

/** \brief For boundaries conditions indicated as type \c periodic, extract the
 corresponding boundary positions.
 \param boundary_conditions Associative
 container of patch names and associated boundary condition types.
 \param mesh
 Mesh object. \param dim Spatial dimensionality.
 \return Container of pairs of periodic boundary positions along
 each periodic dimension.
 \note
 - Periodic boundaries are assumed to be perpendicular to the cartesian axes.
 There must be an even number.
 - Periodic conditions are checked in order across Cartesian dimensions. If
 there are periodic BCs in all dimensions starting from 0 up to some dimension,
 output only those boundaries. If there are gap dimensions between periodic
 dimensions, add {-inf, inf} boundaries.*/
template <typename BCs, typename Mesh>
auto extract_cartesian_periodic_boundaries(BCs const &boundary_conditions,
                                           Mesh const &mesh, std::size_t dim) {
  std::vector<std::vector<double>> boundaries_dim(dim);
  for (auto const &bc : boundary_conditions)
    if (bc.second == "periodic") {
      // Find average and variance of patch face centers.
      auto const &patch_id = mesh.boundaryMesh().findPatchID(bc.first, false);
      auto const &patch = mesh.boundaryMesh()[patch_id];
      auto face_start = patch.start();
      Foam::point average_face_center = Foam::point::zero;
      Foam::point variance_face_center = Foam::point::zero;
      Foam::scalar total_area = 0.;
      for (auto face = face_start; face < face_start + patch.size(); ++face) {
        auto area = face_area(face, mesh);
        auto center = face_center(face, mesh);
        auto weighted_center = area * center;
        total_area += area;
        average_face_center += weighted_center;
        for (std::size_t dd = 0; dd < variance_face_center.size(); ++dd)
          variance_face_center[dd] += weighted_center[dd] * center[dd];
      }
      average_face_center /= total_area;
      variance_face_center /= total_area;
      for (std::size_t dd = 0; dd < variance_face_center.size(); ++dd)
        variance_face_center[dd] -=
            average_face_center[dd] * average_face_center[dd];

      // Choose dimension with least variance.
      auto it_min = std::min_element(variance_face_center.begin(),
                                     variance_face_center.end());
      std::size_t dd = std::distance(variance_face_center.begin(), it_min);
      boundaries_dim[dd].push_back(average_face_center[dd]);
    }

  // Check which dimensions have periodic BCs.
  // If a dimension has periodic BCs they must be exactly two.
  std::vector<bool> dim_has_periodic_bcs(dim, 0);
  for (std::size_t dd = 0; dd < dim; ++dd)
    if (boundaries_dim[dd].size() != 0) {
      if (boundaries_dim[dd].size() != 2)
        throw std::runtime_error{"Extracted " +
                                 std::to_string(boundaries_dim[dd].size()) +
                                 " periodic boundaries "
                                 "along cartesian dimension " +
                                 std::to_string(dd)};
      dim_has_periodic_bcs[dd] = 1;
    }

  // Periodic conditions are checked in order across Cartesian dimensions.
  // If there are periodic BCs in all dimensions starting from 0 up to some
  // dimension, add only those boundaries If there are gap dimensions between
  // periodic dimensions, add {-inf, inf} boundaries.
  bool found_nonperiodic = 0;
  bool nonconsecutive_dims = 0;
  for (bool periodic : dim_has_periodic_bcs) {
    if (found_nonperiodic && periodic) {
      nonconsecutive_dims = 1;
      break;
    }
    if (!periodic)
      found_nonperiodic = 1;
  }
  std::vector<std::pair<double, double>> boundaries;
  boundaries.reserve(boundaries_dim.size());
  for (std::size_t dd = 0; dd < boundaries_dim.size(); ++dd) {
    if (dim_has_periodic_bcs[dd]) {
      boundaries.push_back(
          {std::min(boundaries_dim[dd][0], boundaries_dim[dd][1]),
           std::max(boundaries_dim[dd][0], boundaries_dim[dd][1])});
    } else if (nonconsecutive_dims) {
      boundaries.push_back({-std::numeric_limits<double>::infinity(),
                            std::numeric_limits<double>::infinity()});
    }
  }

  return boundaries;
}
} // namespace ptof

#endif /* PTOF_BOUNDARY_H */
