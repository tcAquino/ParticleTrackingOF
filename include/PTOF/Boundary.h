/**
 * @file   Boundary.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Thu Feb 17 00:00:00 2022
 *
 * @brief Objects and utilities to enforce boundary conditions.
 */

#ifndef PTOF_BOUNDARY_H
#define PTOF_BOUNDARY_H

#include "General/IO.h"
#include "General/Operation.h"
#include "PTOF/Directories.h"
#include "PTOF/DynamicsList.h"
#include "PTOF/Useful.h"
#include <algorithm>
#include <cstddef>
#include <fieldTypes.H>
#include <iterator>
#include <limits>
#include <point.H>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ptof {
/**
 * @brief Apply reflecting boundary condition.
 *
 * @param state Particle state to update.
 *
 * @param contact_point Contact point at boundary.
 *
 * @param unit_normal Unit normal at the contact point, pointing inwards.
 */
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
 * @brief Get the boundary conditions associated with a set of patches from
 *        file.
 *
 * @param filename Name of file holding a patch name followed by a boundary
 *                 condition type in each line.
 *
 * @param delims Possible delimiters between patch names and boundary condition
 *               types.
 */
inline auto get_boundary_conditions(std::string const &filename,
                                    std::string const &delims = "\t,| ") {
  std::unordered_map<std::string, std::string> boundary_conditions;

  auto input = io::open_read(filename);
  std::string in_file = std::string{"In file "} + filename + " : ";
  for (std::vector<std::string> split_line; io::split_line(input, split_line);
       split_line.clear()) {
    std::size_t param_index = 0;
    auto patch_name = io::read<std::string>(
        split_line, param_index, in_file + "Could not parse patch name");
    auto bc_name = io::read<std::string>(
        split_line, param_index,
        in_file + "Could not parse boundary condition type");
    boundary_conditions[patch_name] = bc_name;
  }

  return boundary_conditions;
}

/**
 * @brief Get boundary conditions for dynamics type.
 *
 * @param directories Current case directory information.
 *
 * @param delims Possible delimiters between patch names and boundary condition
 *               types.
 */
template <DynamicsList::Type dynamics>
auto get_boundary_conditions(Directories const &directories,
                             std::string const &delims = "\t,| ") {
  return get_boundary_conditions(directories.dir_boundaryconditions + "/" +
                                     DynamicsList::name(dynamics) + ".bcs",
                                 delims);
}

/**
 * @brief Throw if any patch name does not exist in mesh or associated boundary
 *        condition type is not implemented.
 *
 * @param boundary_conditions Container of pairs of patch names and boundary
 *                            condition types.
 *
 * @param mesh Mesh object.
 */
template <typename BoundaryCondition, typename Mesh,
          typename BoundaryConditionList>
void verify_boundary_conditions(BoundaryCondition const &boundary_conditions,
                                Mesh const &mesh, BoundaryConditionList) {
  for (auto const &bc : boundary_conditions) {
    mesh.boundaryMesh().findPatchID(bc.first, false);
  }

  for (auto const &bc : boundary_conditions) {
    if (!BoundaryConditionList::contains(bc.second)) {
      throw std::runtime_error{std::string{"Boundary condition type "} +
                               bc.second + " : " + "Not supported"};
    }
  }
}

/**
 * @brief Add default type 'empty' to boundary conditions for unspecified
 * patches.
 *
 * @param boundary_conditions Associative container of pairs of patch names and
 *                            boundary condition types, to add to.
 *
 * @param patch_names Container with all patch names.
 */
template <typename BoundaryCondition, typename Container>
void add_unspecified_patches(BoundaryCondition &boundary_conditions,
                             Container const &patch_names) {
  for (auto const &name : patch_names) {
    if (!boundary_conditions.count(name)) {
      boundary_conditions[name] = "empty";
    }
  }
}

/**
 * @brief Add default type 'empty' to boundary conditions for unspecified
 *        patches in mesh.
 *
 * @param boundary_conditions Associative container of pairs of patch names and
 *                            boundary condition types, to add to.
 *
 * @param mesh Mesh object.
 */
template <typename BoundaryCondition, typename Mesh>
void add_unspecified_patches_in_mesh(BoundaryCondition &boundary_conditions,
                                     Mesh const &mesh) {
  add_unspecified_patches(boundary_conditions, mesh.boundaryMesh().names());
}

/** @brief Boundary object with no effect. */
struct Boundary_DoNothing {
  /**
   * @brief Check if position is out of bounds.
   *
   * @param position Position to check.
   *
   * @return \c false (always in bounds).
   */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return false;
  }

  /**
   * @brief Apply boundary condition (do nothing).
   *
   * @param state Current particle state.
   *
   * @param state_old Previous particle state.
   *
   * @param dummy Dummy argument (unused).
   *
   * @return \c false (no effect).
   */
  template <typename State, typename Dummy>
  bool operator()(State const &state, State const &state_old,
                  Dummy const &dummy) const {
    return false;
  }

  /**
   * @brief Apply boundary condition (do nothing).
   *
   * @param state Current particle state.
   *
   * @param dummy Dummy argument (unused).
   *
   * @return \c false (no effect).
   */
  template <typename State, typename Dummy>
  bool operator()(State const &state, Dummy const &dummy) const {
    return false;
  }

  /**
   * @brief Apply boundary condition (do nothing).
   *
   * @param state Current particle state.
   *
   * @return \c false (no effect).
   */
  template <typename State> bool operator()(State const &state) const {
    return false;
  }

  /** @brief Name of boundary condition type. */
  std::string name() const { return "empty"; }
};

/**
 * @brief Find periodic image of intersection point.
 *
 * @param state_outside Out-of-bounds state.
 *
 * @param state_image Periodic image of state, in bounds.
 *
 * @param intersection Boundary intersection info.
 *
 * @param boundary Periodic boundary enforcer.
 *
 * @param locator Object to locate positions in mesh.
 *
 * @return Periodic image of intersection.
 */
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
  if (Foam::magSqr(offset) < Foam::magSqr(displacement)) {
    state_outside.set_position(offset_pos);
  } else {
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
  for (std::size_t iter = 0; iter < max_nr_doublings; ++iter) {
    offset *= 2;
    auto new_intersection = locator.mesh_search().intersection(
        make_point(state_outside.position) - offset,
        make_point(state_outside.position));
    if (new_intersection.hit()) {
      return new_intersection;
    }
  }

  // If finding the intersection failed, take the closest boundary face center
  auto face_id = nearest_boundary_face(state_image.position, locator);

  // This should not happen, but it can for bad position inputs and it can
  // lead to bad memory access, so it is safer to check
  if (face_id == -1) {
    throw std::runtime_error{
        "Could not find image of periodic boundary intersection"};
  }

  return Intersection(true, face_center(face_id, locator.mesh()), face_id);
}

/** @brief Interface for periodic boundary enforcer. */
template <typename Boundary, typename Locator> struct Boundary_Periodic {
  /**
   * @brief Constructor.
   *
   * @param boundary Periodic boundary enforcer.
   *
   * @param locator Object to locate positions in mesh.
   */
  Boundary_Periodic(Boundary &&boundary, Locator &&locator)
      : boundary_periodic{std::forward<Boundary>(boundary)},
        locator{std::forward<Locator>(locator)} {}

  /**
   * @brief Check if position is out of bounds.
   *
   * @param position Position to check
   *
   * @return \c true if out of bounds, \c false otherwise.
   */
  template <typename Position>
  bool out_of_bounds(Position const &position) const {
    return boundary_periodic.out_of_bounds(position);
  }

  /**
   * @brief Enforce boundary condition and place intersection at its periodic
   *        image.
   *
   * @param state Particle state to apply BC to.
   *
   * @param intersection Boundary intersection info, to be replaced by periodic
   *                     image.
   * @return \c true.
   */
  template <typename State, typename Intersection>
  bool operator()(State &state, Intersection &intersection) const {
    auto state_outside = state;
    boundary_periodic(state);
    intersection = periodic_intersection(state_outside, state, intersection,
                                         boundary_periodic, locator);

    return true;
  }

  Boundary boundary_periodic; /**< Periodic boundary type. */
  Locator locator;            /**< Object to locate positions in mesh. */
};
template <typename Boundary, typename Locator>
Boundary_Periodic(Boundary &&, Locator &&)
    -> Boundary_Periodic<Boundary, Locator>;

/**
 * @brief Boundary object that reinjects particles according to prescribed
 *        initial condition.
 */
template <typename InitialCondition> struct Boundary_Reinject {
  /**
   * @brief Constructor.
   *
   * @param initial_condition Initial condition object.
   */
  Boundary_Reinject(InitialCondition &&initial_condition)
      : initial_condition{std::forward<InitialCondition>(initial_condition)} {}

  /**
   * @brief Enforce boundary condition.
   *
   * @param state Particle state to apply BC to.
   *
   * @param state_old Old particle state (unused).
   *
   * @param intersection Boundary intersection info (unused).
   *
   * @return \c true
   */
  template <typename State, typename Intersection>
  bool operator()(State &state, State const &state_old,
                  Intersection const &intersection) const {
    auto [position, cell] = initial_condition.make_position_and_cell();
    state.set_position(position, cell);
    return true;
  }

  /** @brief Name of boundary condition type. */
  std::string name() const { return "reinject"; }

  InitialCondition initial_condition; /**< Initial condition according to which
                                         to reinject. */
};
template <typename InitialCondition>
Boundary_Reinject(InitialCondition &&) -> Boundary_Reinject<InitialCondition>;

/**
 * @brief For boundary conditions indicated as type \c periodic, extract the
 *        corresponding boundary positions.
 *
 * @param boundary_conditions Associative container of patch names and
 *                            associated boundary condition types.
 *
 * @param mesh Mesh object.
 *
 * @param dim Spatial dimensionality.
 *
 * @return Container of pairs of periodic boundary positions along each periodic
 *         dimension.
 *
 * @note
 *
 * - Periodic boundaries are assumed to be perpendicular to the cartesian axes.
 *   There must be an even number.
 *
 *  - Periodic conditions are checked in order across Cartesian dimensions. If
 *    there are periodic BCs in all dimensions starting from 0 up to some
 *    dimension, output only those boundaries. If there are gap dimensions
 *    between periodic dimensions, add <tt>{-inf, inf}</tt> boundaries.
 */
template <typename BoundaryCondition, typename Mesh>
auto extract_cartesian_periodic_boundaries(
    BoundaryCondition const &boundary_conditions, Mesh const &mesh,
    std::size_t dim) {
  std::vector<std::vector<double>> boundaries_dim(dim);
  for (auto const &bc : boundary_conditions) {
    if (bc.second != "periodic") {
      continue;
    }
    // Find average and variance of patch face centers.
    auto patch_id = mesh.boundaryMesh().findPatchID(bc.first, false);
    auto const &patch = mesh.boundaryMesh()[patch_id];
    auto face_start = patch.start();
    Foam::point average_face_center = Foam::zero{};
    std::vector<double> variance_face_center(dim, 0.);
    Foam::scalar total_area = 0.;
    for (auto face = face_start; face < face_start + patch.size(); ++face) {
      auto area = face_area(face, mesh);
      auto center = face_center(face, mesh);
      auto weighted_center = area * center;
      total_area += area;
      average_face_center += weighted_center;
      for (std::size_t dd = 0; dd < variance_face_center.size(); ++dd) {
        variance_face_center[dd] += weighted_center[dd] * center[dd];
      }
    }
    average_face_center /= total_area;
    op::div_scalar_inplace(variance_face_center, total_area);

    for (std::size_t dd = 0; dd < variance_face_center.size(); ++dd) {
      variance_face_center[dd] -=
          average_face_center[dd] * average_face_center[dd];
    }

    // Choose dimension with least variance.
    auto it_min = std::min_element(variance_face_center.begin(),
                                   variance_face_center.end());
    std::size_t dd = std::distance(variance_face_center.begin(), it_min);
    boundaries_dim[dd].push_back(average_face_center[dd]);
  }

  // Check which dimensions have periodic BCs.
  // If a dimension has periodic BCs they must be exactly two.
  std::vector<bool> dim_has_periodic_bcs(dim, 0);
  for (std::size_t dd = 0; dd < dim; ++dd) {
    if (boundaries_dim[dd].size() != 0) {
      if (boundaries_dim[dd].size() != 2) {
        throw std::runtime_error{"Extracted " +
                                 std::to_string(boundaries_dim[dd].size()) +
                                 " periodic boundaries "
                                 "along cartesian dimension " +
                                 std::to_string(dd)};
      }
      dim_has_periodic_bcs[dd] = 1;
    }
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
    if (!periodic) {
      found_nonperiodic = 1;
    }
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
