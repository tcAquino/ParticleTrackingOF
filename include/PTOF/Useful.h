/**
   \file PTOF/Useful.h
   \author Tomas Aquino
   \date 2022/02/17
   \brief Miscelaneous objects and utilities.
*/

#ifndef PTOF_USEFUL_H
#define PTOF_USEFUL_H

#include "CTRW/StateGetter.h"
#include "General/IO.h"
#include "General/Meta.h"
#include "General/Operation.h"
#include "General/Ranges.h"
#include "General/Useful.h"
#include "PTOF/Meta.h"
#include <Vector2D.H>
#include <fieldTypes.H>
#include <iostream>
#include <point.H>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace ptof {
/**
   \class PositionAndCell PTOF/Useful.h "PTOF/Useful.h"
   \brief Object to hold a position and a cell for its location in the mesh.
*/
template <typename Position> struct PositionAndCell {

  PositionAndCell(Position position, Foam::label cell = -1)
      : position{position}, cell{cell} {}

  Position position;
  Foam::label cell;
};

/**
   \brief Print static info for \tparam Class and Class::Parameters, if
   defined.
   \param output Output stream to print info.
   \param notify_if_no_info Print notification if no info is available.
   \param help_option Information about requested info, to print if info is not
   available if notification is requested.
   \return \c true if there is available info, \c false otherwise.
*/
template <typename Class>
bool print_static_info(std::ostream &output, bool notify_if_no_info = false,
                       std::string const &help_option = {}) {
  bool has_info = false;
  if constexpr (meta::has_static_info_v<Class>) {
    output << "\n";
    Class::info(output);
    has_info = true;
  }
  if constexpr (meta::has_parameters_type_v<Class>)
    if constexpr (meta::has_static_info_v<typename Class::Parameters>) {
      output << "\n";
      Class::Parameters::info(output);
      has_info = true;
    }

  if (notify_if_no_info && !has_info) {
    output << "\n"
              "No static info available";
    if (!help_option.empty()) {
      output << " for help option " << help_option;
    }
    if (std::is_same_v<Class, meta::Empty>) {
      output << " : object type is Empty";
      output << "\n";
    }
  }

  return has_info;
}

/**
   \brief Handle help options for main executable.
   \return \c true if help flag was passed as first argument, \c false otherwise
   \note meta::Empty can be passed for non-existent templated types, except
   \c ExecutableInfo, which must define a static <tt>help(std::ostream&)</tt>
   method.
*/
template <typename ExecutableInfo, typename Geometry, typename DirectoriesOF,
          typename Transport, typename Phase, typename Reaction,
          typename Solvers, typename InitialCondition, typename Output>
bool options_help(std::ostream &output, int argc, const char *const *argv) {
  if (!io::check_options_help(argc, argv)) {
    return 0;
  }

  if (argc == 2) {
    output << "\n";
    ExecutableInfo::help(output);
  }

  for (int ii = 2; ii < argc; ++ii) {
    std::string option = std::string{argv[ii]};
    if (option == "-a" || option == "--all") {
      bool info = false;
      info += print_static_info<ExecutableInfo>(output);
      info += print_static_info<Geometry>(output);
      info += print_static_info<DirectoriesOF>(output);
      info += print_static_info<Transport>(output);
      info += print_static_info<Phase>(output);
      info += print_static_info<Reaction>(output);
      info += print_static_info<Solvers>(output);
      info += print_static_info<InitialCondition>(output);
      info += print_static_info<Output>(output);
      if (!info) {
        output << "No help info available\n";
      }
    } else if (option == "-e" || option == "--executable") {
      print_static_info<ExecutableInfo>(output, true, "-m / --main");
    } else if (option == "-g" || option == "--geometry") {
      print_static_info<Geometry>(output, true, "-g / --geometry");
    } else if (option == "-d" || option == "--directories-of") {
      print_static_info<DirectoriesOF>(output, true, "-d / --directories-of");
    } else if (option == "-t" || option == "--transport") {
      print_static_info<Transport>(output, true, "-t / --transport");
    } else if (option == "-p" || option == "--phase") {
      print_static_info<Phase>(output, true, "-p / --phase");
    } else if (option == "-r" || option == "--reaction") {
      print_static_info<Reaction>(output, true, "-r / --reaction");
    } else if (option == "-s" || option == "--solvers") {
      print_static_info<Solvers>(output, true, "-s / --solvers");
    } else if (option == "-i" || option == "--initial-condition") {
      print_static_info<InitialCondition>(output, true,
                                          "-i / --initial-condition");
    } else if (option == "-o" || option == "--output") {
      print_static_info<Output>(output, true, "-o / --output");
    } else {
      output << "\n"
                "Help option "
             << option << " : " << "Not supported\n";
    }
  }

  return 1;
}

/**
   \brief Return whether cell is outside mesh.
   \return true if outside, false otherwise.
   \note \p cell is the actual index of the cell position is in, determined
   elsewhere, not a hint.
*/
bool outside(Foam::label cell) {
  if (cell < 0) {
    return 1;
  }
  return 0;
}

/**
   \brief Check whether cell index is outside mesh.
   \tparam warn_if_outside Output warning if true, do not if false.
   \param position Position to check.
   \param cell Mesh index of cell position is in.
   \param extra_warning_info Additional info if warning is issued.
   \return true if outside, false otherwise.
   \note \p cell is the actual index of the cell position is in, determined
   elsewhere, not a hint.
*/
template <bool warn_if_outside>
bool outside(Foam::label cell, Foam::point const &position,
             std::string const &extra_warning_info = {}) {
  if (cell < 0) {
    if constexpr (warn_if_outside) {
      std::cerr << "Warning: Requested cell at position " << "(" << position[0]
                << ", " << position[1] << ", " << position[2] << ")"
                << " outside mesh.";
      if (!extra_warning_info.empty()) {
        std::cerr << " " << extra_warning_info << "\n";
      }
    }
    return 1;
  }
  return 0;
}

/**
   \brief Check if \c time is between the two particle state times.
   \param particle Particle to check.
   \param time Time to check.
   \return \c true if \c time is between two states, \c false otherwise.
   \note Particle states must define:
   - \c time
*/
template <typename Particle>
bool brackets_time(Particle const &particle, double time) {
  return particle.state_new().time >= time && particle.state_old().time <= time;
}

/**
   \brief Check if a particle is absorbed at time \c time
   \details Checks Particle::State::Info for \c absorbed, returns \c false if it
   does not exist.
   \param particle Particle to interpolate between old and new state.
   \param time Time to check.
   \return \c true if \c absorbed is true, \c false otherwise.
   \note
   - The absorbed state is considered special and meant to be irreversible;
   otherwise, the adsorbed state should be used.
   - Particle states must define:
   -# <tt>Info info</tt>
   -# \c time
*/
template <typename Particle>
bool absorbed(Particle const &particle, double time) {
  using Info = typename Particle::State::Info;
  if constexpr (meta::has_absorbed_v<Info>) {
    if (particle.state_new().info.absorbed &&
        particle.state_new().time <= time) {
      return true;
    }
  }
  return false;
}

/**
   \brief Check if a state is adsorbed.
   \details Checks Particle::State::Info for \c absorbed, returns \c false if
   it does not exist.
   \return \c true if \c absorbed is true, \c false otherwise.
   \note State must define:
   - \c Info
*/
template <typename State> bool adsorbed(State const &state) {
  using Info = typename State::Info;
  if constexpr (meta::has_adsorbed_v<Info>) {
    if (state.info.adsorbed) {
      return true;
    }
  }
  return false;
}

/** \brief Make 3D point from 2D point. */
auto make_point(Foam::Vector2D<Foam::scalar> const &point) {
  return Foam::point{point[0], point[1], 0.};
}

/** \brief Make 3D point from 1D point. */
auto make_point(Foam::scalar point) { return Foam::point{point, 0., 0.}; }

/** \brief Make 3D point from 3D point (simple copy). */
auto make_point(Foam::point const &point) { return point; }

/** \brief Make point from vector. */
Foam::point make_point(std::vector<double> const &point) {
  if (point.size() >= 3) {
    return {point[0], point[1], point[2]};
  }
  if (point.size() == 2) {
    return Foam::point{point[0], point[1], 0.};
  }
  if (point.size() == 1) {
    return Foam::point{point[0], 0., 0.};
  }
  return {0., 0., 0.};
}

/** \brief Oriented area at boundary face, using right-hand-rule (outward for
    boundary faces). */
template <typename Mesh> auto area_outward(Foam::label face, Mesh const &mesh) {
  return mesh.faces()[face].areaNormal(mesh.points());
}

/** \brief Oriented area at boundary face, using left-hand-rule (inward for
  boundary faces). */
template <typename Mesh> auto area_inward(Foam::label face, Mesh const &mesh) {
  return -area_outward(face, mesh);
}

/** \brief Unit normal at boundary face, using right-hand-rule (outward for
    boundary faces). */
template <typename Mesh>
auto unit_normal_outward(Foam::label face, Mesh const &mesh) {
  return mesh.faces()[face].unitNormal(mesh.points());
}

/** \brief Unit normal at boundary face, using left-hand-rule (inward for
  boundary faces). */
template <typename Mesh>
auto unit_normal_inward(Foam::label face, Mesh const &mesh) {
  return -unit_normal_outward(face, mesh);
}

/** \brief Area (magnitude) of a face. */
template <typename Mesh> auto face_area(Foam::label face, Mesh const &mesh) {
  return mesh.faces()[face].mag(mesh.points());
}

/** \brief Center of a face. */
template <typename Mesh> auto face_center(Foam::label face, Mesh const &mesh) {
  return mesh.faces()[face].centre(mesh.points());
}

/** \brief Volume of a cell. */
template <typename Mesh> auto cell_volume(Foam::label cell, Mesh const &mesh) {
  return mesh.cellVolumes()[cell];
}

/** \brief Center of a cell. */
template <typename Mesh> auto cell_center(Foam::label cell, Mesh const &mesh) {
  return mesh.cellCentres()[cell];
}

/** \brief Cell volumes. */
template <typename Container, typename Mesh>
auto cell_volumes(Container const &cell_ids, Mesh const &mesh) {
  std::vector<double> volumes;
  volumes.reserve(cell_ids.size());
  for (auto cell : cell_ids)
    volumes.push_back(cell_volume(cell, mesh));
  return volumes;
}

/** \brief Find patch id in mesh for a face known to be in a boundary patch. */
template <typename Mesh> auto patch_id(Foam::label face, Mesh const &mesh) {
  return mesh.boundaryMesh().whichPatch(face);
}

/** \brief Find patch id in mesh for a patch name. */
template <typename Mesh>
auto patch_id(Foam::word patch_name, Mesh const &mesh) {
  return mesh.boundary().findPatchID(patch_name);
}

/** \brief Face areas (magnitudes). */
template <typename Container, typename Mesh>
auto face_areas(Container const &face_ids, Mesh const &mesh) {
  std::vector<double> areas;
  areas.reserve(face_ids.size());
  for (auto face : face_ids)
    areas.push_back(face_area(face, mesh));
  return areas;
}

/** \brief Total area of faces. */
template <typename Container, typename Mesh>
auto face_area(std::vector<Foam::label> const &face_ids, Mesh const &mesh) {
  double area = 0.;
  for (auto const &face : face_ids)
    area += face_area(face, mesh);
  return area;
}

/** \brief Total area of faces in patch. */
template <typename Mesh>
auto patch_area(std::string const &patch_name, Mesh const &mesh) {
  double area = 0.;
  Foam::label patch_id = mesh.boundary().findPatchID(patch_name);
  auto const &patch = mesh.boundaryMesh()[patch_id];
  Foam::label start = patch.start();
  for (std::size_t idx = 0; idx < static_cast<std::size_t>(patch.size()); ++idx)
    area += face_area(start + idx, mesh);
  return area;
}

/** \brief Total area of faces in patches. */
template <typename Mesh>
auto patch_area(std::vector<std::string> const &patch_names, Mesh const &mesh) {
  double area = 0.;
  for (auto const &patch : patch_names)
    area += patch_area(patch, mesh);
  return area;
}

/** \brief Total areas of faces in each patch. */
template <typename Mesh>
auto patch_areas(std::vector<std::string> const &patch_names,
                 Mesh const &mesh) {
  std::vector<double> areas;
  for (auto const &patch : patch_names)
    areas.push_back(patch_area(patch, mesh));
  return areas;
}

/**
   \brief Cell volume weighted by magnitude of vector field value at cell
   center.
*/
template <typename VectorField, typename Mesh>
auto cell_flux(Foam::label cell_id, VectorField const &vector_field,
               Mesh const &mesh) {
  return cell_volume(cell_id, mesh) *
         Foam::mag(vector_field(cell_center(cell_id, mesh)));
}

/**
   \brief Face area weighted by outward-normal component of vector field at face
   center. */
template <typename VectorField, typename Locator>
auto face_flux_outward(Foam::label face_id, VectorField const &vector_field,
                       Locator const &locator) {
  auto const &mesh = locator.mesh();
  auto patch_id = ptof::patch_id(face_id, mesh);
  auto const &patch = mesh.boundaryMesh()[patch_id];
  Foam::label start = patch.start();

  return area_outward(face_id, mesh) &
         vector_field.boundaryField()[patch_id][face_id - start];
}

/**
   \brief Face area weighted by inward-normal component of vector field at face
   center. */
template <typename VectorField, typename Locator>
auto face_flux_inward(Foam::label face_id, VectorField const &vector_field,
                      Locator const &locator) {
  auto const &mesh = locator.mesh();
  auto patch_id = ptof::patch_id(face_id, mesh);
  auto const &patch = mesh.boundaryMesh()[patch_id];
  Foam::label start = patch.start();

  return area_inward(face_id, mesh) &
         vector_field.boundaryField()[patch_id][face_id - start];
}

/**
   \brief Cell volumes weighted by magnitude of vector field value at cell
   centers.
*/
template <typename Container, typename VectorField, typename Mesh>
auto cell_fluxes(Container const &cell_ids, VectorField const &vector_field,
                 Mesh const &mesh) {
  std::vector<double> fluxes;
  fluxes.reserve(cell_ids.size());
  for (auto cell : cell_ids)
    fluxes.push_back(cell_flux(cell, vector_field, mesh));
  return fluxes;
}

/**
   \brief Face areas weighted by outward-normal component of vector field at
   face center. */
template <typename Container, typename VectorField, typename Locator>
auto face_fluxes_outward(Container const &face_ids,
                         VectorField const &vector_field,
                         Locator const &locator) {
  std::vector<double> fluxes;
  fluxes.reserve(face_ids.size());
  for (auto face : face_ids)
    fluxes.push_back(face_flux_outward(face, vector_field, locator));
  return fluxes;
}

/**
   \brief Face areas weighted by inward-normal component of vector field at
   face center. */
template <typename Container, typename VectorField, typename Locator>
auto face_fluxes_inward(Container const &face_ids,
                        VectorField const &vector_field,
                        Locator const &locator) {
  std::vector<double> fluxes;
  fluxes.reserve(face_ids.size());
  for (auto face : face_ids)
    fluxes.push_back(face_flux_inward(face, vector_field, locator));
  return fluxes;
}

/**
   \brief Sometimes the face center associated with a mesh face is not
   considered within the cell and may be outside the mesh. Verify this.
   \param face Mesh face index.
   \param locator Object to locate positions in mesh.
   \return \c true if face center is in mesh, \c false otherwise.
*/
template <typename Locator>
bool face_center_is_in_mesh(Foam::label face, Locator const &locator) {
  auto const &mesh = locator.mesh();
  return !outside(locator(face_center(face, mesh), mesh.faceOwner()[face]));
}

/**
   \brief Sometimes the face center associated with a mesh face is not
   considered within the owner cell. Verify this.
   \param face Mesh face index.
   \param locator Object to locate positions in mesh.
   \return \c true if face center is in owner cell, \c false otherwise.
*/
template <typename Locator>
bool face_center_is_in_cell(Foam::label face, Locator const &locator) {
  auto const &mesh = locator.mesh();
  auto owner_cell = mesh.faceOwner()[face];
  auto cell = locator(face_center(face, mesh), owner_cell);
  return !outside(cell) && cell == owner_cell;
}

/**
   \brief Compute face center's position if face center is in owner cell, owner
   cell center otherwise.
   \param face Mesh face index.
   \param locator Object to locate positions in mesh.
   \return Position.
*/
template <typename Locator>
auto adjusted_face_center(Foam::label face, Locator const &locator) {
  auto const &mesh = locator.mesh();
  if (face_center_is_in_cell(face, locator)) {
    return face_center(face, mesh);
  } else {
    auto center = face_center(face, mesh);
    std::cerr << "Warning: Face center of face " << face << " at " << "("
              << center[0] << ", " << center[1] << ", " << center[2] << ")"
              << " is not within owner cell. "
              << "Replacing face center by owner cell center\n";
    return cell_center(mesh.faceOwner()[face], mesh);
  }
}

/** \brief Small offset forward given current face and direction. */
template <typename Locator>
Foam::vector offset_face(Foam::label face, Foam::vector const &direction,
                         Locator const &locator) {
  auto const &mesh = locator.mesh();
  Foam::label owner_cell = mesh.faceOwner()[face];
  Foam::point const &cell_center = mesh.cellCentres()[owner_cell];
  Foam::scalar typ_dim = Foam::mag(cell_center - face_center(face, mesh));

  return locator.mesh_search().tol_ * typ_dim * direction /
         Foam::mag(direction);
}

/**
   \brief Small offset along face normal (outward).
   \param face Mesh face index to use face center as beginning point.
   \param locator Object to locate positions in mesh.
   \return Offset vector.
*/
template <typename Locator>
Foam::vector offset_face(Foam::label face, Locator const &locator) {
  auto const &mesh = locator.mesh();
  return offset_face(face, unit_normal_outward(face, mesh), locator);
}

/**
   \brief Small offset forward from a point on a face.
   \param begin Starting point (on a face).
   \param face Mesh face index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
*/
template <typename Locator>
Foam::point offset_forward_face(Foam::point const &begin, Foam::label face,
                                Foam::vector const &direction,
                                Locator const &locator) {
  return begin + offset_face(face, direction, locator);
}

/**
   \brief Small offset forward from a point on a face.
   \details If begin point is in mesh, guarantee offset point in same cell.
   \param begin Starting point (on a face).
   \param face Mesh face index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
*/

template <typename Locator>
Foam::point offset_forward_face_keep_inside(Foam::point const &begin,
                                            Foam::label face,
                                            Foam::vector const &direction,
                                            Locator const &locator) {
  auto offset = offset_face(face, direction, locator);
  auto point = begin + offset;
  Foam::label owner_cell = locator.mesh().faceOwner()[face];
  auto begin_cell = locator(begin, owner_cell);
  if (outside(begin_cell)) {
    return point;
  }
  while (locator(point, begin_cell) != begin_cell) {
    offset /= 2.;
    point = begin + offset;
  }

  return point;
}

/**
   \brief Small offset backward from a point on a face.
   \param begin Starting point (on a face).
   \param face Mesh face index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
*/
template <typename Locator>
Foam::point offset_backward_face(Foam::point const &begin, Foam::label face,
                                 Foam::vector const &direction,
                                 Locator const &locator) {
  return begin - offset_face(begin, face, direction, locator);
}

/**
   \brief Small offset backward from a point on a face.
   \details If begin point is in mesh, guarantee offset point in same cell.
   \param begin Starting point (on a face).
   \param face Mesh face index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
*/
template <typename Locator>
Foam::point offset_backward_face_keep_inside(Foam::point const &begin,
                                             Foam::label face,
                                             Foam::vector const &direction,
                                             Locator const &locator) {
  auto offset = offset_face(face, direction, locator);
  auto point = begin - offset;
  Foam::label owner_cell = locator.mesh().faceOwner()[face];
  auto begin_cell = locator(begin, owner_cell);
  if (outside(begin_cell)) {
    return point;
  }
  while (locator(point, begin_cell) != begin_cell) {
    offset /= 2.;
    point = begin - offset;
  }
  return point;
}

/**
   \brief Small offset outward along face normal starting from face center.
   \details If begin point is in mesh, guarantee offset point in mesh.
   \param face Mesh face index.
   \param locator Object to locate positions in mesh.
   \return Offset point.
*/
template <typename Locator>
Foam::point offset_outward_face(Foam::label face, Locator const &locator) {
  return face_center(face, locator.mesh()) + offset_face(face, locator);
}

/**
   \brief Small offset inward along face normal starting from face center.
   \details If begin point is in mesh, guarantee offset point in mesh.
   \param face Mesh face index.
   \param locator Object to locate positions in mesh.
   \return Offset point.
*/
template <typename Locator>
Foam::point offset_inward_face(Foam::label face, Locator const &locator) {
  return face_center(face, locator.mesh()) - offset_face(face, locator);
}

/**
   \brief Small offset outward along face normal starting from face center.
   \details If begin point is in mesh, guarantee offset point in mesh.
   \param begin Starting point in cell.
   \param cell Mesh cell index.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
*/
template <typename Locator>
Foam::vector offset_cell(Foam::point const &begin, Foam::label cell,
                         Foam::vector const &direction,
                         Locator const &locator) {
  Foam::point const &center = locator.mesh().cellCentres()[cell];
  Foam::scalar typ_dim = Foam::mag(center - begin);

  return locator.mesh_search().tol_ * typ_dim * direction /
         Foam::mag(direction);
}

/**
   \brief Small offset forward.
   \param begin Starting point in cell.
   \param cell Hint for cell index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
*/
template <typename Locator>
Foam::point offset_forward_cell(Foam::point const &begin, Foam::label cell,
                                Foam::vector const &direction,
                                Locator const &locator) {
  return begin + offset_cell(begin, cell, direction, locator);
}

/**
   \brief Small offset forward.
   \details If begin point is in mesh, guarantee offset point in same cell.
   \param begin Starting point in cell.
   \param cell Hint for cell index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
*/
template <typename Locator>
Foam::point offset_forward_cell_keep_inside(Foam::point const &begin,
                                            Foam::label cell,
                                            Foam::vector const &direction,
                                            Locator const &locator) {
  auto offset = offset_cell(begin, cell, direction, locator);
  auto point = begin + offset;
  if (outside(cell)) {
    return point;
  }
  while (locator(point, cell) != cell) {
    offset /= 2.;
    point = begin + offset;
  }
  return point;
}

/**
   \brief Small offset backward.
   \param begin Starting point in cell.
   \param cell Hint for cell index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
*/
template <typename Locator>
Foam::point offset_backward_cell(Foam::point const &begin, Foam::label cell,
                                 Foam::vector const &direction,
                                 Locator const &locator) {
  return begin - offset_cell(begin, cell, direction, locator);
}

/**
   \brief Small offset backward.
   \details If begin point is in mesh, guarantee offset point in same cell.
   \param begin Starting point in cell.
   \param cell Hint for cell index where \p begin is located.
   \param direction Direction of offset.
   \param locator Object to locate positions in mesh.
   \return Offset point.
*/
template <typename Locator>
Foam::point offset_backward_cell_keep_inside(Foam::point const &begin,
                                             Foam::label cell,
                                             Foam::vector const &direction,
                                             Locator const &locator) {
  auto offset = offset_cell(begin, cell, direction, locator);
  auto point = begin - offset;
  if (outside(cell)) {
    return point;
  }
  while (locator(point, cell) != cell) {
    offset /= 2.;
    point = begin - offset;
  }
  return point;
}

/**
   \param mesh Mesh object.
   \return All mesh cell indices.
*/
template <typename Mesh> auto all_cell_ids(Mesh const &mesh) {
  return range::range<std::vector>(0, mesh.nCells());
}

/**
   \brief Append mesh face indices of faces in patch to container.
   \param patch_name Name of patch.
   \param mesh Mesh object.
   \param face_ids Face index container to append to.
*/
template <typename Mesh>
void patch_face_ids(std::string const &patch_name, Mesh const &mesh,
                    std::vector<Foam::label> &face_ids) {
  Foam::label patch_id = ptof::patch_id(patch_name, mesh);
  auto const &patch = mesh.boundaryMesh()[patch_id];
  face_ids.reserve(face_ids.size() + patch.size());
  Foam::label start = patch.start();
  for (std::size_t idx = 0; idx < static_cast<std::size_t>(patch.size()); ++idx)
    face_ids.push_back(start + idx);
}

/**
 \param patch_names Names of patches.
 \param mesh Mesh object.
 \return Mesh face indices of faces in patchs.
*/
template <typename Mesh>
auto patch_face_ids(std::vector<std::string> const &patch_names,
                    Mesh const &mesh) {
  std::vector<Foam::label> face_ids;
  for (auto const &name : patch_names)
    patch_face_ids(name, mesh, face_ids);
  return face_ids;
}

/**
   \param boundaries Container of pairs of lower and upper boundary locations
   along each dimension.
   \param locator Object to locate positions in mesh.
   \return Mesh cell indices within boundaries.
*/
template <typename Locator>
auto cell_ids_region_cartesian(
    std::vector<std::pair<double, double>> const &boundaries,
    Locator const &locator) {
  // Identify degenerate dimensions
  // (if some minimum and maximum boundary values are equal,
  // the Cartesian region has lower dimension than the space)
  std::vector<std::size_t> degenerate_dimensions;
  std::vector<std::size_t> non_degenerate_dimensions;
  for (std::size_t dd = 0; dd < boundaries.size(); ++dd)
    if (boundaries[dd].first == boundaries[dd].second) {
      degenerate_dimensions.push_back(dd);
    } else
      non_degenerate_dimensions.push_back(dd);

  auto const &mesh = locator.mesh();

  std::set<Foam::label> cell_ids;
  for (Foam::label cc = 0; cc < mesh.nCells(); ++cc) {
    auto center = cell_center(cc, mesh);

    bool cell_is_within_non_degenerate_boundaries = 1;
    for (auto dd : non_degenerate_dimensions) {
      if (center[dd] < boundaries[dd].first ||
          center[dd] > boundaries[dd].second) {
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
    auto cell_id = locator(center);
    if (!outside(cell_id))
      cell_ids.insert(cell_id);
  }

  std::vector<Foam::label> cells(cell_ids.begin(), cell_ids.end());
  return cells;
}

/**
   \brief Remove cell indices if mask is below threshold.
   \param cell_ids Container with mesh cell indices.
   \param mask Scalar field.
   \param threshold Threshold.
*/
template <typename Container, typename Mask>
void apply_mask_cells_inplace(Container &cell_ids, Mask const &mask,
                              double threshold = 0.) {
  useful::swap_erase_if(cell_ids, [&mask, threshold](Foam::label face_id) {
    return mask[face_id] < threshold;
  });
}

/**
   \brief Remove face indices if mask is below threshold.
   \param cell_ids Container with mesh cell indices.
   \param mask Scalar field.
   \param threshold Threshold.
   \return \p Container with indices of elements where mask is below threshold
   removed.
*/
template <typename Container, typename Mask>
auto apply_mask_cells(Container const &cell_ids, Mask const &mask,
                      double threshold = 0.) {
  auto masked_ids = cell_ids;
  apply_mask_cells_inplace(masked_ids, mask, threshold);

  return masked_ids;
}

/**
   \brief Remove face indices if mask is below threshold.
   \param face_ids Container with mesh boundary face indices.
   \param mesh Mesh object.
   \param mask Scalar field.
   \param threshold Threshold.
*/
template <typename Container, typename Mesh, typename Mask>
void apply_mask_patch_faces_inplace(Container &face_ids, Mesh const &mesh,
                                    Mask const &mask, double threshold = 0.) {
  useful::swap_erase_if(
      face_ids, [&mask, &mesh, threshold](Foam::label face_id) {
        Foam::label patch_id = mesh.boundaryMesh().whichPatch(face_id);
        return mask.boundaryField()[patch_id][face_id] < threshold;
      });
}

/**
   \brief Remove face indices if mask is below threshold.
   \param face_ids Container with mesh face indices.
   \param mesh Mesh object.
   \param mask Scalar field.
   \param threshold Threshold.
   \return \p Container with indices of elements where mask is below threshold
   removed.
*/
template <typename Container, typename Mesh, typename Mask>
auto apply_mask_patch_faces(Container const &face_ids, Mesh const &mesh,
                            Mask const &mask, double threshold = 0.) {
  auto masked_ids = face_ids;
  apply_mask_patch_faces_inplace(masked_ids, mask, threshold);

  return masked_ids;
}

/**
   \brief Distribution for random number generation of cell indices weighted by
   cell volumes.
   \param cell_ids Container with mesh cell indices.
   \param mesh Mesh object.
*/
template <typename Container, typename Mesh>
auto uniform_cell_distribution(Container const &cell_ids, Mesh const &mesh) {
  auto weights = cell_volumes(cell_ids, mesh);
  return std::discrete_distribution<std::size_t>{weights.begin(),
                                                 weights.end()};
}
/**
   \brief Distribution for random number generation of face indices weighted by
   face areas.
   \param face_ids Container with mesh face indices.
   \param mesh Mesh object.
*/
template <typename Container, typename Mesh>
auto uniform_face_distribution(Container const &face_ids, Mesh const &mesh) {
  auto weights = face_areas(face_ids, mesh);
  return std::discrete_distribution<std::size_t>{weights.begin(),
                                                 weights.end()};
}

/**
   \brief Distribution for random number generation of cell indices weighted by
   cell volumes multiplied by vector field magnitude at cell center.
   \param cell_ids Container with mesh cell indices.
   \param vector_field Vector field.
   \param mesh Mesh object.
*/
template <typename Container, typename VectorField, typename Mesh>
auto fluxweighted_cell_distribution(Container const &cell_ids,
                                    VectorField const &vector_field,
                                    Mesh const &mesh) {
  auto weights = cell_fluxes(cell_ids, vector_field, mesh);
  if (op::sum(weights) == 0.)
    throw std::runtime_error{"Cannot define flux-weighted distribution because "
                             "all cells have zero flux"};
  return std::discrete_distribution<std::size_t>{weights.begin(),
                                                 weights.end()};
}

/**
   \brief Distribution for random number generation of face indices weighted by
   face areas multiplied by vector field component along inward normal.
   \note Faces with vector field pointing outward have weight zero.
*/
template <typename Container, typename VectorField, typename Locator>
auto fluxweighted_face_distribution(Container const &face_ids,
                                    VectorField const &vector_field,
                                    Locator const &locator) {
  auto weights = face_fluxes_inward(face_ids, vector_field, locator);
  op::apply([](auto val) { return std::max(0., val); }, weights);
  if (op::sum(weights) == 0.)
    throw std::runtime_error{"Cannot define flux-weighted distribution because "
                             "all faces have zero flux"};
  return std::discrete_distribution<std::size_t>{weights.begin(),
                                                 weights.end()};
}

/**
   \param subject CTRW object.
   \param time Current time.
   \return Number of absorbed particles.
   \note Particle states must define:
   - <tt>Info info</tt>
   - \c time
*/
template <typename Subject>
std::size_t nr_absorbed(Subject const &subject, double time) {
  std::size_t absorbed = 0;
  for (auto const &part : subject) {
    if (ptof::absorbed(part, time))
      ++absorbed;
  }
  return absorbed;
}

/**
   \param subject CTRW object.
   \param time Current time.
   \return Number of adsorbed particles.
   \note Particle states must define:
   - <tt>Info info</tt>
   - \c time
*/
template <typename Subject>
std::size_t nr_adsorbed(Subject const &subject, double time) {
  std::size_t adsorbed = 0;
  for (auto const &part : subject) {
    if (ptof::adsorbed(part.state_old) && brackets_time(part, time)) {
      ++adsorbed;
    }
  }
  return adsorbed;
}

/**
   \param subject CTRW object.
   \param time Current time.
   \return Total mass.
   \note Particle states must define:
   - \c mass
   - \c time
*/
template <typename Subject> auto mass(Subject const &subject, double time) {
  double mass = 0.;
  for (auto const &part : subject) {
    if (!adsorbed(part.state_old()) && brackets_time(part, time)) {
      mass += ctrw::Get_interp{time, ctrw::Get_mass{}}(part.state_new(),
                                                       part.state_old());
    }
  }
  return mass;
}

/**
   \param subject CTRW object.
   \param time Current time.
   \return Total adsorbed mass.
   \note Particle states must define:
   - \c mass
   - \c time
*/
template <typename Subject>
auto mass_adsorbed(Subject const &subject, double time) {
  double mass = 0.;
  for (auto const &part : subject) {
    if (adsorbed(part.state_old()) && brackets_time(part, time)) {
      mass += ctrw::Get_interp{time, ctrw::Get_mass{}}(part.state_new(),
                                                       part.state_old());
    }
  }
  return mass;
}

/**
   \param subject CTRW object.
   \param time Current time.
   \return Total absorbed mass.
   \note Particle states must define:
   - \c mass
   - \c time
*/
template <typename Subject>
auto mass_absorbed(Subject const &subject, double time) {
  double mass = 0.;
  for (auto const &part : subject) {
    if (absorbed(part, time)) {
      mass += part.state_new().mass;
    }
  }
  return mass;
}

/**
   \param subject CTRW object.
   \param time Current time.
   \param masks Scalar fields.
   \param thresholds Thresholds for each mask.
   \return Total masses in regions where mask is above or equal to threshold.
   \note
   -Particle states must define:
       - mass
       - time
       - cell
   - \c thresholds must have at least the same size as \c masks
*/
template <typename Subject, typename Mask>
auto mass(Subject const &subject, double time,
          std::vector<std::reference_wrapper<const Mask>> masks,
          std::vector<double> thresholds) {
  std::vector<double> masses(masks.size(), 0.);
  for (auto const &part : subject) {
    auto const &state_old = part.state_old();
    if (adsorbed(state_old) || !brackets_time(part, time)) {
      continue;
    }
    auto const &state_new = part.state_new();
    auto cell_new = state_new.cell;
    auto cell_old = state_old.cell;
    if (outside(cell_new) && outside(cell_old)) {
      continue;
    }
    for (std::size_t ii = 0; ii < masks.size(); ++ii) {
      if (masks[ii].get()[cell_new] >= thresholds[ii] &&
          masks[ii].get()[cell_old] >= thresholds[ii]) {
        masses[ii] +=
            ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new, state_old);
      }
    }
  }
  return masses;
}

/**
   \param subject CTRW object.
   \param time Current time.
   \param getter_position Get position from state, gets position directly by
   default.
   \return Mean position (weighted by mass).
   \note Particle states must define:
   - \c position [for default position getter]
   - \c mass
   - \c time
*/
template <typename Subject, typename GetterPosition = ctrw::Get_position>
auto position_mean(Subject const &subject, double time,
                   GetterPosition getter_position = {}) {
  using Position = decltype(getter_position(subject.particles(0).state_new()));
  Position position_mean = Foam::zero{};
  for (auto const &part : subject) {
    if (!adsorbed(part.state_old()) && brackets_time(part, time)) {
      auto const &state_new = part.state_new();
      auto const &state_old = part.state_old();
      position_mean +=
          ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new, state_old) *
          ctrw::Get_interp{time, getter_position}(state_new, state_old);
    }
  }
  return position_mean / mass(subject, time);
}

/**
   \param subject CTRW object.
   \param time Current time.
   \param getter_position Get position from state, gets position directly by
   default.
   \return Second moment of position (weighted by mass).
   \note Particle states must define:
   - \c position [for default positition getter]
   - \c mass
*/
template <typename Subject, typename GetterPosition = ctrw::Get_position>
auto position_second_moment(Subject const &subject, double time,
                            GetterPosition getter_position = {}) {
  using Position = decltype(getter_position(subject.particles(0).state_new()));
  Position second_moment = Foam::zero{};
  for (auto const &part : subject) {
    if (!adsorbed(part.state_old()) && brackets_time(part, time)) {
      auto const &state_new = part.state_new();
      auto const &state_old = part.state_old();
      second_moment +=
          ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new, state_old) *
          op::square(
              ctrw::Get_interp{time, getter_position}(state_new, state_old));
    }
  }
  return second_moment / mass(subject, time);
}

/**
   \param subject CTRW object.
   \param exponents Order of moment or vector with order for each dimension.
   \param time Current time.
   \param getter_position Get position from state, gets \c position directly
   by default. \return Nth moment of position (weighted by mass). \note
   Particle states must define:
   - \c position [for default positition getter]
   - \c mass
*/
template <typename Subject, typename Exponents,
          typename GetterPosition = ctrw::Get_position>
auto position_moment(Subject const &subject, Exponents const &exponents,
                     double time, GetterPosition getter_position = {}) {
  decltype(getter_position(subject.particles(0).state_new())) moment =
      Foam::zero{};
  for (auto const &part : subject) {
    if (!adsorbed(part.state_old()) && brackets_time(part, time)) {
      auto const &state_new = part.state_new();
      auto const &state_old = part.state_old();
      moment +=
          ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new, state_old) *
          op::pow(ctrw::Get_interp{time, getter_position}(state_new, state_old),
                  exponents);
    }
  }
  return moment / mass(subject, time);
}

/**
   \param subject CTRW object.
   \param time Current time.
   \param getter_position Get position from state, gets position directly by
   default.
   \return Position variance (weighted by mass).
   \note Particle states must define:
   - \c position [for default positition getter]
   - \c mass
*/
template <typename Subject, typename GetterPosition = ctrw::Get_position>
auto position_variance(Subject const &subject, double time,
                       GetterPosition getter_position = {}) {
  return position_second_moment(subject, time, getter_position) -
         op::square(position_mean(subject, time, getter_position));
}

/**
   \param subject CTRW object.
   \param time Current time.
   \param field Field to evaluate at particle positions.
   \return Mean field value (weighted by mass).
   \note Particle states must define:
   - \c mass
   - \c time
*/
template <typename Subject, typename Field>
auto mean(Subject const &subject, double time, Field const &field) {
  decltype(field(subject.particles(0).state_new())) field_mean = Foam::zero{};
  for (auto const &part : subject) {
    auto const &state_old = part.state_old();
    auto const &state_new = part.state_new();
    if (!adsorbed(part.state_old()) && brackets_time(part, time)) {
      field_mean +=
          ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new, state_old) *
          ctrw::Get_interp{time, ctrw::Get_property{field}}(state_new,
                                                            state_old);
    }
  }
  return field_mean / mass(subject, time);
}

/**
   \brief Interpolate position between two particle states.
   \details Linearly interpolates between the two state positions according to
   \c time and state times. If interpolated position would be outside mesh, use
   nearest cell center.
   \param state_new New particle state.
   \param state_old Old particle state.
   \param time Time (between state times) to interpolate to.
   \param locator Object to locate positions in mesh.
   \param get_position Get position from state.
   \return Interpolated position.
   \note Particle states must define:
   - \c position
   - \c time
   - \c cell
*/
template <typename State, typename Locator,
          typename Getter = ctrw::Get_position>
auto interpolate_position(State const &state_new, State const &state_old,
                          double time, Locator const &locator,
                          Getter &&get_position = {}) {
  auto position =
      ctrw::Get_interp{time, ctrw::Get_position{}}(state_new, state_old);
  auto cell_id = locator(position, state_new.cell);
  if (outside(cell_id)) {
    auto nearest_cell_id = locator.nearest_cell(position);
    return State::make_position(cell_center(nearest_cell_id, locator.mesh()));
  }
  return position;
}

/**
   \brief Interpolate position between two particle states.
   \details Linearly interpolates between the two state positions according to
   \c time and state times. If interpolated position would be outside mesh, use
   nearest cell center.
   \param state_new New particle state.
   \param state_old Old particle state.
   \param time Time (between state times) to interpolate to.
   \param locator Object to locate positions in mesh.n

   \param get_position Get position from state.
   \return Interpolated position and corresponding cell in mesh.
   \note Particle states must define:
   - \c position
   - \c time
   - \c cell
*/
template <typename State, typename Locator,
          typename Getter = ctrw::Get_position>
auto interpolate_position_with_cell(State const &state_new,
                                    State const &state_old, double time,
                                    Locator const &locator,
                                    Getter &&get_position = {}) {
  auto position = ctrw::Get_interp{time, get_position}(state_new, state_old);
  auto cell_id = locator(position, state_new.cell);
  if (outside(cell_id)) {
    auto nearest_cell_id = locator.nearest_cell(position);
    return PositionAndCell{
        State::make_position(cell_center(nearest_cell_id, locator.mesh())),
        nearest_cell_id};
  }
  return PositionAndCell{position, cell_id};
};

/**
   \brief Interpolate position between two particle states.
   \details Linearly interpolates between the two state positions according to
   \c time and state times. If interpolated position would be outside mesh, use
   nearest cell center.
   \param particle Particle to interpolate between old and new state.
   \param time Time (between state times) to interpolate to.
   \param locator Object to locate positions in mesh.
   \param get_position Get position from state.
   \return Interpolated position.
   \note Particle states must define:
   - \c position
   - \c time
   - \c cell
*/
template <typename Particle, typename Locator,
          typename Getter = ctrw::Get_position>
auto interpolate_position(Particle const &particle, double time,
                          Locator const &locator, Getter &&get_position = {}) {
  return interpolate_position(particle.state_new(), particle.state_old(), time,
                              locator, get_position);
}

/**
   \brief Interpolate position between two particle states.
   \details Linearly interpolates between the two state positions according to
   \c time and state times. If interpolated position would be outside mesh, use
   nearest cell center.
   \param particle Particle to interpolate between old and new state.
   \param time Time (between state times) to interpolate to.
   \param locator Object to locate positions in mesh.n

   \param get_position Get position from state.
   \return Interpolated position and corresponding cell in mesh.
   \note Particle states must define:
   - \c position
   - \c time
   - \c cell
*/
template <typename Particle, typename Locator,
          typename Getter = ctrw::Get_position>
auto interpolate_position_with_cell(Particle const &particle, double time,
                                    Locator const &locator,
                                    Getter &&get_position = {}) {
  return interpolate_position_with_cell(
      particle.state_new(), particle.state_old(), time, locator, get_position);
};

/**
   \brief Compute demand-drived mesh data.
   \note
   - For parallel runs where each thread has access to the whole mesh,
   demand-driven mesh data must be pre-computed, otherwise there can be
   synchronization problems.
   - The cell tree used in searching is not precomputed. Mesh search tools
   should be used for searching instead.
   - Global mesh data is not precomputed. It assumed that the mesh is
   global (not decomposed).
   - Mesh-movement related quantities are not precomputed.
   - Not all mesh connectivity quantities are precomputed (only faceEdges and
   cellCells).
*/
template <typename Mesh> void compute_demand_driven_mesh_data(Mesh &mesh) {
  (void)mesh.V();
  (void)mesh.magSf();
  (void)mesh.C();
  (void)mesh.Cf();
  (void)mesh.faceEdges();
  (void)mesh.cellCells();
}

/**
   \brief Compute demand-drived meshSearch data.
   \note
   - For parallel runs where each thread has access to the whole mesh,
   demand-driven mesh data must be pre-computed, otherwise there can be
   synchronization problems.
   - Non-coupled boundary tree is not precomputed.
*/
template <typename MeshSearch>
void compute_demand_driven_meshSearch_data(MeshSearch &mesh_search) {
  (void)mesh_search.cellTree();
  (void)mesh_search.boundaryTree();
}

/**
   \brief Output time information.
   \param output Output stream.
   \param params Output parameters, holding time_unit_factor to rescale \p
   time. \param time Current time.
 */
template <typename ParametersOutput>
void info_time(std::ostream &output, ParametersOutput const &params,
               double time) {
  output << "Time " << "[" << params.time_units
         << " time units]: " << time / params.time_unit_factor << "\n";
}

/**
   \brief Output information about fraction of particles that have not been
   absorbed.
   \param output Output stream.
   \param subject CTRW object.
   \param time Current time.
   \note Particle states must define:
   - <tt>info.absorbed</tt> [std::size_t]
*/
template <typename Subject>
void info_fraction_not_absorbed(std::ostream &output, Subject const &subject,
                                double time) {
  output << "Fraction not absorbed: "
         << 1. - double(nr_absorbed(subject, time)) / subject.size() << "\n";
}

/** \return Identifier string for output file names. */
inline std::string identifier(std::string const &model_name,
                              std::string const &case_name,
                              std::string const &of_case_name,
                              std::string const &params_transport_name,
                              std::string const &params_reaction_name,
                              std::string const &params_solvers_name,
                              std::string const &params_initial_condition_name,
                              std::string const &params_output_name) {
  return "M_" + model_name + "_C_" + case_name + "_OF_" + of_case_name + "_T_" +
         params_transport_name + "_R_" + params_reaction_name + "_S_" +
         params_solvers_name + "_I_" + params_initial_condition_name + "_O_" +
         params_output_name;
}

/** \return Identifier string for output file names. */
inline std::string identifier(std::string const &model_name,
                              std::string const &case_name,
                              std::string const &of_case_name,
                              std::string const &params_transport_name,
                              std::string const &params_reaction_name,
                              std::string const &params_solvers_name,
                              std::string const &params_initial_condition_name,
                              std::string const &params_output_name,
                              std::size_t run_nr) {
  return identifier(model_name, case_name, of_case_name, params_transport_name,
                    params_reaction_name, params_solvers_name,
                    params_initial_condition_name, params_output_name) +
         "_RUN_" + std::to_string(run_nr);
}

/** \return Identifier string for output file names. */
inline std::string identifier(std::string const &model_name,
                              std::string const &case_name,
                              std::string const &of_case_name,
                              std::string const &params_transport_name,
                              std::string const &params_phase_name,
                              std::string const &params_reaction_name,
                              std::string const &params_solvers_name,
                              std::string const &params_initial_condition_name,
                              std::string const &params_output_name) {
  return "M_" + model_name + "_C_" + case_name + "_OF_" + of_case_name + "_T_" +
         params_transport_name + "_P_" + params_phase_name + "_R_" +
         params_reaction_name + "_S_" + params_solvers_name + "_I_" +
         params_initial_condition_name + "_O_" + params_output_name;
}

/** \return Identifier string for output file names. */
inline std::string identifier(std::string const &model_name,
                              std::string const &case_name,
                              std::string const &of_case_name,
                              std::string const &params_transport_name,
                              std::string const &params_phase_name,
                              std::string const &params_reaction_name,
                              std::string const &params_solvers_name,
                              std::string const &params_initial_condition_name,
                              std::string const &params_output_name,
                              std::size_t run_nr) {
  return identifier(model_name, case_name, of_case_name, params_transport_name,
                    params_phase_name, params_reaction_name,
                    params_solvers_name, params_initial_condition_name,
                    params_output_name) +
         "_RUN_" + std::to_string(run_nr);
}
} // namespace ptof

#endif /* PTOF_USEFUL_H */
