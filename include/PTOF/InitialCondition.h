/**
 \file PTOF/InitialConditions.h
 \author Tomás Aquino
 \date 24/02/2022
*/

#ifndef PTOF_INITIALCONDITIONS_H
#define PTOF_INITIALCONDITIONS_H

#include "CTRW/Meta.h"
#include "PTOF/Boundary.h"
#include "PTOF/Useful.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fieldTypes.H>
#include <map>
#include <numeric>
#include <point.H>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace ptof {
/** \struct InitialConditions PTOF/InitialCondition.h "PTOF/InitialCondition.h"
 * \brief Keep track of names and
 * types of initial conditions. */
struct InitialConditions {
  /** \enum Type
   *  \brief Implemented types. */
  enum class Type {
    uniform,      /**< Homogeneous throughout the domain.                     */
    fluxweighted, /**< Flux-weighted throughout the domain.                   */
    uniform_inlet,      /**< Homogeneous at the inlet.                        */
    fluxweighted_inlet, /**< Flux-weighted at the inlet.                      */
    uniform_solid,      /**< Homogeneous at the solid surface.                */
    uniform_near_solid, /**< Homogeneous at a fixed distance to the solid
                           interface.                                         */
    uniform_region_cartesian, /**< Homogeneous in a Cartesian region          */
    fluxweighted_region_cartesian, /**< Flux-weighted in a Cartesian region   */
    prescribed_positions,          /**< Prescribed particle positions         */
    prescribed_positions_masses,   /**< Prescribed particle positions and masses
                                    */
    prescribed_positions_masses_tags, /**< Prescribed particle positions,
                                         masses, and tags                     */
    uniform_inlet_continuous,      /**< Continuous injection homogeneous at the
                                      inlet.                                  */
    fluxweighted_inlet_continuous, /**< Continuous injection flux-weighted at
                                      the inlet.                              */
  };

  /** Type from name. */
  static auto type(std::string const &name) { return string_to_type.at(name); }

  /** Name from name. */
  static auto name(Type type) { return type_to_string.at(type); }

  /** Check if name exists. */
  static auto contains(std::string const &name) {
    return string_to_type.count(name);
  }

  /** Map names to types. */
  inline static const std::map<std::string, Type> string_to_type{
      {"uniform", Type::uniform},
      {"fluxweighted", Type::fluxweighted},
      {"uniform_inlet", Type::uniform_inlet},
      {"fluxweighted_inlet", Type::fluxweighted_inlet},
      {"uniform_solid", Type::uniform_solid},
      {"uniform_near_solid", Type::uniform_near_solid},
      {"uniform_region_cartesian", Type::uniform_region_cartesian},
      {"fluxweighted_region_cartesian", Type::fluxweighted_region_cartesian},
      {"prescribed_positions", Type::prescribed_positions},
      {"prescribed_positions_masses", Type::prescribed_positions_masses},
      {"prescribed_positions_masses_tags",
       Type::prescribed_positions_masses_tags},
      {"uniform_inlet_continuous", Type::uniform_inlet_continuous},
      {"fluxweighted_inlet_continuous", Type::fluxweighted_inlet_continuous}};

  /** Map types to names. */
  inline static const std::map<Type, std::string> type_to_string{
      {
          Type::uniform,
          "uniform",
      },
      {Type::fluxweighted, "fluxweighted"},
      {Type::uniform_inlet, "uniform_inlet"},
      {Type::fluxweighted_inlet, "fluxweighted_inlet"},
      {Type::uniform_solid, "uniform_solid"},
      {Type::uniform_near_solid, "uniform_near_solid"},
      {Type::uniform_region_cartesian, "uniform_region_cartesian"},
      {Type::fluxweighted_region_cartesian, "fluxweighted_region_cartesian"},
      {Type::prescribed_positions, "prescribed_positions"},
      {Type::prescribed_positions_masses, "prescribed_positions_masses"},
      {Type::prescribed_positions_masses_tags,
       "prescribed_positions_masses_tags"},
      {Type::uniform_inlet_continuous, "uniform_inlet_continuous"},
      {Type::fluxweighted_inlet_continuous, "fluxweighted_inlet_continuous"}};
};

/** \brief Remove cell/face indices if mask is above threshold.
 \param cell_or_face_ids Container with mesh cell or face indices.
 \param mask Scalar field assigning values to mesh cell indices through
 operator[]. \param threshold Threshold for the mask, such that cells where the
 mask is above the threshold are considered. */
template <typename Container, typename Mask>
void apply_mask_in_place(Container &cell_or_face_ids, Mask const &mask,
                         double threshold = 0.) {
  useful::swap_erase_if(cell_or_face_ids, [&mask, threshold](auto ii) {
    return mask[ii] > threshold;
  });
}

/** \brief Remove cell/face indices if mask is above threshold.
\param cell_or_face_ids Container with mesh cell or face indices.
\param mask Scalar field assigning values to mesh cell indices through
operator[]. \param threshold Threshold for the mask, such that cells where the
mask is above the threshold are considered. \return \p cell_or_face_ids with
indices of elements where mask is above threshold removed. */
template <typename Container, typename Mask>
auto apply_mask(Container const &cell_or_face_ids, Mask const &mask,
                double threshold = 0.) {
  auto masked_ids = cell_or_face_ids;
  apply_mask_in_place(masked_ids, mask, threshold);

  return masked_ids;
}

/** \brief Make particles distributed randomly uniformly over given cells.
 \details Place particles at cell centers.
 \param nr_particles Number of particles to make.
 \param cell_ids Mesh cell indices.
 \param mesh Mesh object.
 \param rng Random number generator.
 \param particle_maker Functor to make a particle given a position.
 \return Container with particles. */
template <typename Container, typename Mesh, typename RNG,
          typename ParticleMaker>
auto uniform_cells(std::size_t nr_particles, Container const &cell_ids,
                   Mesh const &mesh, RNG &rng, ParticleMaker &particle_maker) {
  if (cell_ids.size() == 0 && nr_particles != 0)
    throw std::runtime_error{"No available cells to place particles"};

  // Weigh probabilities with cell volumes.
  std::vector<double> weights;
  weights.reserve(cell_ids.size());
  for (auto cell : cell_ids)
    weights.push_back(cell_volume(cell, mesh));
  std::discrete_distribution<Foam::label> dist{weights.begin(), weights.end()};

  using Particle = decltype(particle_maker(Foam::point{}));
  std::vector<Particle> particles;
  particles.reserve(nr_particles);
  for (std::size_t pp = 0; pp < nr_particles; ++pp)
    particles.push_back(particle_maker(cell_center(cell_ids[dist(rng)], mesh)));
  return particles;
}

/** \brief Make particles distributed over given cells with probability
proportional to local velocity magnitude. \details Place particles at cell
centers. \param nr_particles Number of particles to make. \param cell_ids Mesh
cell indices. \param velocity_field Velocity field as a function of state.
\param mesh Mesh object.
\param rng Random number generator.
\param particle_maker Functor to make a particle given a position.
\return Container with particles. */
template <typename Container, typename VelocityField, typename Mesh,
          typename RNG, typename ParticleMaker>
auto fluxweighted_cells(std::size_t nr_particles, Container const &cell_ids,
                        VelocityField const &velocity_field, Mesh const &mesh,
                        RNG &rng, ParticleMaker &particle_maker) {
  // Weigh probabilities with cell volumes
  // and velocity magnitudes at center.
  if (cell_ids.size() == 0 && nr_particles != 0)
    throw std::runtime_error{"No available cells to place particles"};

  std::vector<double> weights;
  weights.reserve(cell_ids.size());
  for (auto cell : cell_ids)
    weights.push_back(cell_volume(cell, mesh) *
                      Foam::mag(velocity_field(cell_center(cell, mesh), cell)));
  std::discrete_distribution<Foam::label> dist{weights.begin(), weights.end()};

  using Particle = decltype(particle_maker(Foam::point{}));
  std::vector<Particle> particles;
  particles.reserve(nr_particles);
  for (std::size_t pp = 0; pp < nr_particles; ++pp)
    particles.push_back(particle_maker(cell_center(cell_ids[dist(rng)], mesh)));
  return particles;
}

/** Make particles randomly uniformly
 *  over given mesh cell faces.
 *  Place particles at center of faces. */
/** \brief Make particles distributed randomly uniformly over given mesh cell
faces. \details Place particles at cell centers. \param nr_particles Number of
particles to make. \param face_ids Mesh face indices. \param locator Object to
locate positions in mesh. \param rng Random number generator. \param
particle_maker Functor to make a particle given a position. \return Container
with particles. */
template <typename Container, typename Locator, typename RNG,
          typename ParticleMaker>
auto uniform_faces(std::size_t nr_particles, Container const &face_ids,
                   Locator const &locator, RNG &rng,
                   ParticleMaker &particle_maker) {
  if (face_ids.size() == 0 && nr_particles != 0)
    throw std::runtime_error{"No available faces to place particles"};

  auto const &mesh = locator.mesh();

  std::vector<double> weights;
  weights.reserve(face_ids.size());
  for (auto face : face_ids)
    weights.push_back(face_area(face, mesh));
  std::discrete_distribution<Foam::label> dist{weights.begin(), weights.end()};

  using Particle = decltype(particle_maker(Foam::point{}));
  std::vector<Particle> particles;
  particles.reserve(nr_particles);
  for (std::size_t pp = 0; pp < nr_particles; ++pp) {
    auto face = face_ids[dist(rng)];
    particles.push_back(particle_maker(
        adjusted_face_center(face, locator, mesh.faceOwner()[face])));
  }
  return particles;
}

/** \brief Make particles distributed over given cell faces with probability
proportional to local velocity magnitude. \details Place particles at face
centers. \param nr_particles Number of particles to make. \param face_ids Mesh
face indices. \param velocity_field Velocity field as a function of state.
\param locator Object to locate positions in mesh.
\param rng Random number generator.
\param particle_maker Functor to make a particle given a position.
\return Container with particles. */
template <typename Container, typename VelocityField, typename Locator,
          typename RNG, typename ParticleMaker>
auto fluxweighted_faces(std::size_t nr_particles, Container const &face_ids,
                        VelocityField const &velocity_field,
                        Locator const &locator, RNG &rng,
                        ParticleMaker &particle_maker) {
  if (face_ids.size() == 0 && nr_particles != 0)
    throw std::runtime_error{"No available faces to place particles"};

  auto const &mesh = locator.mesh();

  std::vector<double> weights;
  weights.reserve(face_ids.size());
  for (auto face : face_ids)
    weights.push_back(face_area(face, mesh) *
                      Foam::mag(velocity_field(adjusted_face_center(
                          face, locator, mesh.faceOwner()[face]))));
  std::discrete_distribution<Foam::label> dist{weights.begin(), weights.end()};

  using Particle = decltype(particle_maker(Foam::point{}));
  std::vector<Particle> particles;
  particles.reserve(nr_particles);
  for (std::size_t pp = 0; pp < nr_particles; ++pp) {
    auto face = face_ids[dist(rng)];
    particles.push_back(particle_maker(
        adjusted_face_center(face, locator, mesh.faceOwner()[face])));
  }
  return particles;
}

/** Make particles randomly uniformly
 *  at given distance to given mesh cell face centers. */
/** \brief Make particles distributed randomly uniformly near given cell face
centers. \param nr_particles Number of particles to make. \param face_ids Mesh
face indices. \param distance Distance from face centers. \param locator Object
to locate positions in mesh. \param rng Random number generator. \param
particle_maker Functor to make a particle given a position. \return Container
with particles. */
template <typename Container, typename Locator, typename RNG,
          typename ParticleMaker>
auto uniform_near_boundary_faces(std::size_t nr_particles,
                                 Container const &face_ids, double distance,
                                 Locator const &locator, RNG &rng,
                                 ParticleMaker &particle_maker) {
  if (face_ids.size() == 0 && nr_particles != 0)
    throw std::runtime_error{"No available faces to place particles"};

  auto const &mesh = locator.mesh();

  std::vector<double> weights;
  weights.reserve(face_ids.size());
  for (auto face : face_ids)
    weights.push_back(face_area(face, mesh));
  std::discrete_distribution<Foam::label> dist{weights.begin(), weights.end()};

  using Particle = decltype(particle_maker(Foam::point{}));
  std::vector<Particle> particles;
  particles.reserve(nr_particles);
  for (std::size_t pp = 0; pp < nr_particles; ++pp) {
    auto rand = dist(rng);
    auto unit_normal = unit_normal_inward(face_ids[rand], mesh);
    auto position = face_center(face_ids[rand], mesh);
    +distance *unit_normal;
    if (outside(locator.mesh_search().findCell(position)) ||
        locator.mesh_search().findNearestBoundaryFace(position) !=
            face_ids[rand])
      continue;
    particles.push_back(particle_maker(position));
  }
  return particles;
}

/** \brief Make particles distributed randomly uniformly over all mesh cells.
\details Place particles at cell centers.
\param nr_particles Number of particles to make.
\param mesh Mesh object.
\param rng Random number generator.
\param particle_maker Functor to make a particle given a position.
\return Container with particles. */
template <typename Mesh, typename RNG, typename ParticleMaker>
auto uniform(std::size_t nr_particles, Mesh const &mesh, RNG &rng,
             ParticleMaker &particle_maker) {
  return uniform_cells(nr_particles, all_cell_ids(mesh), mesh, rng,
                       particle_maker);
}

/** \brief Make particles distributed randomly uniformly over mesh cells where
 mask is above threshold. \details Place particles at cell centers. \param
 nr_particles Number of particles to make. \param mesh Mesh object. \param rng
 Random number generator. \param particle_maker Functor to make a particle given
 a position. \param mask Scalar field assigning values to mesh cell indices
 through operator[]. \param threshold Threshold for the mask, such that cells
 where the mask is above the threshold are considered. \return Container with
 particles. */
template <typename Mesh, typename RNG, typename ParticleMaker, typename Mask>
auto uniform(std::size_t nr_particles, Mesh const &mesh, RNG &rng,
             ParticleMaker &particle_maker, Mask const &mask,
             double threshold = 0.) {
  return uniform_cells(nr_particles,
                       apply_mask(all_cell_ids(mesh), mask, threshold), mesh,
                       rng, particle_maker);
}

/** \brief Make particles distributed with probability proportional to local
velocity magnitude over all mesh cells. \details Place particles at cell
centers. \param nr_particles Number of particles to make. \param velocity_field
Velocity field as a function of state. \param mesh Mesh object. \param rng
Random number generator. \param particle_maker Functor to make a particle given
a position. \return Container with particles. */
template <typename VelocityField, typename Mesh, typename RNG,
          typename ParticleMaker>
auto fluxweighted(std::size_t nr_particles, VelocityField const &velocity_field,
                  Mesh const &mesh, RNG &rng, ParticleMaker &particle_maker) {
  return fluxweighted_cells(nr_particles, all_cell_ids(mesh), velocity_field,
                            mesh, rng, particle_maker);
}

/** \brief Make particles distributed with probability proportional to local
velocity magnitude over mesh cells where mask is above threshold. \details Place
particles at cell centers. \param nr_particles Number of particles to make.
\param velocity_field Velocity field as a function of state.
\param mesh Mesh object.
\param rng Random number generator.
\param particle_maker Functor to make a particle given a position.
\param mask Scalar field assigning values to mesh cell indices through
operator[]. \param threshold Threshold for the mask, such that cells where the
mask is above the threshold are considered. \return Container with particles. */
template <typename VelocityField, typename Mesh, typename RNG,
          typename ParticleMaker, typename Mask>
auto fluxweighted(std::size_t nr_particles, VelocityField const &velocity_field,
                  Mesh const &mesh, RNG &rng, ParticleMaker &particle_maker,
                  Mask const &mask, double threshold = 0.) {
  return fluxweighted_cells(nr_particles,
                            apply_mask(all_cell_ids(mesh), mask, threshold),
                            velocity_field, mesh, rng, particle_maker);
}

/** \brief Make particles distributed randomly uniformly over cell faces in
given patches. \details Place particles at face centers. \param nr_particles
Number of particles to make. \param patch_names Names of patches. \param locator
Object to locate positions in mesh. \param rng Random number generator. \param
particle_maker Functor to make a particle given a position. \return Container
with particles. */
template <typename Locator, typename RNG, typename ParticleMaker>
auto uniform_patches(std::size_t nr_particles,
                     std::vector<std::string> const &patch_names,
                     Locator const &locator, RNG &rng,
                     ParticleMaker &particle_maker) {
  return uniform_faces(nr_particles,
                       patches_face_ids(locator.mesh(), patch_names), locator,
                       rng, particle_maker);
}

/** \brief Make particles distributed randomly uniformly over cell faces in
given patches where mask is above threshold. \details Place particles at face
centers. \param nr_particles Number of particles to make. \param patch_names
Names of patches. \param locator Object to locate positions in mesh. \param rng
Random number generator. \param particle_maker Functor to make a particle given
a position. \param mask Scalar field assigning values to mesh cell indices
through operator[]. \param threshold Threshold for the mask, such that cells
where the mask is above the threshold are considered. \return Container with
particles. */
template <typename Locator, typename RNG, typename ParticleMaker, typename Mask>
auto uniform_patches(std::size_t nr_particles,
                     std::vector<std::string> const &patch_names,
                     Locator const &locator, RNG &rng,
                     ParticleMaker &particle_maker, Mask const &mask,
                     double threshold = 0.) {
  return uniform_faces(nr_particles,
                       apply_mask(patches_face_ids(locator.mesh(), patch_names),
                                  mask, threshold),
                       locator, rng, particle_maker);
}

/** \brief Make particles distributed with probability proportional to local
velocity magnitude over cell faces in given patches. \details Place particles at
face centers. \param nr_particles Number of particles to make. \param
patch_names Names of patches. \param velocity_field Velocity field as a function
of state. \param locator Object to locate positions in mesh. \param rng Random
number generator. \param particle_maker Functor to make a particle given a
position. \return Container with particles. */
template <typename VelocityField, typename Locator, typename RNG,
          typename ParticleMaker>
auto fluxweighted_patches(std::size_t nr_particles,
                          std::vector<std::string> const &patch_names,
                          VelocityField const &velocity_field,
                          Locator const &locator, RNG &rng,
                          ParticleMaker &particle_maker) {
  return fluxweighted_faces(nr_particles,
                            patches_face_ids(locator.mesh(), patch_names),
                            velocity_field, locator, rng, particle_maker);
}

/** \brief Make particles distributed with probability proportional to local
velocity magnitude over cell faces in given patches where mask is above
threshold. \details Place particles at face centers. \param nr_particles Number
of particles to make. \param patch_names Names of patches. \param velocity_field
Velocity field as a function of state. \param locator Object to locate positions
in mesh. \param rng Random number generator. \param particle_maker Functor to
make a particle given a position. \param mask Scalar field assigning values to
mesh cell indices through operator[]. \param threshold Threshold for the mask,
such that cells where the mask is above the threshold are considered. \return
Container with particles. */
template <typename VelocityField, typename Locator, typename RNG,
          typename ParticleMaker, typename Mask>
auto fluxweighted_patches(std::size_t nr_particles,
                          std::vector<std::string> const &patch_names,
                          VelocityField const &velocity_field,
                          Locator const &locator, RNG &rng,
                          ParticleMaker &particle_maker, Mask const &mask,
                          double threshold = 0.) {
  return fluxweighted_faces(
      nr_particles,
      apply_mask(patches_face_ids(locator.mesh(), patch_names), mask,
                 threshold),
      velocity_field, locator, rng, particle_maker);
}

/** \brief Make particles distributed randomly uniformly near cell faces in
given patches. \param nr_particles Number of particles to make. \param
patch_names Names of patches. \param distance Distance from face centers. \param
locator Object to locate positions in mesh. \param rng Random number generator.
\param particle_maker Functor to make a particle given a position.
\return Container with particles. */
template <typename Locator, typename RNG, typename ParticleMaker>
auto uniform_near_boundary_patches(std::size_t nr_particles,
                                   std::vector<std::string> const &patch_names,
                                   double distance, Locator const &locator,
                                   RNG &rng, ParticleMaker particle_maker) {
  return uniform_near_boundary_faces(
      nr_particles, patches_face_ids(locator.mesh(), patch_names), distance,
      locator, rng, particle_maker);
}

/** \brief Make particles distributed randomly uniformly near cell faces in
given patches where mask is above threshold. \param nr_particles Number of
particles to make. \param patch_names Names of patches. \param distance Distance
from face centers. \param locator Object to locate positions in mesh. \param rng
Random number generator. \param particle_maker Functor to make a particle given
a position. \param mask Scalar field assigning values to mesh cell indices
through operator[]. \param threshold Threshold for the mask, such that cells
where the mask is above the threshold are considered. \return Container with
particles. */
template <typename Locator, typename RNG, typename ParticleMaker, typename Mask>
auto uniform_near_boundary_patches(std::size_t nr_particles,
                                   std::vector<std::string> const &patch_names,
                                   double distance, Locator const &locator,
                                   RNG &rng, ParticleMaker particle_maker,
                                   Mask const &mask, double threshold = 0.) {
  return uniform_near_boundary_faces(
      nr_particles,
      apply_mask(patches_face_ids(locator.mesh(), patch_names), mask,
                 threshold),
      distance, locator, rng, particle_maker);
}

/** \brief Maker particles over a time window according to given pulse
properties. \details All particles are made at once but with different initial
times. Particles are injected at times ii \c time_step with ii from 0 to
((time_max - time_min)/time_step - 1). \param pulse Functor that returns a
vector of particles given number of particles and particle maker to make a
particle given position. \param nr_particles Number of particles to make. \param
time_min Injection start time. \param time_step Injection discretization time
step. \param time_max Injection stop time. \param particle_maker Functor to make
a particle given a position. \return Container with particles. */
template <typename Pulse, typename ParticleMaker>
auto continuous_injection(Pulse const &pulse, std::size_t nr_particles,
                          double time_min, double time_step, double time_max,
                          ParticleMaker &particle_maker) {
  std::size_t nr_injections = (time_max - time_min) / time_step;

  using Particle = decltype(particle_maker(Foam::point{}));
  std::vector<Particle> particles;
  particles.reserve(nr_injections * nr_particles);

  for (std::size_t ii = 0; ii < nr_injections; ++ii) {
    particle_maker.time = ii * time_step;
    for (auto const &part : pulse(nr_particles, particle_maker))
      particles.push_back(part);
  }

  return particles;
}

/** \brief Verify if initial condition is implemented.
 \param initial_condition Name of initial condition.
 \param implemented Associative container of implemented initial condition
 names. */
template <typename Implemented>
void verify_initial_condition(std::string const &initial_condition,
                              Implemented const &implemented) {
  if (!implemented.contains(initial_condition))
    throw std::runtime_error{"Initial condition type " + initial_condition +
                             " not supported"};
}

/** \param position_data Container of position vectors for each particle
 \param particle_maker Functor to make a particle given a position.
 \return Container with particles. */
template <typename Positions, typename ParticleMaker>
auto prescribed_positions(Positions const &position_data,
                          ParticleMaker &particle_maker) {
  using Particle = decltype(particle_maker(Foam::point{}));
  std::vector<Particle> particles;
  particles.reserve(position_data.size());
  for (std::size_t pp = 0; pp < position_data.size(); ++pp)
    particles.push_back(particle_maker(make_point(position_data[pp])));

  return particles;
}

/** \param nr_particles Number of particles to make.
\param filename_position_data Name of file with one position per line; the first
\c nr_particles positions are used. \param particle_maker Functor to make a
particle given a position. \return Container with particles. */
template <typename ParticleMaker>
auto prescribed_positions(std::size_t nr_particles,
                          std::string const &filename_data,
                          ParticleMaker &particle_maker) {
  std::string comment_sequence = "#";
  auto input = useful::open_read(filename_data);
  std::vector<std::vector<double>> position_data;
  position_data.reserve(nr_particles);
  for (std::size_t pp = 0; pp < nr_particles; ++pp) {
    std::string line;
    while (std::getline(input, line))
      if (line.find(comment_sequence) != 0)
        break;
    useful::remove_carriage_return_in_place(line);
    useful::clear_escape_in_place(line, comment_sequence);
    if (line.empty())
      throw std::runtime_error{
          std::string{"Could not parse prescribed position for "} +
          std::to_string(pp) + "th particle"};
    auto split_line = useful::split(line);
    position_data.emplace_back();
    position_data.back().reserve(split_line.size());
    for (auto const &position_component : split_line)
      position_data.back().push_back(std::stod(position_component));
  }
  input.close();

  return prescribed_positions(position_data, particle_maker);
}

/**
\param filename_data Name of file with one position per line.
\param particle_maker Functor to make a particle given a position.
\return Container with particles. */
template <typename ParticleMaker>
auto prescribed_positions(std::string const &filename_data,
                          ParticleMaker &particle_maker) {
  std::string comment_sequence = "#";
  auto input = useful::open_read(filename_data);
  std::vector<std::vector<double>> position_data;
  std::string line;
  while (std::getline(input, line)) {
    if (line.find(comment_sequence) == 0)
      continue;
    useful::remove_carriage_return_in_place(line);
    useful::clear_escape_in_place(line, comment_sequence);
    auto split_line = useful::split(line);
    position_data.emplace_back();
    position_data.back().reserve(split_line.size());
    for (auto const &position_component : split_line)
      position_data.back().push_back(std::stod(position_component));
  }
  input.close();

  return prescribed_positions(position_data, particle_maker);
}

/** \param position_data Container of position vectors for each particle
    \param mass_data Container of masses for each particle
    \param particle_maker Functor to make a particle given a position.
    \return Container with particles. */
template <typename Positions, typename Masses, typename ParticleMaker>
auto prescribed_positions_masses(Positions const &position_data,
                                 Masses const &mass_data,
                                 ParticleMaker &particle_maker) {
  using Particle = decltype(particle_maker(Foam::point{}));
  std::vector<Particle> particles;
  particles.reserve(position_data.size());
  for (std::size_t pp = 0; pp < position_data.size(); ++pp) {
    particles.push_back(particle_maker(make_point(position_data[pp])));
    if constexpr (meta::has_mass_v<ParticleMaker>)
      particle_maker.mass = mass_data[pp];
    else
      throw std::runtime_error{
          "prescribed_positions_masses : ParticleMaker does not define mass"};
  }

  return particles;
}

/** \param nr_particles Number of particles to make.
\param filename_data Name of file with one position and one mass per line; the
first \c nr_particles lines are used. \param particle_maker Functor to make a
particle given a position. \return Container with particles. */
template <typename ParticleMaker>
auto prescribed_positions_masses(std::size_t nr_particles,
                                 std::string const &filename_position_mass_data,
                                 ParticleMaker &particle_maker) {
  std::string comment_sequence = "#";
  auto input = useful::open_read(filename_position_mass_data);
  std::vector<std::vector<double>> position_data;
  position_data.reserve(nr_particles);
  std::vector<double> mass_data;
  mass_data.reserve(nr_particles);
  for (std::size_t pp = 0; pp < nr_particles; ++pp) {
    std::string line;
    while (std::getline(input, line))
      if (line.find(comment_sequence) != 0)
        break;
    useful::remove_carriage_return_in_place(line);
    useful::clear_escape_in_place(line, comment_sequence);
    if (line.empty())
      throw std::runtime_error{
          std::string{"Could not parse prescribed position for "} +
          std::to_string(pp) + "th particle"};
    auto split_line = useful::split(line);
    position_data.emplace_back();
    position_data.back().reserve(split_line.size());
    for (std::size_t ii = 0; ii < split_line.size() - 1; ++ii)
      position_data.back().push_back(std::stod(split_line[ii]));
    mass_data.push_back(std::stod(split_line.back()));
  }
  input.close();

  return prescribed_positions_masses(position_data, mass_data, particle_maker);
}

/**
\param filename_data Name of file with one position and one mass per line.
\param particle_maker Functor to make a particle given a position.
\return Container with particles. */
template <typename ParticleMaker>
auto prescribed_positions_masses(std::string const &filename_data,
                                 ParticleMaker &particle_maker) {
  std::string comment_sequence = "#";
  auto input = useful::open_read(filename_data);
  std::vector<std::vector<double>> position_data;
  std::vector<double> mass_data;
  std::string line;
  while (std::getline(input, line)) {
    if (line.find(comment_sequence) == 0)
      continue;
    useful::remove_carriage_return_in_place(line);
    useful::clear_escape_in_place(line, comment_sequence);
    auto split_line = useful::split(line);
    position_data.emplace_back();
    position_data.back().reserve(split_line.size());
    for (std::size_t ii = 0; ii < split_line.size() - 1; ++ii)
      position_data.back().push_back(std::stod(split_line[ii]));
    mass_data.push_back(std::stod(split_line.back()));
  }
  input.close();

  return prescribed_positions_masses(position_data, mass_data, particle_maker);
}

/** \param position_data Container of position vectors for each particle
    \param mass_data Container of masses for each particle
    \param tag_data Container of tags (non-negative integer identifiers) for
   each particle \param particle_maker Functor to make a particle given a
   position. \return Container with particles. */
template <typename Positions, typename Masses, typename Tags,
          typename ParticleMaker>
auto prescribed_positions_masses_tags(Positions const &position_data,
                                      Masses const &mass_data,
                                      Tags const &tag_data,
                                      ParticleMaker &particle_maker) {
  using Particle = decltype(particle_maker(Foam::point{}));
  std::vector<Particle> particles;
  particles.reserve(position_data.size());
  for (std::size_t pp = 0; pp < position_data.size(); ++pp) {
    particles.push_back(particle_maker(make_point(position_data[pp])));
    if constexpr (meta::has_mass_v<ParticleMaker>)
      particle_maker.mass = mass_data[pp];
    else
      throw std::runtime_error{"prescribed_positions_masses_tags : "
                               "ParticleMaker does not define mass"};
    if constexpr (meta::has_tag_v<ParticleMaker>)
      particle_maker.tag = tag_data[pp];
    else
      throw std::runtime_error{"prescribed_positions_masses_tags : "
                               "ParticleMaker does not define tag"};
  }

  return particles;
}

/** \param nr_particles Number of particles to make.
\param filename_data Name of file with one position, one mass, and one tag per
line; the first \c nr_particles lines are used. \param particle_maker Functor to
make a particle given a position. \return Container with particles. */
template <typename ParticleMaker>
auto prescribed_positions_masses_tags(
    std::size_t nr_particles, std::string const &filename_position_mass_data,
    ParticleMaker &particle_maker) {
  std::string comment_sequence = "#";
  auto input = useful::open_read(filename_position_mass_data);
  std::vector<std::vector<double>> position_data;
  position_data.reserve(nr_particles);
  std::vector<double> mass_data;
  mass_data.reserve(nr_particles);
  std::vector<std::size_t> tag_data;
  tag_data.reserve(nr_particles);
  for (std::size_t pp = 0; pp < nr_particles; ++pp) {
    std::string line;
    while (std::getline(input, line))
      if (line.find(comment_sequence) != 0)
        break;
    useful::remove_carriage_return_in_place(line);
    useful::clear_escape_in_place(line, comment_sequence);
    if (line.empty())
      throw std::runtime_error{
          std::string{"Could not parse prescribed position for "} +
          std::to_string(pp) + "th particle"};
    auto split_line = useful::split(line);
    position_data.emplace_back();
    position_data.back().reserve(split_line.size());
    for (std::size_t ii = 0; ii < split_line.size() - 2; ++ii)
      position_data.back().push_back(std::stod(split_line[ii]));
    mass_data.push_back(std::stod(split_line[split_line.size() - 2]));
    tag_data.push_back(std::stod(split_line.back()));
  }
  input.close();

  return prescribed_positions_masses(position_data, mass_data, particle_maker);
}

/**
\param filename_data Name of file with one position, one mass, and one tag per
line. \param particle_maker Functor to make a particle given a position. \return
Container with particles. */
template <typename ParticleMaker>
auto prescribed_positions_masses_tags(std::string const &filename_data,
                                      ParticleMaker &particle_maker) {
  std::string comment_sequence = "#";
  auto input = useful::open_read(filename_data);
  std::vector<std::vector<double>> position_data;
  std::vector<double> mass_data;
  std::vector<std::size_t> tag_data;
  std::string line;
  while (std::getline(input, line)) {
    if (line.find(comment_sequence) == 0)
      continue;
    useful::remove_carriage_return_in_place(line);
    useful::clear_escape_in_place(line, comment_sequence);
    auto split_line = useful::split(line);
    position_data.emplace_back();
    position_data.back().reserve(split_line.size());
    for (std::size_t ii = 0; ii < split_line.size() - 2; ++ii)
      position_data.back().push_back(std::stod(split_line[ii]));
    mass_data.push_back(std::stod(split_line[split_line.size() - 2]));
    tag_data.push_back(std::stod(split_line.back()));
  }
  input.close();

  return prescribed_positions_masses(position_data, mass_data, particle_maker);
}

/** \class InitialConditionParameters_Cases PTOF/InitialCondition.h
 * "PTOF/InitialCondition.h" \brief Initial condition condition parameters to
 * handle all initial condition types in \c InitialCondition_Cases. */
struct InitialConditionParameters_Cases {
private:
  std::string comment_sequence =
      "#"; /**< Sequence of characters marking comment for file parsing. */

public:
  InitialConditions::Type type;
  double distance_wall;
  std::vector<std::pair<double, double>> region_boundaries;
  std::string filename_data;
  double initial_mass;
  double time_min;
  double time_step_accuracy_adv;
  double time_step_accuracy_diff;
  double time_step_accuracy_react;
  double time_max;
  double time_step;

  template <typename Geometry, typename TransportParameters,
            typename ReactionParameters>
  InitialConditionParameters_Cases(Directories const &directories,
                                   std::string const &name,
                                   Geometry const &geometry,
                                   TransportParameters const &params_transport,
                                   ReactionParameters const &params_reaction) {
    std::string ic_name;
    auto input =
        useful::open_read(directories.dir_parameters +
                          "/parameters_initial_condition_" + name + ".dat");
    useful::read_first_from_line(input, ic_name, comment_sequence);
    verify_initial_condition(ic_name, InitialConditions{});
    type = InitialConditions::type(ic_name);
    if (type == InitialConditions::Type::uniform_near_solid)
      useful::read_first_from_line(input, distance_wall, comment_sequence);
    if (type == InitialConditions::Type::uniform_region_cartesian ||
        type == InitialConditions::Type::fluxweighted_region_cartesian) {
      std::string line;
      while (std::getline(input, line))
        if (line.find(comment_sequence) != 0)
          break;
      useful::remove_carriage_return_in_place(line);
      useful::clear_escape_in_place(line, comment_sequence);
      auto split_line = useful::split(line);
      if (split_line.size() % 2 != 0)
        std::runtime_error{
            "Initial condition region boundaries must come in pairs"};
      for (std::size_t ii = 0; ii < split_line.size(); ii += 2) {
        region_boundaries.push_back(
            {std::stod(split_line[ii]), std::stod(split_line[ii + 1])});
      }
    }
    if (type == InitialConditions::Type::prescribed_positions) {
      useful::read_first_from_line(input, filename_data, comment_sequence);
      useful::expand_env_in_place(
          useful::expand_home_dir_in_place(filename_data));
      std::cout << filename_data << std::endl;
    }
    if (type != InitialConditions::Type::prescribed_positions_masses &&
        type == InitialConditions::Type::prescribed_positions_masses_tags)
      useful::read_first_from_line(input, initial_mass, comment_sequence);
    useful::read_first_from_line(input, time_min, comment_sequence);
    if (type == InitialConditions::Type::uniform_inlet_continuous ||
        type == InitialConditions::Type::fluxweighted_inlet_continuous) {
      useful::read_first_from_line(input, time_max, comment_sequence);
      useful::read_first_from_line(input, time_step_accuracy_adv,
                                   comment_sequence);
      useful::read_first_from_line(input, time_step_accuracy_diff,
                                   comment_sequence);
      useful::read_first_from_line(input, time_step_accuracy_react,
                                   comment_sequence);
      time_step =
          std::min({time_step_accuracy_adv * params_transport.advection_time,
                    time_step_accuracy_diff * params_transport.diffusion_time,
                    time_step_accuracy_react * params_reaction.reaction_time});
    }
  }

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Initial condition parameters\n"
           "--------------------------------------------------------------\n"
           "- Initial condition type:\n"
           "\tuniform: Homogeneous throughout the domain\n"
           "\tfluxweighted: Flux-weighted throughout the domain\n"
           "\tuniform_inlet: Homogeneous at the inlet\n"
           "\tfluxweighted_inlet: Flux-weighted at the inlet\n"
           "\tuniform_solid: Homogeneous at the solid surface\n"
           "\tuniform_near_solid: Homogeneous at a fixed distance to the\n"
           "\t                    solid interface\n"
           "\tuniform_region_cartesian: Homogeneous in a prescribed\n"
           "\t                          cartesian region\n"
           "\tfluxweighted_region_cartesian: Flux-weighted in a prescribed\n"
           "\t                               cartesian region\n"
           "\tprescribed_positions: Prescribed particle positions\n"
           "\tprescribed_positions_masses: Prescribed particle positions and\n"
           "                               masses\n"
           "\tprescribed_positions_masses_tags: Prescribed particle\n"
           "                                    positions, masses, and tags\n"
           "\tuniform_inlet_continuous: Continuous injection homogeneous at\n"
           "\t                          the inlet\n"
           "\tfluxweighted_inlet_continuous: Continuous injection\n"
           "\t                               flux-weighted at the inlet\n"
           "- Initial distance from solid phase (pass only for type\n"
           "  uniform_near_solid)\n"
           "- Region boundaries (pass only for types uniform_region_cartesian\n"
           "  or fluxweighted_region_cartesian)\n"
           "- Filename (with full path) for prescribed positions, positions\n"
           "  and masses, or positions, masses, and tags (pass only for type\n"
           "  prescribed_positions, prescribed_position_masses, or\n"
           "  prescribed_positions_masses_tags, respectively)\n"
           "- Total injected mass in each injection step (pass for all types\n"
           "  except prescribed_position_masses and\n"
           "  prescribed_positions_masses_tags)\n"
           "- Initial injection time\n"
           "- Final injection time (pass only for types\n"
           "  uniform_inlet_continuous or fluxweighted_inlet_continuous)\n"
           "- Maximum timestep for continuous injection discretization in\n"
           "  units of advection time (pass only for types\n"
           "  uniform_inlet_continuous or fluxweighted_inlet_continuous)\n"
           "- Maximum timestep for continuous injection discretization in\n"
           "  units of diffusion time (pass only for types\n"
           "  uniform_inlet_continuous or fluxweighted_inlet_continuous)\n"
           "- Maximum timestep for continuous injection discretization in\n"
           "  units of reaction times (pass only for types\n"
           "  uniform_inlet_continuous or fluxweighted_inlet_continuous)\n"
           "--------------------------------------------------------------\n";
  }
};

/** \class InitialCondition_Cases PTOF/InitialCondition.h
 * "PTOF/InitialCondition.h" \brief InitialCondition object to handle
 * implemented initial condition types. */
template <typename Geometry, typename VelocityField, typename ParticleMaker,
          typename Mask = useful::Empty>
class InitialCondition_Cases {
public:
  using Parameters =
      InitialConditionParameters_Cases; /**< Initial condition parameters.*/
  const Parameters parameters;          /**< Initial condition parameters. */

private:
  typename Geometry::Locator const &_locator; /**< Object to search mesh. */
  VelocityField const &_velocity_field; /**< Underlying flow field.       */
  ParticleMaker _particle_maker; /**< Makes a particle given a position.*/
  std::size_t _nr_particles;     /**< Number of particles per injection step. */
  std::mt19937 _rng{std::random_device{}()}; /**< Random number generator. */
  Mask _mask{};          /**< Mask to disallow cells/faces        */
  double _threshold{0.}; /**< Threshold for mask        */

  struct PositionMaker {
    double time = 0.;
    auto operator()(Foam::point const &position) { return position; }
  } _position_maker; /** Functor to forward positions and hold time value.*/

public:
  /** Constructor.
   \param geometry Domain geometry info and utilities.
   \param velocity_field Velocity field as a function of state.
   \param particle_maker Functor to make a particle given a position.
   \param nr_particles Number of particles per injection step.
   \param parameters Output parameters. */
  InitialCondition_Cases(Geometry const &geometry,
                         VelocityField const &velocity_field,
                         ParticleMaker particle_maker, std::size_t nr_particles,
                         Parameters const &parameters)
      : parameters{parameters}, _locator{geometry.locator},
        _velocity_field{velocity_field}, _particle_maker{particle_maker},
        _nr_particles{nr_particles} {}

  /** Constructor.
  \param geometry Domain geometry info and utilities.
  \param velocity_field Velocity field as a function of state.
  \param particle_maker Functor to make a particle given a position.
  \param nr_particles Number of particles per injection step.
  \param parameters Output parameters.
  \param mask Scalar field assigning values to mesh cell indices through
  operator[]. \param threshold Threshold for the mask, such that cells where the
  mask is above the threshold are considered. */
  InitialCondition_Cases(Geometry const &geometry,
                         VelocityField const &velocity_field,
                         ParticleMaker particle_maker, std::size_t nr_particles,
                         Parameters const &parameters, Mask &&mask,
                         double threshold = 0.)
      : parameters{parameters}, _locator{geometry.locator},
        _velocity_field{velocity_field}, _particle_maker{particle_maker},
        _nr_particles{nr_particles}, _mask{std::forward<Mask>(mask)},
        _threshold{threshold} {}

  /** \brief Make particles according to prescribed initial condition and
   particle maker. \param nr_particles Number of particles to make \return
   Container with particles.  */
  auto operator()(std::size_t nr_particles) {
    return make_particles(nr_particles, parameters.type, _particle_maker);
  }

  /** \brief Make particles according to prescribed initial condition and
  particle maker, with \c nr_particles particles as specified in \c parameters.
  \return Container with particles.  */
  auto operator()() { return this->operator()(_nr_particles); }

  /** \brief Make a single particle according to prescribed initial condition
  and particle maker. \return One particle.  */
  auto make_particle() { return this->operator()(1)[0]; }

  /** \brief Make positions according to prescribed initial condition.
  \param nr_positions Number of positions to make
  \return Container with position.  */
  auto make_positions(std::size_t nr_positions) {
    // This works by slightly abusing make_particles with a position maker
    // instead of a particle maker, in order to produce only positions instead
    // of particles, but according to the same prescription as for producing
    // particles.
    return make_particles(nr_positions, parameters.type, _position_maker);
  }

  /** \brief Make a single position according to prescribed initial condition
  and particle maker. \return One position.  */
  auto make_position() { return make_positions(1)[0]; }

private:
  /** \brief Make given number of particles.
   \param nr_particles Number of particles to make.
   \param type Initial condition type.
   \param particle_maker Functor to make a particle given a position.
   \return Container with particles.
   */
  template <typename ParticleMakerOther>
  auto make_particles(std::size_t nr_particles, InitialConditions::Type type,
                      ParticleMakerOther &particle_maker) {
    switch (type) {
    case InitialConditions::Type::uniform:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return uniform(nr_particles, _locator.mesh(), _rng, particle_maker);
      else
        return uniform(nr_particles, _locator.mesh(), _rng, particle_maker,
                       _mask, _threshold);
    case InitialConditions::Type::fluxweighted:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return fluxweighted(nr_particles, _velocity_field, _locator.mesh(),
                            _rng, particle_maker);
      else
        return fluxweighted(nr_particles, _velocity_field, _locator.mesh(),
                            _rng, particle_maker, _mask, _threshold);
    case InitialConditions::Type::uniform_inlet:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return uniform_patches(nr_particles, {"inlet"}, _locator, _rng,
                               particle_maker);
      else
        return uniform_patches(nr_particles, {"inlet"}, _locator, _rng,
                               particle_maker, _mask, _threshold);
    case InitialConditions::Type::fluxweighted_inlet:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return fluxweighted_patches(nr_particles, {"inlet"}, _velocity_field,
                                    _locator, _rng, particle_maker);
      else
        return fluxweighted_patches(nr_particles, {"inlet"}, _velocity_field,
                                    _locator, _rng, particle_maker, _mask,
                                    _threshold);
    case InitialConditions::Type::uniform_solid:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return uniform_patches(nr_particles, {"wallFluidSolid"}, _locator, _rng,
                               particle_maker);
      else
        return uniform_patches(nr_particles, {"wallFluidSolid"}, _locator, _rng,
                               particle_maker, _mask, _threshold);
    case InitialConditions::Type::uniform_near_solid:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return uniform_near_boundary_patches(nr_particles, {"wallFluidSolid"},
                                             parameters.distance_wall, _locator,
                                             _rng, particle_maker);
      else
        return uniform_near_boundary_patches(
            nr_particles, {"wallFluidSolid"}, parameters.distance_wall,
            _locator, _rng, particle_maker, _mask, _threshold);
    case InitialConditions::Type::uniform_region_cartesian:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return uniform_cells(
            nr_particles,
            cell_ids_region_cartesian(parameters.region_boundaries, _locator),
            _locator.mesh(), _rng, particle_maker);
      else
        return uniform_cells(
            nr_particles,
            apply_mask(cell_ids_region_cartesian(parameters.region_boundaries,
                                                 _locator),
                       _mask, _threshold),
            _locator.mesh(), _rng, particle_maker);
    case InitialConditions::Type::fluxweighted_region_cartesian:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return fluxweighted_cells(
            nr_particles,
            cell_ids_region_cartesian(parameters.region_boundaries, _locator),
            _velocity_field, _locator.mesh(), _rng, particle_maker);
      else
        return fluxweighted_cells(
            nr_particles,
            apply_mask(cell_ids_region_cartesian(parameters.region_boundaries,
                                                 _locator),
                       _mask, _threshold),
            _velocity_field, _locator.mesh(), _rng, particle_maker);
    case InitialConditions::Type::prescribed_positions:
      return prescribed_positions(nr_particles, parameters.filename_data,
                                  particle_maker);
    case InitialConditions::Type::prescribed_positions_masses:
      return prescribed_positions_masses(nr_particles, parameters.filename_data,
                                         particle_maker);
    case InitialConditions::Type::prescribed_positions_masses_tags:
      return prescribed_positions_masses_tags(
          nr_particles, parameters.filename_data, particle_maker);
    case InitialConditions::Type::uniform_inlet_continuous:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return continuous_injection(
            [this](std::size_t nr_particles,
                   ParticleMakerOther &particle_maker) {
              return uniform_patches(nr_particles, {"inlet"}, _locator, _rng,
                                     particle_maker);
            },
            nr_particles, parameters.time_min, parameters.time_max,
            parameters.time_step, particle_maker);
      else
        return continuous_injection(
            [this](std::size_t nr_particles,
                   ParticleMakerOther &particle_maker) {
              return uniform_patches(nr_particles, {"inlet"}, _locator, _rng,
                                     particle_maker, _mask, _threshold);
            },
            nr_particles, parameters.time_min, parameters.time_max,
            parameters.time_step, particle_maker);
    case InitialConditions::Type::fluxweighted_inlet_continuous:
      if constexpr (std::is_same_v<Mask, useful::Empty>)
        return continuous_injection(
            [this](std::size_t nr_particles,
                   ParticleMakerOther &particle_maker) {
              return fluxweighted_patches(nr_particles, {"inlet"},
                                          _velocity_field, _locator, _rng,
                                          particle_maker);
            },
            nr_particles, parameters.time_min, parameters.time_max,
            parameters.time_step, particle_maker);
      else
        return continuous_injection(
            [this](std::size_t nr_particles,
                   ParticleMakerOther &particle_maker) {
              return fluxweighted_patches(nr_particles, {"inlet"},
                                          _velocity_field, _locator, _rng,
                                          particle_maker, _mask, _threshold);
            },
            nr_particles, parameters.time_min, parameters.time_max,
            parameters.time_step, particle_maker);
    default:
      throw std::runtime_error{std::string{"Initial condition type "} +
                               InitialConditions::name(parameters.type) +
                               " not supported"};
    }
  }

public:
  /** \brief Output information about current object. */
  template <typename OStream> void info_runtime(OStream &output) const {
    output
        << "--------------------------------------------------------------\n"
           "Initial condition\n"
           "--------------------------------------------------------------\n"
           "Type: " +
               InitialConditions::name(parameters.type) + "\n"
        << "--------------------------------------------------------------\n";
  }
};
template <typename Geometry, typename VelocityField, typename ParticleMaker,
          typename Parameters, typename Mask>
InitialCondition_Cases(Geometry const &, VelocityField const &, ParticleMaker,
                       std::size_t, Parameters const &, Mask &&, double)
    -> InitialCondition_Cases<Geometry, VelocityField, ParticleMaker, Mask>;
} // namespace ptof

#endif /* PTOF_INITIALCONDITIONS_H */
