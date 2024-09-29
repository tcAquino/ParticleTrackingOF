/**
 \file PTOF/InitialCondition.h
 \author Tomás Aquino
 \date 29/09/2024
*/

#ifndef PTOF_INITIALCONDITION_H
#define PTOF_INITIALCONDITION_H

#include "CTRW/Meta.h"
#include "General/Useful.h"
#include <cstddef>
#include <fieldTypes.H>
#include <map>
#include <point.H>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace ptof {
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
 \param rng Random number param.
 \generator particle_maker Functor to make a particle given a position.
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

/** \brief Make particles distributed randomly uniformly over given mesh cell
faces. \details Place particles at face centers. \param nr_particles Number of
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
 mask is above threshold.
 \details Place particles at cell centers.
 \param nr_particles Number of particles to make.
 \param mesh Mesh object.
 \param rng Random number generator.
 \param particle_maker Functor to make a particle given a position.
 \param mask Scalar field assigning values to mesh cell indices through
 operator[].
 \param threshold Threshold for the mask, such that cells where the
 mask is above the threshold are considered.
 \return Container with particles.
 */
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
operator[].
\param threshold Threshold for the mask, such that cells where the
mask is above the threshold are considered.
\return Container with particles. */
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
\c nr_particles positions are used.
\param particle_maker Functor to make a particle given a position.
\return Container with particles.
\note
- Empty lines are ignored.
- Unspecified position components are set to zero. */
template <typename ParticleMaker>
auto prescribed_positions(std::size_t nr_particles,
                          std::string const &filename_data,
                          ParticleMaker &particle_maker) {
  auto input = useful::open_read(filename_data);
  std::vector<std::vector<double>> position_data;
  position_data.reserve(nr_particles);
  for (std::size_t pp = 0; pp < nr_particles; ++pp) {
    auto split_line = useful::split_line(input);
    if (split_line.empty())
      continue;
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
\return Container with particles.
\note
- Empty lines are ignored.
- Unspecified position components are set to zero. */
template <typename ParticleMaker>
auto prescribed_positions(std::string const &filename_data,
                          ParticleMaker &particle_maker) {
  auto input = useful::open_read(filename_data);
  std::vector<std::vector<double>> position_data;
  std::vector<std::string> split_line;
  while (useful::split_line(input, split_line)) {
    if (split_line.empty())
      continue;
    position_data.emplace_back();
    position_data.back().reserve(split_line.size());
    for (auto const &position_component : split_line)
      position_data.back().push_back(std::stod(position_component));
    split_line.clear();
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
first \c nr_particles lines are used.
\param particle_maker Functor to make a particle given a position.
\return Container with particles.
\note
- Empty lines are ignored.
- Unspecified position components are set to zero. */
template <typename ParticleMaker>
auto prescribed_positions_masses(std::size_t nr_particles,
                                 std::string const &filename_position_mass_data,
                                 ParticleMaker &particle_maker) {
  auto input = useful::open_read(filename_position_mass_data);
  std::vector<std::vector<double>> position_data;
  position_data.reserve(nr_particles);
  std::vector<double> mass_data;
  mass_data.reserve(nr_particles);
  for (std::size_t pp = 0; pp < nr_particles; ++pp) {
    auto split_line = useful::split_line(input);
    if (split_line.empty())
      continue;
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
\return Container with particles.
\note
- Empty lines are ignored.
- Unspecified position components are set to zero. */
template <typename ParticleMaker>
auto prescribed_positions_masses(std::string const &filename_data,
                                 ParticleMaker &particle_maker) {
  auto input = useful::open_read(filename_data);
  std::vector<std::vector<double>> position_data;
  std::vector<double> mass_data;
  std::vector<std::string> split_line;
  while (useful::split_line(input, split_line)) {
    if (split_line.empty())
      continue;
    position_data.emplace_back();
    position_data.back().reserve(split_line.size());
    for (std::size_t ii = 0; ii < split_line.size() - 1; ++ii)
      position_data.back().push_back(std::stod(split_line[ii]));
    mass_data.push_back(std::stod(split_line.back()));
    split_line.clear();
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
      throw std::runtime_error{"Failed to prescribe mass because "
                               "ParticleMaker does not define mass"};
    if constexpr (meta::has_tag_v<ParticleMaker>)
      particle_maker.tag = tag_data[pp];
    else
      throw std::runtime_error{"Failed to prescribe tag because "
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
  auto input = useful::open_read(filename_position_mass_data);
  std::vector<std::vector<double>> position_data;
  position_data.reserve(nr_particles);
  std::vector<double> mass_data;
  mass_data.reserve(nr_particles);
  std::vector<std::size_t> tag_data;
  tag_data.reserve(nr_particles);
  for (std::size_t pp = 0; pp < nr_particles; ++pp) {
    auto split_line = useful::split_line(input);
    if (split_line.empty())
      continue;
    if (split_line.size() < 2)
      throw std::runtime_error{
          std::string{"Could not parse prescribed tag for "} +
          useful::ordinal(pp + 1) + " particle"};
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
line.
\param particle_maker Functor to make a particle given a position.
\return Container with particles.
\note
- Empty lines are ignored.
- Unspecified position components are set to zero. */
template <typename ParticleMaker>
auto prescribed_positions_masses_tags(std::string const &filename_data,
                                      ParticleMaker &particle_maker) {
  auto input = useful::open_read(filename_data);
  std::vector<std::vector<double>> position_data;
  std::vector<double> mass_data;
  std::vector<std::size_t> tag_data;
  std::vector<std::string> split_line;
  while (useful::split_line(input, split_line)) {
    if (split_line.empty())
      continue;
    if (split_line.size() < 2)
      throw std::runtime_error{
          std::string{"Could not parse prescribed tag for "} +
          useful::ordinal(position_data.size() + 1) + " particle"};
    position_data.emplace_back();
    position_data.back().reserve(split_line.size());
    for (std::size_t ii = 0; ii < split_line.size() - 2; ++ii)
      position_data.back().push_back(std::stod(split_line[ii]));
    mass_data.push_back(std::stod(split_line[split_line.size() - 2]));
    tag_data.push_back(std::stod(split_line.back()));
    split_line.clear();
  }
  input.close();

  return prescribed_positions_masses(position_data, mass_data, particle_maker);
}
} // namespace ptof

#endif /* PTOF_INITIALCONDITION_H */
