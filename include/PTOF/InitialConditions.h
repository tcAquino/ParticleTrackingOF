/**
 \file PTOF/InitialConditions.h
 \author Tomás Aquino
 \date 24/02/2022
*/

#ifndef PTOF_INITIALCONDITIONS_H
#define PTOF_INITIALCONDITIONS_H

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <exception>
#include <numeric>
#include <random>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <fvMesh.H>
#include "PTOF/Boundary.h"
#include "PTOF/Useful.h"

namespace ptof
{
  /** \struct InitialConditions PTOF/InitialConditions.h "PTOF/InitialConditions.h"
   * \brief Keep track of names and
   * types of initial conditions. */
  struct InitialConditions
  {
   /** \enum Type
     *  \brief Implemented types. */
    enum class Type
    {
      uniform,                        /**< Homogeneous throughout the domain.                     */
      fluxweighted,                   /**< Flux-weighted throughout the domain.                   */
      uniform_inlet,                  /**< Homogeneous at the inlet.                              */
      fluxweighted_inlet,             /**< Flux-weighted at the inlet.                            */
      uniform_solid,                  /**< Homogeneous at the solid surface.                      */
      uniform_near_solid,             /**< Homogeneous at a fixed distance to the solid interface.*/
      uniform_region_cartesian,       /**< Homogeneous in a Cartesian region        */
      fluxweighted_region_cartesian,  /**< Flux-weighted in a Cartesian region        */
      prescribed_positions,           /**< Prescribed positions                               */
      uniform_inlet_continuous,       /**< Continuous injection homogeneous at the inlet.         */
      fluxweighted_inlet_continuous,  /**< Continuous injection flux-weighted at the inlet.       */
    };
    
    /** Type from name. */
    static auto type(std::string const& name)
    { return string_to_type.at(name); }
    
    /** Name from name. */
    static auto name(Type type)
    { return type_to_string.at(type); }
    
    /** Check if name exists. */
    static auto contains(std::string const& name)
    { return string_to_type.count(name); }
    
    /** Map names to types. */
    inline static const
    std::unordered_map<std::string, Type> string_to_type
    {
      { "uniform", Type::uniform },
      { "fluxweighted", Type::fluxweighted },
      { "uniform_inlet", Type::uniform_inlet },
      { "fluxweighted_inlet", Type::fluxweighted_inlet },
      { "uniform_solid", Type::uniform_solid },
      { "uniform_near_solid", Type::uniform_near_solid },
      { "uniform_region_cartesian", Type::uniform_region_cartesian },
      { "fluxweighted_region_cartesian", Type::fluxweighted_region_cartesian },
      { "prescribed_positions", Type::prescribed_positions },
      { "uniform_inlet_continuous", Type::uniform_inlet_continuous },
      { "fluxweighted_inlet_continuous", Type::fluxweighted_inlet_continuous }
    };
    
    /** Map types to names. */
    inline static const
    std::unordered_map<Type, std::string> type_to_string
    {
      { Type::uniform, "uniform", },
      { Type::fluxweighted, "fluxweighted" },
      { Type::uniform_inlet, "uniform_inlet" },
      { Type::fluxweighted_inlet, "fluxweighted_inlet" },
      { Type::uniform_solid, "uniform_solid" },
      { Type::uniform_near_solid, "uniform_near_solid" },
      { Type::uniform_region_cartesian, "uniform_region_cartesian" },
      { Type::fluxweighted_region_cartesian, "fluxweighted_region_cartesian" },
      { Type::fluxweighted_region_cartesian, "prescribed_positions" },
      { Type::uniform_inlet_continuous, "uniform_inlet_continuous" },
      { Type::fluxweighted_inlet_continuous, "fluxweighted_inlet_continuous" }
    };
  };
  
  /** \brief Remove cell/face indices if mask is above threshold.
   \param cell_or_face_ids Container with mesh cell or face indices.
   \param mask Scalar field assigning values to mesh cell indices through operator[].
   \param threshold Threshold for the mask, such that cells where the mask is above the threshold are considered. */
  template
  <typename Container, typename Mask>
  void apply_mask_in_place
  (Container& cell_or_face_ids,
   Mask const& mask,
   double threshold = 0.)
  {
    useful::swap_erase_if(cell_or_face_ids,
                          [&mask, threshold](auto ii){ return mask[ii] > threshold; });
  }
  
  /** \brief Remove cell/face indices if mask is above threshold.
  \param cell_or_face_ids Container with mesh cell or face indices.
  \param mask Scalar field assigning values to mesh cell indices through operator[].
  \param threshold Threshold for the mask, such that cells where the mask is above the threshold are considered.
  \return \p cell_or_face_ids with indices of elements where mask is above threshold removed. */
  template
  <typename Container, typename Mask>
  auto apply_mask
  (Container const& cell_or_face_ids,
   Mask const& mask,
   double threshold = 0.)
  {
    auto masked_ids = cell_or_face_ids;
    apply_mask_in_place(masked_ids, mask, threshold);

    return masked_ids;
  }
  
  /** \brief Make particles distributed randomly uniformly over given cells.
   \details Place particles at cell centers.
   \param nr_particles Number of particles to make.
   \param cell_ids Mesh cell indices.
   \param mesh OpenFOAM mesh.
   \param rng Random number generator.
   \param particle_maker Functor to make a particle given a position.
   \return Container with particles. */
  template
  <typename Container,
  typename Mesh,
  typename RNG,
  typename ParticleMaker>
  auto uniform_cells
  (std::size_t nr_particles,
   Container const& cell_ids,
   Mesh const& mesh,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    if (cell_ids.size() == 0 && nr_particles != 0)
      throw std::runtime_error{
        "No available cells to place particles" };
    
    // Weigh probabilities with cell volumes.
    std::vector<double> weights;
    weights.reserve(cell_ids.size());
    for (auto cell : cell_ids)
      weights.push_back(cell_volume(cell, mesh));
    std::discrete_distribution<Foam::label> dist{
      weights.begin(), weights.end() };
    
    using Particle = decltype(particle_maker(Foam::point{}));
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    for (std::size_t pp = 0; pp < nr_particles; ++pp)
      particles.push_back(particle_maker(cell_center(cell_ids[dist(rng)],
                                                     mesh)));
    return particles;
  }
  
  /** \brief Make particles distributed over given cells with probability proportional to local velocity magnitude.
  \details Place particles at cell centers.
  \param nr_particles Number of particles to make.
  \param cell_ids Mesh cell indices.
  \param velocity_field Velocity field as a function of state.
  \param mesh OpenFOAM mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template
  <typename Container,
  typename VelocityField,
  typename Mesh,
  typename RNG,
  typename ParticleMaker>
  auto fluxweighted_cells
  (std::size_t nr_particles,
   Container const& cell_ids,
   VelocityField const& velocity_field,
   Mesh const& mesh,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    // Weigh probabilities with cell volumes
    // and velocity magnitudes at center.
    if (cell_ids.size() == 0 && nr_particles != 0)
      throw std::runtime_error{
        "No available cells to place particles" };
    
    std::vector<double> weights;
    weights.reserve(cell_ids.size());
    for (auto cell : cell_ids)
      weights.push_back(cell_volume(cell, mesh)
                        *Foam::mag(velocity_field(cell_center(cell, mesh),
                                                  cell)));
    std::discrete_distribution<Foam::label> dist{
      weights.begin(), weights.end() };
    
    using Particle = decltype(particle_maker(Foam::point{}));
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    for (std::size_t pp = 0; pp < nr_particles; ++pp)
      particles.push_back(particle_maker(cell_center(cell_ids[dist(rng)], mesh)));
    return particles;
  }
  
  /** \brief Make particle at given mesh face's center.
  \param face Mesh face index.
  \param mesh OpenFOAM mesh.
  \param locator Object to locate positions in mesh.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template
  <
  typename Mesh,
  typename Locator,
  typename ParticleMaker>
  auto make_particle_at_face_center
  (Foam::label face,
   Mesh const& mesh,
   Locator const& locator,
   ParticleMaker& particle_maker)
  {
    if (face_center_is_in_cell(face, mesh, locator))
      return particle_maker(face_center(face, mesh));
    else
    {
      auto center = face_center(face, mesh);
      auto owner_cell = mesh.owner()[face];
      std::cerr << "Warning: Face center of cell "
                << face << " at "
                << "("
                << center[0] << ", "
                << center[1] << ", "
                << center[2] << ")"
                << " is not within owner cell "
                << owner_cell << ". "
                << "Trying small displacement toward cell center... ";
      auto offset_face_center = offset_inward_face(face, locator);
      if (locator.mesh_search().findCell(offset_face_center) == owner_cell)
      {
        std::cerr << "Success\n";
        return particle_maker(offset_face_center);
      }
      else
      {
        std::cerr << "Failed. "
                  << "Placing particle at center of owner face\n";
        return particle_maker(cell_center(owner_cell, mesh));
      }
    }
  }

  /** Make particles randomly uniformly
   *  over given mesh cell faces.
   *  Place particles at center of faces. */
  /** \brief Make particles distributed randomly uniformly over given mesh cell faces.
  \details Place particles at cell centers.
  \param nr_particles Number of particles to make.
  \param face_ids Mesh face indices.
  \param mesh OpenFOAM mesh.
  \param locator Object to locate positions in mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template
  <typename Container,
  typename Mesh,
  typename Locator,
  typename RNG,
  typename ParticleMaker>
  auto uniform_faces
  (std::size_t nr_particles,
   Container const& face_ids,
   Mesh const& mesh,
   Locator const& locator,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    if (face_ids.size() == 0 && nr_particles != 0)
      throw std::runtime_error{
        "No available faces to place particles" };
    
    std::vector<double> weights;
    weights.reserve(face_ids.size());
    for (auto face : face_ids)
      weights.push_back(face_area(face, mesh));
    std::discrete_distribution<Foam::label> dist{
      weights.begin(), weights.end() };
    
    using Particle = decltype(particle_maker(Foam::point{}));
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    for (std::size_t pp = 0; pp < nr_particles; ++pp)
      particles.push_back(make_particle_at_face_center(face_ids[dist(rng)],
                                                       mesh, locator,
                                                       particle_maker));
    return particles;
  }
  
  /** \brief Make particles distributed over given cell faces with probability proportional to local velocity magnitude.
  \details Place particles at face centers.
  \param nr_particles Number of particles to make.
  \param face_ids Mesh face indices.
  \param velocity_field Velocity field as a function of state.
  \param mesh OpenFOAM mesh.
  \param locator Object to locate positions in mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template
  <typename Container,
  typename VelocityField,
  typename Mesh,
  typename Locator,
  typename RNG,
  typename ParticleMaker>
  auto fluxweighted_faces
  (std::size_t nr_particles,
   Container const& face_ids,
   VelocityField const& velocity_field,
   Mesh const& mesh,
   Locator const& locator,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    if (face_ids.size() == 0 && nr_particles != 0)
      throw std::runtime_error{
        "No available faces to place particles" };
    
    std::vector<double> weights;
    weights.reserve(face_ids.size());
    for (auto face : face_ids)
      weights.push_back(face_area(face, mesh)*
                        Foam::mag(velocity_field(face_center(face, mesh),
                                                  mesh.owner()[face])));
    std::discrete_distribution<Foam::label> dist{
      weights.begin(), weights.end() };
    
    using Particle = decltype(particle_maker(Foam::point{}));
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    for (std::size_t pp = 0; pp < nr_particles; ++pp)
      particles.push_back(make_particle_at_face_center(face_ids[dist(rng)],
                                                       mesh, locator,
                                                       particle_maker));
    return particles;
  }
  
  /** Make particles randomly uniformly
   *  at given distance to given mesh cell face centers. */
  /** \brief Make particles distributed randomly uniformly near given cell face centers.
  \param nr_particles Number of particles to make.
  \param face_ids Mesh face indices.
  \param distance Distance from face centers.
  \param mesh OpenFOAM mesh.
  \param locator Object to locate positions in mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template
  <typename Container,
  typename Mesh,
  typename Locator,
  typename RNG,
  typename ParticleMaker>
  auto uniform_near_boundary_faces
  (std::size_t nr_particles,
   Container const& face_ids,
   double distance,
   Mesh const& mesh,
   Locator const& locator,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    if (face_ids.size() == 0 && nr_particles != 0)
      throw std::runtime_error{
        "No available faces to place particles" };
    
    std::vector<double> weights;
    weights.reserve(face_ids.size());
    for (auto face : face_ids)
      weights.push_back(face_area(face, mesh));
    std::discrete_distribution<Foam::label> dist{
      weights.begin(), weights.end() };
    
    using Particle = decltype(particle_maker(Foam::point{}));
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    for (std::size_t pp = 0; pp < nr_particles; ++pp)
    {
      auto rand = dist(rng);
      auto unit_normal = unit_normal_inward(face_ids[rand], mesh);
      auto position = face_center(face_ids[rand], mesh);
        + distance*unit_normal;
      if (outside(locator.mesh_search().findCell(position))
          || locator.mesh_search().findNearestBoundaryFace(position)
            != face_ids[rand])
        continue;
      particles.push_back(particle_maker(position));
    }
    return particles;
  }
  
  /** \brief Make particles distributed randomly uniformly over all mesh cells.
  \details Place particles at cell centers.
  \param nr_particles Number of particles to make.
  \param mesh OpenFOAM mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template
  <typename Mesh,
  typename RNG,
  typename ParticleMaker>
  auto uniform
  (std::size_t nr_particles,
   Mesh const& mesh,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    return uniform_cells(nr_particles,
                         all_cell_ids(mesh),
                         mesh, rng, particle_maker);
  }
  
  /** \brief Make particles distributed randomly uniformly over mesh cells where mask is above threshold.
   \details Place particles at cell centers.
   \param nr_particles Number of particles to make.
   \param mesh OpenFOAM mesh.
   \param rng Random number generator.
   \param particle_maker Functor to make a particle given a position.
   \param mask Scalar field assigning values to mesh cell indices through operator[].
   \param threshold Threshold for the mask, such that cells where the mask is above the threshold are considered.
   \return Container with particles. */
  template
  <typename Mesh,
  typename RNG,
  typename ParticleMaker,
  typename Mask>
  auto uniform
  (std::size_t nr_particles,
   Mesh const& mesh,
   RNG& rng,
   ParticleMaker& particle_maker,
   Mask const& mask,
   double threshold = 0.)
  {
    return uniform_cells(nr_particles,
                         apply_mask(all_cell_ids(mesh),
                                    mask, threshold),
                         mesh, rng, particle_maker);
  }
  
  /** \brief Make particles distributed with probability proportional to local velocity magnitude over all mesh cells.
  \details Place particles at cell centers.
  \param nr_particles Number of particles to make.
  \param velocity_field Velocity field as a function of state.
  \param mesh OpenFOAM mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template
  <typename VelocityField,
  typename Mesh,
  typename RNG,
  typename ParticleMaker>
  auto fluxweighted
  (std::size_t nr_particles,
   VelocityField const& velocity_field,
   Mesh const& mesh,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    return fluxweighted_cells(nr_particles,
                               all_cell_ids(mesh),
                               velocity_field,
                               mesh, rng, particle_maker);
  }
  
  /** \brief Make particles distributed with probability proportional to local velocity magnitude over mesh cells where mask is above threshold.
  \details Place particles at cell centers.
  \param nr_particles Number of particles to make.
  \param velocity_field Velocity field as a function of state.
  \param mesh OpenFOAM mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \param mask Scalar field assigning values to mesh cell indices through operator[].
  \param threshold Threshold for the mask, such that cells where the mask is above the threshold are considered.
  \return Container with particles. */
  template
  <typename VelocityField,
  typename Mesh,
  typename RNG,
  typename ParticleMaker,
  typename Mask>
  auto fluxweighted
  (std::size_t nr_particles,
   VelocityField const& velocity_field,
   Mesh const& mesh,
   RNG& rng,
   ParticleMaker& particle_maker,
   Mask const& mask,
   double threshold = 0.)
  {
    return fluxweighted_cells(nr_particles,
                               apply_mask(all_cell_ids(mesh),
                                          mask, threshold),
                               velocity_field,
                               mesh, rng, particle_maker);
  }
  
  /** \brief Make particles distributed randomly uniformly over cell faces in given patches.
  \details Place particles at face centers.
  \param nr_particles Number of particles to make.
  \param patch_names Names of patches.
  \param mesh OpenFOAM mesh.
  \param locator Object to locate positions in mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template
  <typename Mesh,
  typename Locator,
  typename RNG,
  typename ParticleMaker>
  auto uniform_patches
  (std::size_t nr_particles,
   std::vector<std::string> const& patch_names,
   Mesh const& mesh,
   Locator const& locator,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    return uniform_faces(nr_particles,
                         patches_face_ids(mesh, patch_names),
                         mesh, locator,
                         rng, particle_maker);
  }
  
  /** \brief Make particles distributed randomly uniformly over cell faces in given patches where mask is above threshold.
  \details Place particles at face centers.
  \param nr_particles Number of particles to make.
  \param patch_names Names of patches.
  \param mesh OpenFOAM mesh.
  \param locator Object to locate positions in mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \param mask Scalar field assigning values to mesh cell indices through operator[].
  \param threshold Threshold for the mask, such that cells where the mask is above the threshold are considered.
  \return Container with particles. */
  template
  <typename Mesh,
  typename Locator,
  typename RNG,
  typename ParticleMaker,
  typename Mask>
  auto uniform_patches
  (std::size_t nr_particles,
   std::vector<std::string> const& patch_names,
   Mesh const& mesh,
   Locator const& locator,
   RNG& rng,
   ParticleMaker& particle_maker,
   Mask const& mask,
   double threshold = 0.)
  {
    return uniform_faces(nr_particles,
                         apply_mask(patches_face_ids(mesh, patch_names),
                                    mask, threshold),
                         mesh, locator,
                         rng, particle_maker);
  }
  
  /** \brief Make particles distributed with probability proportional to local velocity magnitude over cell faces in given patches.
  \details Place particles at face centers.
  \param nr_particles Number of particles to make.
  \param patch_names Names of patches.
  \param velocity_field Velocity field as a function of state.
  \param mesh OpenFOAM mesh.
  \param locator Object to locate positions in mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template
  <typename VelocityField,
  typename Mesh,
  typename Locator,
  typename RNG,
  typename ParticleMaker>
  auto fluxweighted_patches
  (std::size_t nr_particles,
   std::vector<std::string> const& patch_names,
   VelocityField const& velocity_field,
   Mesh const& mesh,
   Locator const& locator,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    return fluxweighted_faces(nr_particles,
                               patches_face_ids(mesh, patch_names),
                               velocity_field,
                               mesh, locator,
                               rng, particle_maker);
  }
  
  /** \brief Make particles distributed with probability proportional to local velocity magnitude over cell faces in given patches where mask is above threshold.
  \details Place particles at face centers.
  \param nr_particles Number of particles to make.
  \param patch_names Names of patches.
  \param velocity_field Velocity field as a function of state.
  \param mesh OpenFOAM mesh.
  \param locator Object to locate positions in mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \param mask Scalar field assigning values to mesh cell indices through operator[].
  \param threshold Threshold for the mask, such that cells where the mask is above the threshold are considered.
  \return Container with particles. */
  template
  <typename VelocityField,
  typename Mesh,
  typename Locator,
  typename RNG,
  typename ParticleMaker,
  typename Mask>
  auto fluxweighted_patches
  (std::size_t nr_particles,
   std::vector<std::string> const& patch_names,
   VelocityField const& velocity_field,
   Mesh const& mesh,
   Locator const& locator,
   RNG& rng,
   ParticleMaker& particle_maker,
   Mask const& mask,
   double threshold = 0.)
  {
    return fluxweighted_faces(nr_particles,
                               apply_mask(patches_face_ids(mesh, patch_names),
                                          mask, threshold),
                               velocity_field,
                               mesh, locator,
                               rng, particle_maker);
  }
  
  /** \brief Make particles distributed randomly uniformly near cell faces in given patches.
  \param nr_particles Number of particles to make.
  \param patch_names Names of patches.
  \param distance Distance from face centers.
  \param mesh OpenFOAM mesh.
  \param locator Object to locate positions in mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template
  <typename Mesh,
  typename Locator,
  typename RNG,
  typename ParticleMaker>
  auto uniform_near_boundary_patches
  (std::size_t nr_particles,
   std::vector<std::string> const& patch_names,
   double distance,
   Mesh const& mesh,
   Locator const& locator,
   RNG& rng,
   ParticleMaker particle_maker)
  {
    return uniform_near_boundary_faces(nr_particles,
                                       patches_face_ids(mesh, patch_names),
                                       distance,
                                       mesh, locator,
                                       rng, particle_maker);
  }
  
  /** \brief Make particles distributed randomly uniformly near cell faces in given patches where mask is above threshold.
  \param nr_particles Number of particles to make.
  \param patch_names Names of patches.
  \param distance Distance from face centers.
  \param mesh OpenFOAM mesh.
  \param locator Object to locate positions in mesh.
  \param rng Random number generator.
  \param particle_maker Functor to make a particle given a position.
  \param mask Scalar field assigning values to mesh cell indices through operator[].
  \param threshold Threshold for the mask, such that cells where the mask is above the threshold are considered.
  \return Container with particles. */
  template
  <typename Mesh,
  typename Locator,
  typename RNG,
  typename ParticleMaker,
  typename Mask>
  auto uniform_near_boundary_patches
  (std::size_t nr_particles,
   std::vector<std::string> const& patch_names,
   double distance,
   Mesh const& mesh,
   Locator const& locator,
   RNG& rng,
   ParticleMaker particle_maker,
   Mask const& mask,
   double threshold = 0.)
  {
    return uniform_near_boundary_faces(nr_particles,
                                       apply_mask(patches_face_ids(mesh, patch_names),
                                                  mask, threshold),
                                       distance,
                                       mesh, locator,
                                       rng, particle_maker);
  }
  
  /** \brief Maker particles over a time window according to given pulse properties.
  \details All particles are made at once but with different initial times.
   Particles are injected at times ii \c time_step with ii from 0 to ((time_max - time_min)/time_step - 1).
  \param pulse Functor that returns a vector of particles given number of particles and particle maker to make a particle given position.
  \param nr_particles Number of particles to make.
  \param time_min Injection start time.
  \param time_step Injection discretization time step.
  \param time_max Injection stop time.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template
  <typename Pulse,
  typename ParticleMaker>
  auto continuous_injection
  (Pulse const& pulse,
   std::size_t nr_particles,
   double time_min,
   double time_step,
   double time_max,
   ParticleMaker& particle_maker)
  {
    std::size_t nr_injections = (time_max - time_min)/time_step;
    
    using Particle = decltype(particle_maker(Foam::point{}));
    std::vector<Particle> particles;
    particles.reserve(nr_injections*nr_particles);
    
    for (std::size_t ii = 0; ii < nr_injections; ++ii)
    {
      particle_maker.time = ii*time_step;
      for (auto const& part : pulse(nr_particles, particle_maker))
        particles.push_back(part);
    }
    
    return particles;
  }
  
  /** \brief Verify if initial condition is implemented.
   \param initial_condition Name of initial condition.
   \param implemented Associative container of implemented initial condition names. */
  template <typename Implemented>
  void verify_initial_condition
  (std::string const& initial_condition,
   Implemented const& implemented)
  {
    if (!implemented.contains(initial_condition))
      throw std::runtime_error{
        "Initial condition type "
        + initial_condition
        + " not supported"
      };
  }
  
  /** \param position_data Container of position vectors for each particle
   \param particle_maker Functor to make a particle given a position.
   \return Container with particles. */
  template <typename Positions, typename ParticleMaker>
  auto prescribed_positions
  (Positions const& position_data,
   ParticleMaker& particle_maker)
  {
    using Particle = decltype(particle_maker(Foam::point{}));
    std::vector<Particle> particles;
    particles.reserve(position_data.size());
    for (std::size_t pp = 0; pp < position_data.size(); ++pp)
      particles.push_back(particle_maker(make_point(position_data[pp])));
    
    return particles;
  }
  
  /** \param nr_particles Number of particles to make.
  \param filename_position_data Name of file with one position per line; the first \c nr_particles positions are used.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template <typename ParticleMaker>
  auto prescribed_positions
  (std::size_t nr_particles,
   std::string const& filename_position_data,
   ParticleMaker& particle_maker)
  {
    auto input = useful::open_read(filename_position_data);
    
    std::vector<std::vector<double>> position_data;
    position_data.reserve(nr_particles);
    for (std::size_t pp = 0; pp < nr_particles; ++pp)
    {
      std::string line;
      std::getline(input, line);
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of("\t,|\r "));
      boost::algorithm::split(split_line, line, boost::is_any_of("\t,| "),
                              boost::token_compress_on);
      position_data.emplace_back();
      position_data.back().reserve(split_line.size());
      for (auto const& position_component : split_line)
        position_data.back().push_back(std::stod(position_component));
    }
    input.close();
    
    if (nr_particles == 0)
      nr_particles = position_data.size();
    
    return prescribed_positions(position_data, particle_maker);
  }
  
  /**
  \param filename_position_data Name of file with one position per line.
  \param particle_maker Functor to make a particle given a position.
  \return Container with particles. */
  template <typename ParticleMaker>
  auto prescribed_positions
  (std::string const& filename_position_data,
   ParticleMaker& particle_maker)
  {
    auto input = useful::open_read(filename_position_data);
    std::vector<std::vector<double>> position_data;
    std::string line;
    while (std::getline(input, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of("\t,|\r "));
      boost::algorithm::split(split_line, line, boost::is_any_of("\t,| "),
                              boost::token_compress_on);
      position_data.emplace_back();
      position_data.back().reserve(split_line.size());
      for (auto const& position_component : split_line)
        position_data.back().push_back(std::stod(position_component));
    }
    input.close();
    
    return prescribed_positions(position_data, particle_maker);
  }
  
  /** \class InitialCondition_Cases PTOF/InitialConditions.h "PTOF/InitialConditions.h"
   *  \brief InitialCondition object to handle implemented initial condition types. */
  template
  <typename Geometry,
  typename VelocityField,
  typename ParticleMaker,
  typename Parameters,
  typename Mask = useful::Empty>
  class InitialCondition_Cases
  {
    typename Geometry::Mesh const& mesh;                /**< Underlying mesh.                  */
    typename Geometry::Locator const& locator;          /**< Object to search mesh.            */
    VelocityField const& velocity_field;                /**< Underlying flow field.            */
    ParticleMaker particle_maker;                       /**< Makes a particle given a position.*/
    Parameters params;                                  /**< Initial condition parameters. */
    std::mt19937 rng{ std::random_device{}() };         /**< Random number generator.          */
    Mask mask{ useful::Empty{} };                       /**< Mask to disallow cells/faces        */
    double threshold{ 0. };                        /**< Threshold for mask        */
    
  private:
    struct PositionMaker
    {
      double time = 0.;
      auto operator()(Foam::point const& position)
      { return position; }
      
    } _position_maker;                                /** Functor to forward positions and hold time value.*/
 
  public:
    /** Constructor.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field as a function of state.
     \param particle_maker Functor to make a particle given a position.
     \param params Output parameters. */
    InitialCondition_Cases
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     ParticleMaker particle_maker,
     Parameters const& params)
    : mesh{ geometry.mesh() }
    , locator{ geometry.locator }
    , velocity_field{ velocity_field }
    , particle_maker{ particle_maker }
    , params{ params }
    {}
    
    /** Constructor.
    \param geometry Domain geometry info and utilities.
    \param velocity_field Velocity field as a function of state.
    \param particle_maker Functor to make a particle given a position.
    \param params Output parameters.
    \param mask Scalar field assigning values to mesh cell indices through operator[].
    \param threshold Threshold for the mask, such that cells where the mask is above the threshold are considered. */
    InitialCondition_Cases
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     ParticleMaker particle_maker,
     Parameters const& params,
     Mask const& mask,
     double threshold = 0.)
    : mesh{ geometry.mesh() }
    , locator{ geometry.locator }
    , velocity_field{ velocity_field }
    , particle_maker{ particle_maker }
    , params{ params }
    , mask{ mask }
    , threshold{ threshold }
    {}
    
    /** \brief Make particles according to prescribed initial condition and particle maker.
     \param nr_particles Number of particles to make
     \return Container with particles.  */
    auto operator()(std::size_t nr_particles)
    {
      return make_particles(nr_particles, params.type, particle_maker);
    }
    
    /** \brief Make particles according to prescribed initial condition and particle maker, with \c nr_particles particles as specified in \c params.
    \return Container with particles.  */
    auto operator()()
    {
      return this->operator()(params.nr_particles);
    }
    
    /** \brief Make a single particle according to prescribed initial condition and particle maker.
    \return One particle.  */
    auto make_particle()
    {
      return this->operator()(1)[0];
    }
    
    /** \brief Make positions according to prescribed initial condition.
    \param nr_positions Number of positions to make
    \return Container with position.  */
    auto make_positions(std::size_t nr_positions)
    {
      // This works by slightly abusing make_particles with a position maker
      // instead of a particle maker, in order to produce only positions instead of particles,
      // but according to the same prescription as for producing particles.
      return make_particles(nr_positions, params.type,
                            _position_maker);
    }
    
    /** \brief Make a single position according to prescribed initial condition and particle maker.
    \return One position.  */
    auto make_position()
    {
      return make_positions(1)[0];
    }
    
  private:
    /** \brief Make given number of particles.
     \param nr_particles Number of particles to make.
     \param type Initial condition type.
     \param particle_maker Functor to make a particle given a position.
     \return Container with particles.
     */
    template <typename ParticleMakerOther>
    auto make_particles
    (std::size_t nr_particles,
     InitialConditions::Type type,
     ParticleMakerOther& particle_maker)
    {
      switch (type)
      {
        case InitialConditions::Type::uniform:
          if constexpr (std::is_same_v<Mask, useful::Empty>)
            return uniform(nr_particles, mesh, rng, particle_maker);
          else
            return uniform(nr_particles, mesh, rng, particle_maker,
                           mask, threshold);
        case InitialConditions::Type::fluxweighted:
          if constexpr (std::is_same_v<Mask, useful::Empty>)
            return fluxweighted(nr_particles,
                                 velocity_field,
                                 mesh, rng, particle_maker);
          else
            return fluxweighted(nr_particles,
                                 velocity_field,
                                 mesh, rng, particle_maker,
                                 mask, threshold);
        case InitialConditions::Type::uniform_inlet:
          if constexpr (std::is_same_v<Mask, useful::Empty>)
            return uniform_patches(nr_particles,
                                   { "inlet" },
                                   mesh, locator,
                                   rng, particle_maker);
          else
            return uniform_patches(nr_particles,
                                   { "inlet" },
                                   mesh, locator,
                                   rng, particle_maker,
                                   mask, threshold);
        case InitialConditions::Type::fluxweighted_inlet:
          if constexpr (std::is_same_v<Mask, useful::Empty>)
            return fluxweighted_patches(nr_particles,
                                         { "inlet" },
                                         velocity_field,
                                         mesh, locator,
                                         rng, particle_maker);
          else
            return fluxweighted_patches(nr_particles,
                                         { "inlet" },
                                         velocity_field,
                                         mesh, locator,
                                         rng, particle_maker,
                                         mask, threshold);
        case InitialConditions::Type::uniform_solid:
          if constexpr (std::is_same_v<Mask, useful::Empty>)
            return uniform_patches(nr_particles,
                                   { "wallFluidSolid" },
                                   mesh, locator,
                                   rng, particle_maker);
          else
            return uniform_patches(nr_particles,
                                   { "wallFluidSolid" },
                                   mesh, locator,
                                   rng, particle_maker,
                                   mask, threshold);
        case InitialConditions::Type::uniform_near_solid:
          if constexpr (std::is_same_v<Mask, useful::Empty>)
            return uniform_near_boundary_patches(nr_particles,
                                                 { "wallFluidSolid" },
                                                 params.distance_wall,
                                                 mesh, locator,
                                                 rng, particle_maker);
          else
            return uniform_near_boundary_patches(nr_particles,
                                                 { "wallFluidSolid" },
                                                 params.distance_wall,
                                                 mesh, locator,
                                                 rng, particle_maker,
                                                 mask, threshold);
        case InitialConditions::Type::uniform_region_cartesian:
          if constexpr (std::is_same_v<Mask, useful::Empty>)
            return uniform_cells(nr_particles,
                                 cell_ids_region_cartesian(params.region_boundaries,
                                                           mesh, locator),
                                 mesh, rng, particle_maker);
          else
            return uniform_cells(nr_particles,
                                 apply_mask(cell_ids_region_cartesian(params.region_boundaries,
                                                                      mesh, locator),
                                            mask, threshold),
                                 mesh, rng, particle_maker);
        case InitialConditions::Type::fluxweighted_region_cartesian:
          if constexpr (std::is_same_v<Mask, useful::Empty>)
            return fluxweighted_cells(nr_particles,
                                       cell_ids_region_cartesian(params.region_boundaries,
                                                                 mesh, locator),
                                       velocity_field,
                                       mesh,
                                       rng, particle_maker);
          else
            return fluxweighted_cells(nr_particles,
                                       apply_mask(cell_ids_region_cartesian(params.region_boundaries,
                                                                            mesh, locator),
                                                  mask, threshold),
                                       velocity_field,
                                       mesh,
                                       rng, particle_maker);
        case InitialConditions::Type::prescribed_positions:
          return prescribed_positions(nr_particles,
                                      params.position_data,
                                      particle_maker);
        case InitialConditions::Type::uniform_inlet_continuous:
          if constexpr (std::is_same_v<Mask, useful::Empty>)
            return continuous_injection([this]
                                        (std::size_t nr_particles,
                                         ParticleMakerOther& particle_maker)
                                        { return uniform_patches(nr_particles,
                                                                 { "inlet" },
                                                                 mesh, locator,
                                                                 rng,
                                                                 particle_maker); },
                                        nr_particles,
                                        params.time_min,
                                        params.time_max,
                                        params.time_step,
                                        particle_maker);
          else
            return continuous_injection([this]
                                        (std::size_t nr_particles,
                                         ParticleMakerOther& particle_maker)
                                        { return uniform_patches(nr_particles,
                                                                 { "inlet" },
                                                                 mesh, locator,
                                                                 rng,
                                                                 particle_maker,
                                                                 mask, threshold); },
                                        nr_particles,
                                        params.time_min,
                                        params.time_max,
                                        params.time_step,
                                        particle_maker);
        case InitialConditions::Type::fluxweighted_inlet_continuous:
          if constexpr (std::is_same_v<Mask, useful::Empty>)
            return continuous_injection([this]
                                        (std::size_t nr_particles,
                                         ParticleMakerOther& particle_maker)
                                        { return fluxweighted_patches(nr_particles,
                                                                       { "inlet" },
                                                                       velocity_field,
                                                                       mesh, locator,
                                                                       rng,
                                                                       particle_maker); },
                                        nr_particles,
                                        params.time_min,
                                        params.time_max,
                                        params.time_step,
                                        particle_maker);
          else
            return continuous_injection([this]
                                        (std::size_t nr_particles,
                                         ParticleMakerOther& particle_maker)
                                        { return fluxweighted_patches(nr_particles,
                                                                       { "inlet" },
                                                                       velocity_field,
                                                                       mesh, locator,
                                                                       rng,
                                                                       particle_maker,
                                                                       mask, threshold); },
                                        nr_particles,
                                        params.time_min,
                                        params.time_max,
                                        params.time_step,
                                        particle_maker);
        default:
          throw std::runtime_error{
            std::string{ "Initial condition type " }
            + InitialConditions::name(params.type)
            + " not supported" };
      }
    }
    
  public:
    /** \brief Output information about current object. */
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Initial condition\n"
        "--------------------------------------------------\n"
        "Type: " + InitialConditions::name(params.type) + "\n"
        "--------------------------------------------------\n";
    }
  };
  template
  <typename Geometry,
  typename VelocityField,
  typename ParticleMaker,
  typename Parameters,
  typename Mask>
  InitialCondition_Cases
  (Geometry const&,
   VelocityField const&,
   ParticleMaker,
   Parameters const&,
   Mask const&,
   double) ->
  InitialCondition_Cases
  <Geometry, VelocityField, ParticleMaker, Parameters, Mask const&>;
  template
  <typename Geometry,
  typename VelocityField,
  typename ParticleMaker,
  typename Parameters,
  typename Mask>
  InitialCondition_Cases
  (Geometry const&,
   VelocityField const&,
   ParticleMaker,
   Parameters const&,
   Mask const&) ->
  InitialCondition_Cases
  <Geometry, VelocityField, ParticleMaker, Parameters, Mask const&>;
}

#endif /* PTOF_INITIALCONDITIONS_H */
