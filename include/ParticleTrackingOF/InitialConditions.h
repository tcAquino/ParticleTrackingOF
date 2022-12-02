//
//  InitialConditions.h
//
//  Created by Tomás Aquino on 24/02/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef InitialConditions_OF_h
#define InitialConditions_OF_h

#include <cstddef>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <unordered_map>
#include <vector>
#include <fvMesh.H>
#include "ParticleTrackingOF/Boundary.h"
#include "ParticleTrackingOF/useful.h"

namespace ptof
{
  // Keep track of names and
  // types of initial conditions
  struct InitialConditions
  {
    enum class Type
    {
      uniform,
      flux_weighted,
      uniform_inlet,
      flux_weighted_inlet,
      uniform_solid,
      uniform_near_solid,
      uniform_inlet_continuous,
      flux_weighted_inlet_continuous,
    };
    
    // Type from name
    static auto type(std::string const& name)
    { return string_to_type.at(name); }
    
    // Name from name
    static auto name(Type type)
    { return type_to_string.at(type); }
    
    // Check if name exists
    static auto contains(std::string const& name)
    { return string_to_type.count(name); }
    
    // Map names to types
    inline static const
    std::unordered_map<std::string, Type> string_to_type
    {
      { "uniform", Type::uniform },
      { "flux_weighted", Type::flux_weighted },
      { "uniform_inlet", Type::uniform_inlet },
      { "flux_weighted_inlet", Type::flux_weighted_inlet },
      { "uniform_solid", Type::uniform_solid },
      { "uniform_near_solid", Type::uniform_near_solid },
      { "uniform_inlet_continuous", Type::uniform_inlet_continuous },
      { "flux_weighted_inlet_continuous", Type::flux_weighted_inlet_continuous }
    };
    
    // Map types to names
    inline static const
    std::unordered_map<Type, std::string> type_to_string
    {
      { Type::uniform, "uniform", },
      { Type::flux_weighted, "flux_weighted" },
      { Type::uniform_inlet, "uniform_inlet" },
      { Type::flux_weighted_inlet, "flux_weighted_inlet" },
      { Type::uniform_solid, "uniform_solid" },
      { Type::uniform_near_solid, "uniform_near_solid" },
      { Type::uniform_inlet_continuous, "uniform_inlet_continuous" },
      { Type::flux_weighted_inlet_continuous, "flux_weighted_inlet_continuous" }
    };
  };
  
  // Make particles randomly uniformly
  // over given mesh cells
  // Place particles at center of cells
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
    // Weigh probabilities with cell volumes
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
  
  // Make particles with probability proportional
  // to velocity magnitude over given mesh cells
  // Place particles at center of cells
  template
  <typename Container,
  typename VelocityField,
  typename Mesh,
  typename RNG,
  typename ParticleMaker>
  auto flux_weighted_cells
  (std::size_t nr_particles,
   Container const& cell_ids,
   VelocityField const& velocity_field,
   Mesh const& mesh,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    // Weigh probabilities with cell volumes
    // and velocity magnitudes at center
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
  
  // Make particle at face center
  // If face center is not within cell owner,
  // warn and try
  // placing particle slight offset towards center
  // or place particle at center
  template
  <typename Mesh,
  typename MeshSearch,
  typename ParticleMaker>
  auto make_particle_at_face_center
  (Foam::label face,
   Mesh const& mesh,
   MeshSearch const& mesh_search,
   ParticleMaker& particle_maker)
  {
    if (face_center_is_in_cell(face, mesh, mesh_search))
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
      auto offset_face_center = offset_inward_face(face, mesh_search);
      if (mesh_search.findCell(offset_face_center) == owner_cell)
      {
        std::cerr << "Success" << std::endl;
        return particle_maker(offset_face_center);
      }
      else
      {
        std::cerr << "Failed. "
                  << "Placing particle at center of owner face"
                  << std::endl;
        return particle_maker(cell_center(owner_cell, mesh));
      }
    }
  }

  // Make particles randomly uniformly
  // over given mesh cell faces
  // Place particles at center of faces
  template
  <typename Container,
  typename Mesh,
  typename MeshSearch,
  typename RNG,
  typename ParticleMaker>
  auto uniform_faces
  (std::size_t nr_particles,
   Container const& face_ids,
   Mesh const& mesh,
   MeshSearch const& mesh_search,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
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
                                                       mesh, mesh_search,
                                                       particle_maker));
    return particles;
  }
  
  // Make particles with probability proportional
  // to velocity magnitude
  // over given mesh cell faces
  // Place particles at center of faces
  template
  <typename Container,
  typename VelocityField,
  typename Mesh,
  typename MeshSearch,
  typename RNG,
  typename ParticleMaker>
  auto flux_weighted_faces
  (std::size_t nr_particles,
   Container const& face_ids,
   VelocityField const& velocity_field,
   Mesh const& mesh,
   MeshSearch const& mesh_search,
   RNG& rng,
   ParticleMaker& particle_maker)
  { 
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
                                                       mesh, mesh_search,
                                                       particle_maker));
    return particles;
  }
  
  // Make particles randomly uniformly
  // at given distance to given mesh cell face centers
  template
  <typename Container,
  typename Mesh,
  typename MeshSearch,
  typename RNG,
  typename ParticleMaker>
  auto uniform_near_boundary_faces
  (std::size_t nr_particles,
   Container const& face_ids,
   double distance,
   Mesh const& mesh,
   MeshSearch const& mesh_search,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
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
      if (mesh_search.findCell(position) == -1
          || mesh_search.findNearestBoundaryFace(position)
            != face_ids[rand])
        continue;
      particles.push_back(particle_maker(position));
    }
    return particles;
  }
  
  // Make particles randomly uniformly
  // over all mesh cells
  // Place particles at center of cells
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
  
  // Make particles with probability proportional
  // to velocity magnitude
  // over all mesh cells
  // Place particles at center of cells
  template
  <typename VelocityField,
  typename Mesh,
  typename RNG,
  typename ParticleMaker>
  auto flux_weighted
  (std::size_t nr_particles,
   VelocityField const& velocity_field,
   Mesh const& mesh,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    return flux_weighted_cells(nr_particles,
                               all_cell_ids(mesh),
                               velocity_field,
                               mesh, rng, particle_maker);
  }
  
  // Make particles randomly uniformly
  // over mesh cell faces in given patches
  // Place particles at center of faces
  template
  <typename Mesh,
  typename MeshSearch,
  typename RNG,
  typename ParticleMaker>
  auto uniform_patches
  (std::size_t nr_particles,
   std::vector<std::string> const& patch_names,
   Mesh const& mesh,
   MeshSearch const& mesh_search,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    return uniform_faces(nr_particles,
                         patches_face_ids(mesh, patch_names),
                         mesh, mesh_search,
                         rng, particle_maker);
  }
  
  // Make particles with probability proportional
  // to velocity magnitude
  // over mesh cell faces in given patches
  // Place particles at center of faces
  template
  <typename VelocityField,
  typename Mesh,
  typename MeshSearch,
  typename RNG,
  typename ParticleMaker>
  auto flux_weighted_patches
  (std::size_t nr_particles,
   std::vector<std::string> const& patch_names,
   VelocityField const& velocity_field,
   Mesh const& mesh,
   MeshSearch const& mesh_search,
   RNG& rng,
   ParticleMaker& particle_maker)
  {
    return flux_weighted_faces(nr_particles,
                               patches_face_ids(mesh, patch_names),
                               velocity_field,
                               mesh, mesh_search,
                               rng, particle_maker);
  }
  
  // Make particles randomly uniformly
  // at given distance to mesh cell face centers
  // over given patches
  template
  <typename Mesh,
  typename MeshSearch,
  typename RNG,
  typename ParticleMaker>
  auto uniform_near_boundary_patches
  (std::size_t nr_particles,
   std::vector<std::string> const& patch_names,
   double distance,
   Mesh const& mesh,
   MeshSearch const& mesh_search,
   RNG& rng, ParticleMaker particle_maker)
  {
    return uniform_near_boundary_faces(nr_particles,
                                       patches_face_ids(mesh, patch_names),
                                       distance,
                                       mesh, mesh_search,
                                       rng, particle_maker);
  }
  
  // Generate pulse injections
  // according to Pulse object
  // within a given time window
  // at a constant time step
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
  
  // Verify if initial condition is implemented
  template <typename Implemented>
  void verify_initial_condition
  (std::string const& initial_condition,
   Implemented const& implemented)
  {
    if (!implemented.contains(initial_condition))
      throw std::runtime_error{
        "Boundary condition type "
        + initial_condition
        + " not supported"
      };
  }
  
  // InitialCondition object to handle implemented
  // initial condition types
  template
  <typename Geometry,
  typename VelocityField,
  typename ParticleMaker,
  typename Parameters>
  class InitialCondition_Cases
  {
    typename Geometry::Mesh const& mesh;                // Underlying mesh
    typename Geometry::MeshSearch const& mesh_search;   // Object to search mesh
    VelocityField const& velocity_field;                // Underlying flow field
    ParticleMaker particle_maker;                       // Makes a particle given a position
    std::mt19937 rng{ std::random_device{}() };         // Random number generator
    
  private:
    struct PositionMaker
    {
      double time = 0.;
      auto operator()(Foam::point const& position)
      { return position; }
    } position_maker;
 
  public:
    // Construct given
    // geometry information,
    // velocity field,
    // particle maker given position,
    // initial condition parameters
    InitialCondition_Cases
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     ParticleMaker particle_maker,
     Parameters const& params)
    : mesh{ geometry.mesh }
    , mesh_search{ geometry.mesh_search }
    , velocity_field{ velocity_field }
    , particle_maker{ particle_maker }
    , params{ params }
    {}
    
    // Make given number of particles
    // according to prescribed initial condition and particle maker
    auto operator()(std::size_t nr_particles)
    {
      return make_particles(nr_particles, params.type, particle_maker);
    }
    
    // Make the number of particles prescribed in the parameters
    auto operator()()
    {
      return this->operator()(params.nr_particles);
    }
    
    // Make one particle
    auto make_particle()
    {
      return this->operator()(1)[0];
    }
    
    // Make given number of positions
    auto make_positions(std::size_t nr_particles)
    {
      return make_particles(nr_particles, params.type,
                            position_maker);
    }
    
    // Make one position
    auto make_position()
    {
      return make_positions(1)[0];
    }
    
  private:
    Parameters params;  // Initial condition parameters
            
    // Make given number of particles
    // according to given initial condition and particle maker
    template <typename ParticleMakerOther>
    auto make_particles
    (std::size_t nr_particles,
     InitialConditions::Type type,
     ParticleMakerOther& particle_maker)
    {
      switch (type)
      {
        case InitialConditions::Type::uniform:
          return uniform(nr_particles, mesh, rng, particle_maker);
        case InitialConditions::Type::flux_weighted:
          return flux_weighted(nr_particles,
                               velocity_field,
                               mesh, rng, particle_maker);
        case InitialConditions::Type::uniform_inlet:
          return uniform_patches(nr_particles,
                                 { "inlet" },
                                 mesh, mesh_search,
                                 rng, particle_maker);
        case InitialConditions::Type::flux_weighted_inlet:
          return flux_weighted_patches(nr_particles,
                                       { "inlet" },
                                       velocity_field,
                                       mesh, mesh_search,
                                       rng, particle_maker);
        case InitialConditions::Type::uniform_solid:
          return uniform_patches(nr_particles,
                                 { "wallFluidSolid" },
                                 mesh, mesh_search,
                                 rng, particle_maker);
        case InitialConditions::Type::uniform_near_solid:
          return uniform_near_boundary_patches(nr_particles,
                                               { "wallFluidSolid" },
                                               params.distance_wall,
                                               mesh, mesh_search,
                                               rng, particle_maker);
        case InitialConditions::Type::uniform_inlet_continuous:
          return continuous_injection([this]
                                      (std::size_t nr_particles,
                                       ParticleMakerOther& particle_maker)
                                      { return uniform_patches(nr_particles,
                                                               { "inlet" },
                                                               mesh, mesh_search,
                                                               rng,
                                                               particle_maker); },
                                      nr_particles,
                                      params.time_min,
                                      params.time_max,
                                      params.time_step,
                                      particle_maker);
        case InitialConditions::Type::flux_weighted_inlet_continuous:
          return continuous_injection([this]
                                      (std::size_t nr_particles,
                                       ParticleMakerOther& particle_maker)
                                      { return flux_weighted_patches(nr_particles,
                                                                     { "inlet" },
                                                                     velocity_field,
                                                                     mesh, mesh_search,
                                                                     rng,
                                                                     particle_maker); },
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
    // Output information about current object
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
}

#endif /* InitialConditions_OF_h */
