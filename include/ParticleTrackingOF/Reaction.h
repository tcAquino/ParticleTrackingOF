//
//  Reaction.h
//  ParticleTracking_OpenFOAM
//
//  Created by Tomás Aquino on 10/03/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Reaction_OF_h
#define Reaction_OF_h

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>
#include <fieldTypes.H>
#include "ParticleTrackingOF/useful.h"

namespace ptof
{
  // Make map of homogeneous surface concentrations
  // over given face indices
  template
  <typename Container,
  typename Mesh>
  auto uniform_solid_reactant
  (double surface_concentration,
   Container const& face_ids,
   Mesh const& mesh)
  {                       
    std::unordered_map<Foam::label, double> surface_concentrations;
    for (auto face : face_ids)
      surface_concentrations[face] = surface_concentration;
    return surface_concentrations;
  }
  
  // Make map of uniformly random surface concentrations
  // with given mean
  // over faces in given patches
  template <typename Mesh>
  auto uniform_solid_reactant_patches
  (double surface_concentration,
   std::vector<std::string> const& patch_names,
   Mesh const& mesh)
  {
    return uniform_solid_reactant(surface_concentration,
                                  patches_face_ids(mesh, patch_names),
                                  mesh);
  }
  
  // Find nearest boundary face and distance to it
  template <typename MeshSearch>
  std::pair<Foam::label, double> nearest_boundary_face_dist
  (Foam::vector const& position, MeshSearch const& mesh_search)
  {
    // Do not use hint because this can get stuck in local minima
    auto const& mesh = mesh_search.mesh();
    auto face_id = mesh_search.findNearestBoundaryFace(position);
    auto dist = Foam::mag(position - mesh.faces()[face_id].centre(mesh.points()));

    return { face_id, dist };
  }
  
  // A_F + A_S -> A_S surface reaction
  template <typename MeshSearch>
  class Reaction_AFluidPlusASolidtoASolid
  {
  public:
    // Type of container to hold reactant surface concentrations
    using SurfaceConcentrations = std::unordered_map<Foam::label, double>;
    
    // Construct give well-mixed volumetric reaction rate,
    // threshold distance, map of initial reactive surface concentrations,
    // and mesh search object
    Reaction_AFluidPlusASolidtoASolid
    (double reaction_rate,
     double distance,
     SurfaceConcentrations surface_concentrations,
     MeshSearch&& mesh_search)
    : surface_reaction_rate{ reaction_rate/distance }
    , distance{ distance }
    , surface_concentrations{ surface_concentrations }
    , mesh_search{ std::forward<MeshSearch>(mesh_search) }
    {}
    
    // React over exposure time
    template <typename State>
    auto operator()(State& state, double exposure_time)
    {
      auto nearest = nearest_boundary_face_dist(make_point(state.position),
                                                mesh_search);
      if (nearest.second < distance
          && surface_concentrations.count(nearest.first))
        state.mass *= std::exp(-surface_reaction_rate*
                               surface_concentrations.at(nearest.first)*
                               exposure_time);
    }
    
    // Output generic information about object
    template <typename OStream>
    static void info(OStream& output)
    {
      output <<
        "--------------------------------------------------\n"
        "Reaction\n"
        "--------------------------------------------------\n"
        "A_Fluid + A_Solid -> A_Solid surface reaction\n"
        "--------------------------------------------------\n";
    }
    
    double surface_reaction_rate; // Well-mixed volumetric rate by threshold distance
    double distance;              // Threshold distance for reaction to occur
    SurfaceConcentrations surface_concentrations; // Map of concentrations over reactive faces
    MeshSearch mesh_search;                       // To find mesh faces
  };
  template <typename MeshSearch>
  Reaction_AFluidPlusASolidtoASolid
  (double,
   double,
   std::unordered_map<Foam::label, double>,
   MeshSearch&&) ->
  Reaction_AFluidPlusASolidtoASolid<MeshSearch>;
  
  // A_F + A_S -> Nothing surface reaction
  template <typename MeshSearch>
  class Reaction_AFluidPlusASolidtoNothing
  {
  public:
    // Type of container to hold reactant surface concentrations
    using SurfaceConcentrations = std::unordered_map<Foam::label, double>;
    
    // Construct give well-mixed volumetric reaction rate,
    // threshold distance, map of initial reactive surface concentrations,
    // and mesh search object
    Reaction_AFluidPlusASolidtoNothing
    (double reaction_rate,
     double distance,
     SurfaceConcentrations surface_concentrations,
     MeshSearch&& mesh_search)
    : surface_reaction_rate{ reaction_rate/distance }
    , distance{ distance }
    , surface_concentrations{ surface_concentrations }
    , mesh_search{ std::forward<MeshSearch>(mesh_search) }
    {}
    
    // React over exposure time
    template <typename State>
    auto operator()(State& state, double exposure_time)
    {
      auto nearest = nearest_boundary_face_dist(make_point(state.position),
                                                mesh_search);
      if (nearest.second < distance
          && surface_concentrations.count(nearest.first))
      {
        double aux = surface_reaction_rate*
          surface_concentrations.at(nearest.first)*
          exposure_time;
        state.mass *= std::exp(-aux);
        surface_concentrations[nearest.first] *= std::exp(-aux/
          face_area(nearest.first, mesh_search.mesh()));
      }
    }
    
    // Output generic information about object
    template <typename OStream>
    static void info(OStream& output)
    {
      output <<
        "--------------------------------------------------\n"
        "Reaction\n"
        "--------------------------------------------------\n"
        "A_Fluid + A_Solid -> Nothing surface reaction\n"
        "--------------------------------------------------\n";
    }
    
    double surface_reaction_rate; // Well-mixed volumetric rate by threshold distance
    double distance;              // Threshold distance for reaction to occur
    SurfaceConcentrations surface_concentrations; // Map of concentrations over reactive faces
    MeshSearch mesh_search;                       // To find mesh faces
  };
  template <typename MeshSearch>
  Reaction_AFluidPlusASolidtoNothing
  (double,
   double,
   std::unordered_map<Foam::label, double>,
   MeshSearch&&) ->
  Reaction_AFluidPlusASolidtoNothing<MeshSearch>;
}

#endif /* Reaction_OF_h */
