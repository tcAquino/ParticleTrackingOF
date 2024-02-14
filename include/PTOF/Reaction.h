/**
* \file PTOF/Reaction.h
* \author Tomás Aquino
* \date 10/03/2022
*/

#ifndef PTOF_REACTION_H
#define PTOF_REACTION_H

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>
#include <fieldTypes.H>
#include "General/Constants.h"
#include "PTOF/Locator.h"
#include "PTOF/Useful.h"

namespace ptof
{
  /** Make map of homogeneous surface concentrations
   * over given face indices. */
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
  
  /** Make map of uniformly random surface concentrations
   * with given mean
   * over faces in given patches. */
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
  
  /** Find nearest boundary face and diff_coeff to it. */
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
  
  /** \class SurfaceReaction_AFluidPlusASolidtoASolid PTOF/Reaction.h "PTOF/Reaction.h"
   * \brief A_F + A_S -> A_S surface reaction. */
  class SurfaceReaction_AFluidPlusASolidtoASolid
  {
  public:
    /** Type of container to hold reactant surface concentrations. */
    using SurfaceConcentrations = std::unordered_map<Foam::label, double>;
    
    /** Construct give well-mixed volumetric reaction rate,
     * threshold diff_coeff, map of initial reactive surface concentrations,
     * and mesh search object. */
    SurfaceReaction_AFluidPlusASolidtoASolid
    (double rate_constant,
     double diff_coeff,
     SurfaceConcentrations surface_concentrations)
    : rate_constant{ rate_constant }
    , diff_coeff{ diff_coeff }
    , surface_concentrations{ surface_concentrations }
    {}
    
    /** React over exposure time. */
    template <typename State>
    void operator()(State& state, State const& state_old, Foam::label face)
    {
      double exposure_time = state.time - state_old.time;
      double probability_first_order =
        rate(face)*std::sqrt(constants::pi*exposure_time/diff_coeff);
      double probability_second_order = probability_first_order/(1.+probability_first_order/2.);
      state.mass *= std::exp(-probability_second_order);
    }
      
    /** Surface reaction rate at face */
    double rate(Foam::label face) const
    {
      if (surface_concentrations.count(face))
        return rate_constant*
          surface_concentrations.at(face);
        
      return 0.;
    }
    
    /** Output generic information about object. */
    template <typename OStream>
    static void info(OStream& output)
    {
      output <<
        "--------------------------------------------------\n"
        "Surface reaction\n"
        "--------------------------------------------------\n"
        "A_Fluid + A_Solid -> A_Solid\n"
        "--------------------------------------------------\n";
    }
    
    double rate_constant;                 /**< Surface reaction rate per solid concentrations.*/
    double diff_coeff;                            /**< Diffusion coefficient.*/
    SurfaceConcentrations surface_concentrations; /**< Map of concentrations over reactive faces.*/
  };
  
  /** \class SurfaceReaction_DoNothing
   * \brief No reaction. */
  class SurfaceReaction_DoNothing
  {
  public:
    /** React over exposure time. */
    template <typename State>
    void operator()(State& state, State const& state_old, Foam::label face) const
    {}
    
    /** Surface reaction rate at face */
    double rate(Foam::label face) const
    { return 0.; }
    
    /** Output generic information about object. */
    template <typename OStream>
    static void info(OStream& output)
    {
      output <<
        "--------------------------------------------------\n"
        "Surface reaction\n"
        "--------------------------------------------------\n"
        "None\n"
        "--------------------------------------------------\n";
    }
  };
  
  /** \class Reaction_DoNothing
   * \brief No reaction. */
  class Reaction_DoNothing
  {
  public:
    /** React over exposure time. */
    template <typename State>
    void operator()(State&, double) const
    {}
    
    /** Output generic information about object. */
    template <typename OStream>
    static void info(OStream& output)
    {
      output <<
        "--------------------------------------------------\n"
        "Reaction\n"
        "--------------------------------------------------\n"
        "None\n"
        "--------------------------------------------------\n";
    }
  };
}

#endif /* PTOF_REACTION_H */
