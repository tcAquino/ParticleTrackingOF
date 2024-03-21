/**
 \file PTOF/Reaction.h
 \author Tomás Aquino
 \date 10/03/2022
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
  /**
   \param surface_concentration Homogeneous surface concentration value.
   \param face_ids Mesh face indices where to assign concentration.
   \param mesh OpenFOAM mesh.
   \return Map of homogeneous surface concentrations over given face indices. */
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
  
  /**
  \param surface_concentration Mean surface concentration value.
  \param patch_names Mesh patch names where to assign concentration.
  \param mesh OpenFOAM mesh.
  \return Map of homogeneous surface concentrations over given face indices. */
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
  
  /**
   \param position Spatial position.
   \param locator Object to locate positions in mesh.
   \return Pair of nearest boundary face index and distance to it. */
  template <typename Locator>
  std::pair<Foam::label, double> nearest_boundary_face_dist
  (Foam::vector const& position, Locator const& locator)
  {
    auto const& mesh = locator.mesh_search().mesh();
    // Do not use hint because this can get stuck in local minima
    auto face_id = locator.mesh_search().findNearestBoundaryFace(position);
    auto dist = Foam::mag(position - mesh.faces()[face_id].centre(mesh.points()));

    return { face_id, dist };
  }
  
  /** \class SurfaceReaction_AFluidPlusASolidtoASolid PTOF/Reaction.h "PTOF/Reaction.h"
   * \brief A_F + A_S -> A_S surface reaction. */
  class SurfaceReaction_AFluidPlusASolidtoASolid
  {
  public:
    /** Container to hold reactant surface concentrations. */
    using SurfaceConcentrations = std::unordered_map<Foam::label, double>;
    
    /** Constructor.
      \param rate_constant Surface reaction rate per solid concentration.
      \param diff_coeff Diffusion coefficient.
      \param surface_concentrations Map of initial reactive surface concentrations over mesh faces. */
    SurfaceReaction_AFluidPlusASolidtoASolid
    (double rate_constant,
     double diff_coeff,
     SurfaceConcentrations surface_concentrations)
    : rate_constant{ rate_constant }
    , diff_coeff{ diff_coeff }
    , surface_concentrations{ surface_concentrations }
    {}
    
    /** \param state State of particle to react.
     \param state_old Previous state of particle.
     \param face Mesh face index where reaction occurs.
     \brief React over exposure time. */
    template <typename State>
    void operator()(State& state, State const& state_old, Foam::label face)
    {
      double exposure_time = state.time - state_old.time;
      double probability_first_order =
        rate(face)*std::sqrt(constants::pi*exposure_time/diff_coeff);
      double probability_second_order = probability_first_order/(1.+probability_first_order/2.);
      state.mass *= std::exp(-probability_second_order);
    }
      
    /** \param face Mesh face index.
     \return Surface reaction rate at face */
    double rate(Foam::label face) const
    {
      if (surface_concentrations.count(face))
        return rate_constant*
          surface_concentrations.at(face);
        
      return 0.;
    }
    
    /** \brief Output generic information about object. */
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
    
    double rate_constant;                         /**< Surface reaction rate per solid concentration. */
    double diff_coeff;                            /**< Diffusion coefficient. */
    SurfaceConcentrations surface_concentrations; /**< Map of concentrations over reactive faces. */
  };
  
  /** \class SurfaceReaction_DoNothing
   * \brief Surface reaction that does nothing. */
  class SurfaceReaction_DoNothing
  {
  public:
    /** \param state State of particle to react (unused).
    \param state_old Previous state of particle (unused).
    \param face Mesh face index where reaction occurs (unused).
    \brief React over exposure time (do nothing). */
    template <typename State>
    void operator()(State& state, State const& state_old, Foam::label face) const
    {}
    
    /** \param face Mesh face index (unused).
     \return Surface reaction rate at face (always 0.). */
    double rate(Foam::label face) const
    { return 0.; }
    
    /** \brief Output generic information about object. */
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
  
  /** \class BulkReaction_DoNothing PTOF/Reaction.h "PTOF/Reaction.h"
   * \brief Bulk reaction that does nothing. */
  class BulkReaction_DoNothing
  {
  public:
    /** \param state State of particle to react (unused).
    \param exposure_time Exposure time over which to react (unused).
    \brief React over exposure time (do nothing). */
    template <typename State>
    void operator()(State& state, double exposure_time) const
    {}
    
    /** \brief Output generic information about object. */
    template <typename OStream>
    static void info(OStream& output)
    {
      output <<
        "--------------------------------------------------\n"
        "Bulk reaction\n"
        "--------------------------------------------------\n"
        "None\n"
        "--------------------------------------------------\n";
    }
  };
}

#endif /* PTOF_REACTION_H */
