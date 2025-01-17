/**
   \file PTOF/Reaction.h
   \author Tomás Aquino
   \date 10/03/2022
   \brief Objects and utilities for surface and bulk reactions.
*/

#ifndef PTOF_REACTION_H
#define PTOF_REACTION_H

#include "General/Constants.h"
#include "PTOF/Locator.h"
#include "PTOF/Useful.h"
#include <cmath>
#include <fieldTypes.H>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace ptof {
/**
   \brief Append homogeneous surface concentrations over patch faces to map.
   \param face_ids Mesh face indices where to assign concentration.
   \param surface_concentration Homogeneous surface concentration value.
   \param mesh Mesh object.
   \param surface_concentrations Map of surface concentrations over given face
   indices, to append result.
*/
template <typename Container, typename Mesh>
void uniform_solid_reactant(
    Container const &face_ids, double surface_concentration, Mesh const &mesh,
    std::unordered_map<Foam::label, double> &surface_concentrations) {
  for (auto face : face_ids)
    surface_concentrations[face] = surface_concentration;
}

/**
   \param patch_names Mesh patch names where to assign concentration.
   \param surface_concentration_values Surface concentration value in each
   patch.
   \param mesh Mesh object.
   \return Map of homogeneous surface concentrations over patch face indices.
*/
template <typename Mesh>
auto uniform_solid_reactant(
    std::vector<std::string> const &patch_names,
    std::vector<double> const &surface_concentration_values, Mesh const &mesh) {
  std::unordered_map<Foam::label, double> surface_concentrations;
  for (std::size_t ii = 0; ii < patch_names.size(); ++ii)
    uniform_solid_reactant(patch_face_ids({patch_names[ii]}, mesh),
                           surface_concentration_values[ii], mesh,
                           surface_concentrations);
  return surface_concentrations;
}

/**
   \param position Spatial position.
   \param locator Object to locate positions in mesh.
   \return Pair of nearest boundary face index and distance to it.
*/
template <typename Locator>
std::pair<Foam::label, double>
nearest_boundary_face_dist(Foam::vector const &position,
                           Locator const &locator) {
  auto const &mesh = locator.mesh();
  // Do not use hint because this can get stuck in local minima
  auto face_id = locator.mesh_search().findNearestBoundaryFace(position);
  auto dist = Foam::mag(position - mesh.faces()[face_id].centre(mesh.points()));

  return {face_id, dist};
}

/**
   \class SurfaceReaction_AFluidPlusASolidtoASolid PTOF/Reaction.h
   "PTOF/Reaction.h"
   \brief \f$A_F + A_S \to A_S\f$ surface reaction.
*/
class SurfaceReaction_AFluidPlusASolidtoASolid {
public:
  /** Container to hold reactant surface concentrations. */
  using SurfaceConcentrations = std::unordered_map<Foam::label, double>;

  /**
     \brief Constructor.
     \param rate_constant Surface reaction rate per solid concentration.
     \param diff_coeff Diffusion coefficient.
     \param surface_concentrations Map of initial reactive surface
     concentrations over mesh faces.
  */
  SurfaceReaction_AFluidPlusASolidtoASolid(
      double rate_constant, double diff_coeff,
      SurfaceConcentrations surface_concentrations)
      : rate_constant{rate_constant}, diff_coeff{diff_coeff},
        surface_concentrations{surface_concentrations} {}

  /**
     \brief React over exposure time.
     \param state State of particle to react.
     \param state_old Previous state of particle.
     \param face Mesh face index where reaction occurs.
  */
  template <typename State>
  void operator()(State &state, State const &state_old, Foam::label face) {
    double exposure_time = state.time - state_old.time;
    double probability_first_order =
        rate(face) * std::sqrt(constants::pi * exposure_time / diff_coeff);
    double probability_second_order =
        probability_first_order / (1. + probability_first_order / 2.);
    state.mass *= std::exp(-probability_second_order);
  }

  /**
     \param face Mesh face index.
     \return Surface reaction rate at face.
  */
  double rate(Foam::label face) const {
    if (surface_concentrations.count(face))
      return rate_constant * surface_concentrations.at(face);

    return 0.;
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static void info(std::ostream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Surface reaction\n"
           "--------------------------------------------------------------\n"
           "A_Fluid + A_Solid -> A_Solid\n"
           "--------------------------------------------------------------\n";
  }

  double rate_constant; /**< Surface reaction rate per solid concentration. */
  double diff_coeff;    /**< Diffusion coefficient. */
  SurfaceConcentrations
      surface_concentrations; /**< Map of concentrations over reactive faces. */
};

/**
   \class SurfaceReaction_DoNothing PTOF/Reaction.h "PTOF/Reaction.h"
   \brief Surface reaction that does nothing.
*/
class SurfaceReaction_DoNothing {
public:
  /**
     \param state State of particle to react (unused).
     \param state_old Previous state of particle (unused).
     \param face Mesh face index where reaction occurs (unused).
     \brief React over exposure time (do nothing).
  */
  template <typename State>
  void operator()(State &state, State const &state_old,
                  Foam::label face) const {}

  /**
     \param face Mesh face index (unused).
     \return Surface reaction rate at face (always 0.).
  */
  double rate(Foam::label face) const { return 0.; }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static void info(std::ostream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Surface reaction\n"
           "--------------------------------------------------------------\n"
           "None\n"
           "--------------------------------------------------------------\n";
  }
};

/**
   \class BulkReaction_DoNothing PTOF/Reaction.h "PTOF/Reaction.h"
   \brief Bulk reaction that does nothing.
*/
class BulkReaction_DoNothing {
public:
  /**
     \brief React over exposure time (do nothing).
     \param state State of particle to react (unused).
     \param exposure_time Exposure time over which to react (unused).
  */
  template <typename State>
  void operator()(State &state, double exposure_time) const {}

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static void info(std::ostream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Bulk reaction\n"
           "--------------------------------------------------------------\n"
           "None\n"
           "--------------------------------------------------------------\n";
  }
};
} // namespace ptof

#endif /* PTOF_REACTION_H */
