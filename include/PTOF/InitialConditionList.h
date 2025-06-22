/**
   \file PTOF/InitialConditionList.h
   \author Tomas Aquino
   \date 29/09/2024
   \brief Initial condition types.
*/

#ifndef PTOF_INITIALCONDITIONLIST_H
#define PTOF_INITIALCONDITIONLIST_H

#include <map>
#include <string>

namespace ptof {
/**
   \struct InitialConditionList PTOF/InitialConditionList.h
   "PTOF/InitialConditionList.h"
   \brief Keep track of names and types of initial conditions.
*/
struct InitialConditionList {
  /**
     \enum Type
     \brief Implemented types.
  */
  enum class Type {
    point,                    /**< At specified position. */
    uniform_all_cells,        /**< Homogeneous throughout the domain. */
    fluxweighted_all_cells,   /**< Flux-weighted throughout the domain. */
    uniform_patch_faces,      /**< Homogeneous over boundary patch_faces. */
    fluxweighted_patch_faces, /**< Flux-weighted over boundary patch_faces. */
    uniform_near_patch,       /**< Homogeneous at a fixed distance to boundary
                                   patch_faces. */
    fluxweighted_near_patch,  /**< Flux-weighted at a fixed distance to
                              boundary patch_faces. */
    uniform_cartesian_cells,  /**< Homogeneous in a Cartesian region. */
    fluxweighted_cartesian_cells, /**< Flux-weighted in a Cartesian region. */
    prescribed_positions,         /**< Prescribed particle positions. */
    prescribed_positions_masses,  /**< Prescribed particle positions and masses
                                   */
    prescribed_positions_masses_tags,       /**< Prescribed particle positions,
                                               masses, and tags. */
    prescribed_positions_masses_tags_times, /**< Prescribed particle positions,
                                               masses, tags, and times. */
    uniform_cylinder,       /** Uniform within a cylinder, not tied to cells. */
    uniform_parallelipiped, /** Uniform within a parallelipiped, not tied to
                               cells. */
    uniform_sphere          /** Uniform within a sphere, not tied to cells. */
  };

  /**
     \brief Type from name.
     \param name Boundary condition name.
     \return Boundary condition type.
  */
  static auto type(std::string const &name) { return name_to_type.at(name); }

  /**
     \brief Name from type.
     \param type Boundary condition type.
     \return Boundary condition name.
  */
  static auto name(Type type) { return type_to_name.at(type); }

  /**
     \brief Check if name exists.
     \param name condition name.
     \return \c true if name exists, \c false otherwise.
  */
  static bool contains(std::string const &name) {
    return name_to_type.count(name);
  }

  /** Map names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"point", Type::point},
      {"uniform_all_cells", Type::uniform_all_cells},
      {"fluxweighted_all_cells", Type::fluxweighted_all_cells},
      {"uniform_patch_faces", Type::uniform_patch_faces},
      {"fluxweighted_patch_faces", Type::fluxweighted_patch_faces},
      {"uniform_near_patch", Type::uniform_near_patch},
      {"fluxweighted_near_patch", Type::fluxweighted_near_patch},
      {"uniform_cartesian_cells", Type::uniform_cartesian_cells},
      {"fluxweighted_cartesian_cells", Type::fluxweighted_cartesian_cells},
      {"prescribed_positions", Type::prescribed_positions},
      {"prescribed_positions_masses", Type::prescribed_positions_masses},
      {"prescribed_positions_masses_tags",
       Type::prescribed_positions_masses_tags},
      {"prescribed_positions_masses_tags_times",
       Type::prescribed_positions_masses_tags_times}};

  /** Map types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::point, "point"},
      {Type::uniform_all_cells, "uniform_all_cells"},
      {Type::fluxweighted_all_cells, "fluxweighted_all_cells"},
      {Type::uniform_patch_faces, "uniform_patch_faces"},
      {Type::fluxweighted_patch_faces, "fluxweighted_patch_faces"},
      {Type::uniform_near_patch, "uniform_near_patch"},
      {Type::fluxweighted_near_patch, "fluxweighted_near_patch"},
      {Type::uniform_cartesian_cells, "uniform_cartesian_cells"},
      {Type::fluxweighted_cartesian_cells, "fluxweighted_cartesian_cells"},
      {Type::prescribed_positions, "prescribed_positions"},
      {Type::prescribed_positions_masses, "prescribed_positions_masses"},
      {Type::prescribed_positions_masses_tags,
       "prescribed_positions_masses_tags"},
      {Type::prescribed_positions_masses_tags_times,
       "Prescribed_Positions_Masses_Tags_Times"}};
};
} // namespace ptof

#endif /* PTOF_INITIALCONDITIONLIST_H */
