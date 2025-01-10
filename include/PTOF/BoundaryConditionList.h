/**
   \file PTOF/BoundaryConditionList.h
   \author Tomás Aquino
   \date 29/09/2024
   \brief Boundary condition types.
*/

#ifndef PTOF_BOUNDARYCONDITIONLIST_H
#define PTOF_BOUNDARYCONDITIONLIST_H

#include <map>
#include <string>

namespace ptof {
/**
   \struct BoundaryConditionList PTOF/BoundaryConditionList.h
   "PTOF/BoundaryConditionList.h"
   \brief Names and types of boundary conditions.
*/
struct BoundaryConditionList {
  /**
     \enum Type
     \brief Implemented types.
  */
  enum class Type {
    reflecting, /**< Reflecting                                         */
    reacting_reflecting, /**< Reacting and reflecting                   */
    periodic,  /**< Periodic                                            */
    absorbing, /**< Absorbing                                           */
    info,      /**< Information upon hitting                            */
    inlet, /**< Reflecting for non-negative flow velocities, absorbing otherwise
            */
    custom, /**< Enforced by custom Boundary object                  */
    empty   /**< No effect (default for unspecified patches in mesh) */
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

  /** \brief Map names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"reflecting", Type::reflecting},
      {"reacting_reflecting", Type::reacting_reflecting},
      {"periodic", Type::periodic},
      {"absorbing", Type::absorbing},
      {"inlet", Type::inlet},
      {"custom", Type::custom},
      {"empty", Type::empty},
  };

  /** \brief Map types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::reflecting, "reflecting"},
      {Type::reacting_reflecting, "reacting_reflecting"},
      {Type::periodic, "periodic"},
      {Type::absorbing, "absorbing"},
      {Type::inlet, "inlet"},
      {Type::custom, "custom"},
      {Type::empty, "empty"}};
};
} // namespace ptof

#endif /* PTOF_BOUNDARYCONDITIONLIST_H */
