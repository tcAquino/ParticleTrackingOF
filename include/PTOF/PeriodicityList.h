/**
   \file PTOF/PeriodicityList.h
   \author Tomas Aquino
   \date 2022/02/17
   \brief Helpers for choosing
*/

#ifndef PTOF_PERIODICITYLIST_H
#define PTOF_PERIODICITYLIST_H

#include <map>
#include <string>

namespace ptof {
/**
   \struct PeriodicityList PTOF/PeriodicityList.h "PTOF/PeriodicityList.h"
   \brief Keep track of names of periodicity types for boundary conditions.
*/
struct PeriodicityList {
  /** \enum Type PTOF/PeriodicityList.h "PTOF/PeriodicityList.h"
   *  \brief Implemented periodicity types. */
  enum class Type {
    cartesian,     /**< Cartesian periodicity.                         */
    symmetryplanes /**< Periodicity according to symmetry planes.      */
  };

  /** \brief Get type from name. */
  static auto type(std::string const &name) { return name_to_type.at(name); }

  /** \brief Get name from type. */
  static auto name(Type type) { return type_to_name.at(type); }

  /** \brief Check if name exists. */
  static bool contains(std::string const &name) {
    return name_to_type.count(name);
  }

  /** \brief Map names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"cartesian", Type::cartesian}, {"symmetryplanes", Type::symmetryplanes}};

  /** \brief Map types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::cartesian, "cartesian"}, {Type::symmetryplanes, "symmetryplanes"}};
};
} // namespace ptof

#endif /* PTOF_PERIODICITYLIST_H */
