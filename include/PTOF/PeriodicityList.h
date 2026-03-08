/**
 * @file   PeriodicityList.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Thu Feb 17 00:00:00 2022
 *
 * @brief Helpers for choosing periodicity type.
 */

#ifndef PTOF_PERIODICITYLIST_H
#define PTOF_PERIODICITYLIST_H

#include <map>
#include <string>

namespace ptof {
/** @brief Keep track of names of periodicity types for boundary conditions. */
struct PeriodicityList {
  /** @brief Implemented periodicity types. */
  enum class Type {
    cartesian,     /**< Cartesian periodicity. */
    symmetryplanes /**< Periodicity according to symmetry planes. */
  };

  /** @brief Get type from name. */
  static auto type(std::string const &name) { return name_to_type.at(name); }

  /** @brief Get name from type. */
  static auto name(Type type) { return type_to_name.at(type); }

  /** @brief Check if name exists. */
  static bool contains(std::string const &name) {
    return name_to_type.count(name);
  }

  /** Map of names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"cartesian", Type::cartesian}, {"symmetryplanes", Type::symmetryplanes}};

  /** Map of types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::cartesian, "cartesian"}, {Type::symmetryplanes, "symmetryplanes"}};
};
} // namespace ptof

#endif /* PTOF_PERIODICITYLIST_H */
