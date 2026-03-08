/**
 * @file   EndCriterionList.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Thu May  8 00:00:00 2025
 *
 * @brief Criteria to stop dynamics.
 */

#include <map>
#include <string>

#ifndef PTOF_ENDCRITERIONLIST_H
#define PTOF_ENDCRITERIONLIST_H

namespace ptof {
/**
 * @brief Keep track of names and types of end criteria.
 */
struct EndCriterionList {
  /** @brief Implemented types. */
  enum class Type {
    time,                 /**< Specified time. */
    time_max,             /**< Maximum output time. */
    mass_below,           /**< Total mass below value. */
    mass_above,           /**< Total mass above value. */
    all_absorbed,         /**< All particles absorbed. */
    one_absorbed,         /**< One particle absorbed. */
    fraction_not_absorbed /**< Particle fraction not absorbed. */
  };

  /** @return Type from name. */
  static auto type(std::string const &name) { return name_to_type.at(name); }

  /** @return Name from type. */
  static auto name(Type type) { return type_to_name.at(type); }

  /** @return \c true if name exists, \c false otherwise. */
  static bool contains(std::string const &name) {
    return name_to_type.count(name);
  }

  /** Map of names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"time", Type::time},
      {"time_max", Type::time_max},
      {"mass_below", Type::mass_below},
      {"mass_above", Type::mass_above},
      {"all_absorbed", Type::all_absorbed},
      {"one_absorbed", Type::one_absorbed},
      {"fraction_not_absorbed", Type::fraction_not_absorbed}};

  /** Map of types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::time, "time"},
      {Type::time_max, "time_max"},
      {Type::mass_below, "mass_below"},
      {Type::mass_above, "mass_above"},
      {Type::all_absorbed, "all_absorbed"},
      {Type::one_absorbed, "one_absorbed"},
      {Type::fraction_not_absorbed, "fraction_not_absorbed"}};
};
} // namespace ptof

#endif /* PTOF_CRITERIONLIST_H */
