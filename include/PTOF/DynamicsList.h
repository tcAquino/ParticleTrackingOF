/**
 * @file   DynamicsList.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Thu Feb 17 00:00:00 2022
 *
 * @brief Helpers for choosing dynamics type.
 */

#ifndef PTOF_DYNAMICSLIST_H
#define PTOF_DYNAMICSLIST_H

#include <map>
#include <string>

namespace ptof {
/**
 * @brief Keep track of names of dynamics types for different types of
 *        simulations.
 */
struct DynamicsList {
  /** @brief Implemented dynamics types. */
  enum class Type {
    transport,    /**< No parameters, custom boundary does nothing. */
    firstpassage, /**< Parameters should hold InitialCondition object. */
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
      {"transport", Type::firstpassage},
      {"transportema", Type::firstpassage},
  };

  /** Map of types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::transport, "transport"},
      {Type::firstpassage, "firstpassage"},
  };
};
} // namespace ptof

#endif /* PTOF_DYNAMICSLIST_H */
