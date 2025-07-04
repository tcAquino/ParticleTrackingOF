/**
   \file PTOF/MeasurementSpacingList.h
   \author Tomas Aquino
   \date 07/03/2022
   \brief Measurement spacing types condition types.
*/

#ifndef PTOF_MEASUREMENTSPACINGLIST_H
#define PTOF_MEASUREMENTSPACINGLIST_H

#include <map>
#include <string>

namespace ptof {
/** \struct MeasureSpacing PTOF/MeasurementSpacingList.h
 * "PTOF/MeasurementSpacingList.h" \brief Keep track of names and types of
 * measure spacing */
struct MeasurementSpacingList {
  /** \enum Type
   *  \brief Implemented types. */
  enum class Type {
    linear, /**< Fixed number of linearly-spaced measurements, starting and
               ending at specified  times */
    log,    /**< Fixed number of log-spaced measurements, starting and ending at
               specified times.      */
    linear_step, /**< At fixed time intervals, starting at specified time.    */
    log_step /**< With fixed factor between times, starting at specified time.
              */
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

  /** Map of names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"linear", Type::linear},
      {"log", Type::log},
      {"linear_step", Type::linear_step},
      {"log_step", Type::linear_step}};

  /** Map of types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::linear, "linear"},
      {Type::log, "log"},
      {Type::linear_step, "linear_step"},
      {Type::log_step, "logstep"}};
};
} // namespace ptof

#endif /* PTOF_MEASUREMENTLIST_H */
