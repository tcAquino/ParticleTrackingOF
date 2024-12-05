/**
 \file PTOF/MeasurementSpacing.h
 \author Tomás Aquino
 \date 07/03/2022
*/

#ifndef PTOF_MEASUREMENTSPACING_H
#define PTOF_MEASUREMENTSPACING_H

#include <map>
#include <string>

namespace ptof {
/** \struct MeasureSpacing PTOF/MeasurementSpacing.h "PTOF/MeasurementSpacing.h"
 *  \brief Keep track of names and types of measure spacing */
struct MeasurementSpacing {
  /** \enum Type
   *  \brief Implemented types. */
  enum class Type {
    linear, /**< Fixed number of linearly-spaced measurements, starting and
               ending at specified  times */
    log,    /**< Fixed number of log-spaced measurements, starting and ending at
               specified times.      */
    step    /**< At fixed time intervals, starting at specified time.    */
  };

  /** \return Type from name. */
  static auto type(std::string const &name) { return name_to_type.at(name); }

  /** \return Name from type. */
  static auto name(Type type) { return type_to_name.at(type); }

  /** \return \c true if name exists, \c false otherwise. */
  static bool exists(std::string const &name) {
    return name_to_type.count(name);
  }

  /** Map of names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"linear", Type::linear}, {"log", Type::log}, {"step", Type::step}};

  /** Map of types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::linear, "linear"}, {Type::log, "log"}, {Type::step, "step"}};
};

/** \struct MeasurementSpacingUnits PTOF/MeasurementSpacing.h "PTOF/MeasurementSpacing.h"
 * \brief Keep track of names and types of measure spacing units (scaling) */
struct MeasurementSpacingUnits {
  /** \enum Type
   *  \brief Implemented types. */
  enum class Type {
    diffusion, /**< Units of diffusion scale */
    advection, /**< Units of advection scale */
    reaction,  /**<  Units of reaction scale*/
    arbitrary  /**< Arbitrary units*/
  };

  /** \return Type from name. */
  static auto type(std::string const &name) { return name_to_type.at(name); }

  /** \return Name from type. */
  static auto name(Type type) { return type_to_name.at(type); }

  /** \return \c true if criterion name exists, \c false otherwise. */
  static bool exists(std::string const &name) {
    return name_to_type.count(name);
  }

  /** Map of names to types. */
  inline static const std::map<std::string, Type> name_to_type{
      {"diffusion", Type::diffusion},
      {"advection", Type::advection},
      {"reaction", Type::reaction},
      {"arbitrary", Type::arbitrary}};

  /** Map of types to names. */
  inline static const std::map<Type, std::string> type_to_name{
      {Type::diffusion, "diffusion"},
      {Type::advection, "advection"},
      {Type::reaction, "reaction"},
      {Type::arbitrary, "arbitrary"}};
};
} // namespace ptof

#endif /* PTOF_MEASUREMENT_H */
