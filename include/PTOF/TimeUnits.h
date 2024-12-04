/**
   \file PTOF/TimeUnits.h
   \author Tomás Aquino
   \date 30/09/2024
   \brief Time unit types.
*/

#ifndef PTOF_TIMEUNITS_H
#define PTOF_TIMEUNITS_H

#include <map>
#include <stdexcept>
#include <string>

namespace ptof {
/**
   \struct TimeUnits PTOF/TimeUnits.h "PTOF/TimeUnits.h"
   \brief Keep track of names and types of time units.
*/
struct TimeUnits {
  /**
     \enum Type
     \brief Implemented types.
  */
  enum class Type {
    diffusion, /**< Units of diffusion scale */
    advection, /**< Units of advection scale */
    reaction,  /**<  Units of reaction scale*/
    arbitrary  /**< Arbitrary units*/
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

/** \brief Get time units factor. */
template <typename TransportParameters, typename ReactionParameters>
double time_unit_factor(std::string const &time_units,
                        TransportParameters const &params_transport,
                        ReactionParameters const &params_reaction) {
  switch (TimeUnits::type(time_units)) {
  case TimeUnits::Type::diffusion: {
    return params_transport.diffusion_time;
  }
  case TimeUnits::Type::advection: {
    return params_transport.advection_time;
  }
  case TimeUnits::Type::reaction: {
    return params_reaction.reaction_time;
  }
  case TimeUnits::Type::arbitrary: {
    return 1.;
  }
  default:
    throw std::runtime_error{std::string("Time units ") + time_units +
                             " not supported"};
  }
}
} // namespace ptof

#endif /* PTOF_TIMEUNITS_H */
