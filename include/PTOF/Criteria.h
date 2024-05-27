/**
 \file PTOF/Criteria.h
 \author Tomás Aquino
 \date 07/03/2022
*/

#ifndef PTOF_CRITERIA_H
#define PTOF_CRITERIA_H

#include "PTOF/Useful.h"
#include <cstddef>
#include <string>
#include <unordered_map>

namespace ptof {
/** \class Criterion PTOF/Output.h "PTOF/Output.h"
 * \brief  Polymorphic functor to implement criteria
 \note To be used only through derived classes. */
template <typename Subject> struct Criterion {
  virtual ~Criterion() {}

  virtual bool operator()(double) const = 0;

protected:
  Criterion(Subject const &subject) : _subject{subject} {}

  Subject const &_subject;
};

/** \class Criterion_time PTOF/Output.h "PTOF/Output.h
 *  \brief Check if time is greater than value. */
template <typename Subject> struct Criterion_time final : Criterion<Subject> {
  Criterion_time(Subject const &subject, double end_value)
      : Criterion<Subject>(subject), end_value{end_value} {}

  bool operator()(double time) const override { return time >= end_value; }

  double end_value;
};

/** \class Criterion_mass_below PTOF/Output.h "PTOF/Output.h"
 *  \brief Check if mass is less than or equal to value. */
template <typename Subject>
struct Criterion_mass_below final : Criterion<Subject> {
  Criterion_mass_below(Subject const &subject, double end_value)
      : Criterion<Subject>(subject), end_value{end_value} {}

  bool operator()(double time) const override {
    return mass(Criterion<Subject>::_subject, time) <= end_value;
  }

  double end_value;
};

/** \class Criterion_mass_above PTOF/Output.h "PTOF/Output.h"
 *  \brief Check if mass is greater than or equal to value. */
template <typename Subject>
struct Criterion_mass_above final : Criterion<Subject> {
  Criterion_mass_above(Subject const &subject, double end_value)
      : Criterion<Subject>(subject), end_value{end_value} {}

  bool operator()(double time) const override {
    return mass(Criterion<Subject>::_subject, time) >= end_value;
  }

  double end_value;
};

/** \class Criterion_all_absorbed PTOF/Output.h "PTOF/Output.h"
 * \brief Check if all particles have been absorbed. */
template <typename Subject>
struct Criterion_all_absorbed final : Criterion<Subject> {
  Criterion_all_absorbed(Subject const &subject)
      : Criterion<Subject>(subject) {}

  bool operator()(double time) const override {
    return nr_absorbed(Criterion<Subject>::_subject, time) ==
           Criterion<Subject>::_subject.size();
  }
};

/** \class Criterion_one_absorbed PTOF/Output.h "PTOF/Output.h"
 *  \brief Check if at least one particle has been absorbed. */
template <typename Subject>
struct Criterion_one_absorbed final : Criterion<Subject> {
  Criterion_one_absorbed(Subject const &subject)
      : Criterion<Subject>(subject) {}

  bool operator()(double time) const {
    return nr_absorbed(Criterion<Subject>::_subject, time) > 0;
  }
};

/** \class Criterion_fraction_not_absorbed PTOF/Output.h "PTOF/Output.h"
 *  \brief Check if the fraction of particles that have not been absorbed is at
 * most a certain value. */
template <typename Subject>
struct Criterion_fraction_not_absorbed final : Criterion<Subject> {
  Criterion_fraction_not_absorbed(Subject const &subject, double end_value)
      : Criterion<Subject>(subject), end_value{end_value} {}

  bool operator()(double time) const override {
    return Criterion<Subject>::_subject.size() -
               nr_absorbed(Criterion<Subject>::_subject, time) <=
           std::size_t(end_value * Criterion<Subject>::_subject.size());
  }

  double end_value;
};

/** \struct EndCriterion PTOF/Output.h "PTOF/Output.h"
 * \brief Keep track of names and types of end criteria. */
struct EndCriterion {
  /** \enum Type
   *  \brief Implemented types. */
  enum class Type {
    time,                 /**< Specified time. */
    time_max,             /**< Maximum output time. */
    mass_below,           /**< Total mass below value. */
    mass_above,           /**< Total mass above value. */
    all_absorbed,         /**< All particles absorbed. */
    one_absorbed,         /**< One particle absorbed. */
    fraction_not_absorbed /**< Particle fraction not absorbed. */
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
  inline static const std::unordered_map<std::string, Type> name_to_type{
      {"time", Type::time},
      {"time_max", Type::time_max},
      {"mass_below", Type::mass_below},
      {"mass_above", Type::mass_above},
      {"all_absorbed", Type::all_absorbed},
      {"one_absorbed", Type::one_absorbed},
      {"fraction_not_absorbed", Type::fraction_not_absorbed}};

  /** Map of types to names. */
  inline static const std::unordered_map<Type, std::string> type_to_name{
      {Type::time, "time"},
      {Type::time_max, "time_max"},
      {Type::mass_below, "mass_below"},
      {Type::mass_above, "mass_above"},
      {Type::all_absorbed, "all_absorbed"},
      {Type::one_absorbed, "one_absorbed"},
      {Type::fraction_not_absorbed, "fraction_not_absorbed"}};
};
} // namespace ptof

#endif /* PTOF_CRITERIA_H */
