/**
   \file PTOF/Criterion.h
   \author Tomás Aquino
   \date 07/03/2022
   \brief Criteria to check for different scenarios.
*/

#ifndef PTOF_CRITERION_H
#define PTOF_CRITERION_H

#include "PTOF/Useful.h"
#include <algorithm>
#include <cstddef>
#include <memory>
#include <vector>

namespace ptof {
/**
   \class Criterion PTOF/Criterion.h "PTOF/Criterion.h"
   \brief Abstract polymorphic functor to implement criteria.
*/
template <typename Subject> struct Criterion {
  virtual ~Criterion() {}

  virtual bool operator()(double) const = 0;

protected:
  Criterion(Subject const &subject) : _subject{subject} {}

  Subject const &_subject;
};

/**
   \class Criterion_time PTOF/Criterion.h "PTOF/Criterion.h
   \brief Check if time is greater than value.
*/
template <typename Subject> struct Criterion_time final : Criterion<Subject> {
  Criterion_time(Subject const &subject, double end_value)
      : Criterion<Subject>(subject), end_value{end_value} {}

  bool operator()(double time) const override { return time >= end_value; }

  double end_value;
};

/**
   \class Criterion_mass_below PTOF/Criterion.h "PTOF/Criterion.h"
   \brief Check if mass is less than or equal to value.
*/
template <typename Subject>
struct Criterion_mass_below final : Criterion<Subject> {
  Criterion_mass_below(Subject const &subject, double end_value)
      : Criterion<Subject>(subject), end_value{end_value} {}

  bool operator()(double time) const override {
    return mass(this->_subject, time) <= end_value;
  }

  double end_value;
};

/**
   \class Criterion_mass_above PTOF/Criterion.h "PTOF/Criterion.h"
   \brief Check if mass is greater than or equal to value.
*/
template <typename Subject>
struct Criterion_mass_above final : Criterion<Subject> {
  Criterion_mass_above(Subject const &subject, double end_value)
      : Criterion<Subject>(subject), end_value{end_value} {}

  bool operator()(double time) const override {
    return mass(this->_subject, time) >= end_value;
  }

  double end_value;
};

/**
   \class Criterion_all_absorbed PTOF/Criterion.h "PTOF/Criterion.h"
   \brief Check if all particles have been absorbed.
*/
template <typename Subject>
struct Criterion_all_absorbed final : Criterion<Subject> {
  Criterion_all_absorbed(Subject const &subject)
      : Criterion<Subject>(subject) {}

  bool operator()(double time) const override {
    return nr_absorbed(this->_subject, time) == this->_subject.size();
  }
};

/**
   \class Criterion_one_absorbed PTOF/Criterion.h "PTOF/Criterion.h"
   \brief Check if at least one particle has been absorbed.
*/
template <typename Subject>
struct Criterion_one_absorbed final : Criterion<Subject> {
  Criterion_one_absorbed(Subject const &subject)
      : Criterion<Subject>(subject) {}

  bool operator()(double time) const {
    return nr_absorbed(this->_subject, time) > 0;
  }
};

/**
   \class Criterion_fraction_not_absorbed PTOF/Criterion.h "PTOF/Criterion.h"
   \brief Check if the fraction of particles that have not been absorbed is at
   most a certain value.
*/
template <typename Subject>
struct Criterion_fraction_not_absorbed final : Criterion<Subject> {
  Criterion_fraction_not_absorbed(Subject const &subject, double end_value)
      : Criterion<Subject>(subject), end_value{end_value} {}

  bool operator()(double time) const override {
    return this->_subject.size() - nr_absorbed(this->_subject, time) <=
           std::size_t(end_value * this->_subject.size());
  }

  double end_value;
};

/**
   \class Criterion_and PTOF/Criterion.h "PTOF/Criterion.h"
   \brief Check if all criteria are true.
*/
template <typename Subject> struct Criterion_and final : Criterion<Subject> {
  using Criteria = std::vector<std::unique_ptr<Criterion<Subject>>>;

  Criterion_and(Subject const &subject, Criteria criteria)
      : Criterion<Subject>{subject}, _criteria{std::move(criteria)} {}

  bool operator()(double time) const override {
    return std::all_of(_criteria.begin(), _criteria.end(),
                       [this, time](auto &cc) { return (*cc)(time); });
  }

private:
  Criteria _criteria;
};

/**
   \class Criterion_or PTOF/Criterion.h "PTOF/Criterion.h"
   \brief Check if at least one criterion is true.
*/
template <typename Subject> struct Criterion_or final : Criterion<Subject> {
  using Criteria = std::vector<std::unique_ptr<Criterion<Subject>>>;

  Criterion_or(Subject const &subject, Criteria criteria)
      : Criterion<Subject>{subject}, _criteria{std::move(criteria)} {}

  bool operator()(double time) const override {
    return std::any_of(_criteria.begin(), _criteria.end(),
                       [this, time](auto &cc) { return (*cc)(time); });
  }

private:
  Criteria _criteria;
};

} // namespace ptof

#endif /* PTOF_CRITERION_H */
