/**
 \file General/Ranges.h
 \author Tomas Aquino
 \date 08/06/2016
*/

// Initialize and set containers to ranges

#ifndef GENERAL_RANGES_H
#define GENERAL_RANGES_H

#include "General/Meta.h"
#include <cmath>
#include <stdexcept>
#include <type_traits>
#include <vector>

/**
   \namespace range Numerical sequences with different inter-element spacings.
*/
namespace range {
/** \return Linear partition of [\p xx1, \p xx2] with \p nr_points points. */
template <template <typename...> typename Container = std::vector,
          typename Scalar, typename... Args>
Container<Scalar, Args...> linspace(Scalar xx1, Scalar xx2,
                                    std::size_t nr_points) {
  if (nr_points == 0) {
    return {};
  }

  Scalar inc = (xx2 - xx1) / (nr_points - 1);

  if constexpr (meta::has_push_back_v<Container<Scalar, Args...>>) {
    Container<Scalar, Args...> output;
    if constexpr (meta::has_reserve_v<Container<Scalar, Args...>>) {
      output.reserve(nr_points);
    }
    output.push_back(xx1);
    for (std::size_t ii = 0; ii < nr_points - 1; ++ii) {
      output.push_back(xx1 + ii * inc);
    }
    output.push_back(xx2);
    return output;
  } else {
    Container<Scalar, Args...> output(nr_points);
    output[0] = xx1;
    for (std::size_t ii = 1; ii < output.size() - 1; ++ii) {
      output[ii] = xx1 + ii * inc;
    }
    output[nr_points - 1] = xx2;

    return output;
  }
}

/** \return log partition of [\p xx1, \p xx2] with \p nr_points points. **/
template <template <typename...> typename Container = std::vector,
          typename Scalar, typename... Args>
Container<Scalar, Args...> logspace(Scalar xx1, Scalar xx2,
                                    std::size_t nr_points) {
  if (nr_points == 0) {
    return {};
  }

  double inc = std::pow(xx2 / xx1, 1. / (nr_points - 1));

  if constexpr (meta::has_push_back_v<Container<Scalar, Args...>>) {
    Container<Scalar, Args...> output;
    if constexpr (meta::has_reserve_v<Container<Scalar, Args...>>) {
      output.reserve(nr_points);
    }
    output.push_back(xx1);
    for (std::size_t ii = 1; ii < nr_points - 1; ++ii) {
      output.push_back(output.back() * inc);
    }
    output.push_back(xx2);
    return output;
  } else {
    Container<Scalar, Args...> output(nr_points);
    output[0] = xx1;
    for (std::size_t ii = 1; ii < output.size() - 1; ++ii) {
      output[ii] = output[ii - 1] * inc;
    }
    output[nr_points - 1] = xx2;
    return output;
  }
}

/**
   \return Linear partition of \c xx1 up to \c xx2 (excluded), with spacing
   given by \p increment.
   \note For signed \c Scalar types, the last element can be smaller than the
   first. The sign of the increment is adjusted automatically.
 */
template <template <typename...> typename Container = std::vector,
          typename Scalar, typename... Args>
Container<Scalar, Args...> range(Scalar xx1, Scalar increment, Scalar xx2) {
  if constexpr (std::is_unsigned_v<Scalar>) {
    if (xx2 < xx1) {
      throw std::runtime_error{
          "Last element must be greater than first for unisgned types"};
    }
  }
  increment = xx2 < xx1 ? -std::abs(increment) : std::abs(increment);

  std::size_t nr_points = (xx2 - xx1) / increment;
  if (nr_points == 0) {
    return {};
  }

  if constexpr (meta::has_push_back_v<Container<Scalar, Args...>>) {
    Container<Scalar, Args...> output;
    if constexpr (meta::has_reserve_v<Container<Scalar, Args...>>) {
      output.reserve(nr_points);
    }
    output.push_back(xx1);
    for (std::size_t ii = 0; ii < nr_points; ++ii) {
      output.push_back(output.back() + increment);
    }
    return output;
  } else {
    Container<Scalar, Args...> output(nr_points);
    for (std::size_t ii = 1; ii < output.size(); ++ii) {
      output[ii] = output[ii - 1] + increment;
    }
    return output;
  }
}

/** \return Linear partition of [\p xx1, \p xx2] with unit spacing. */
template <template <typename...> typename Container = std::vector,
          typename Scalar, typename... Args>
Container<Scalar, Args...> range(Scalar xx1, Scalar xx2) {
  return range(xx1, Scalar{1}, xx2);
}
} // namespace range

#endif /* Ranges_h */
