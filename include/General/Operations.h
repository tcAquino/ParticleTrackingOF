/**
 \file General/Operations.h
 \author Tomás Aquino
 \date 08/06/2019

\brief  Miscelaneous operations on containers.

 \note
 - Many methods assume containers with consistent sizes are passed in.
 - Some methods require random access with operator[].
 - In most cases, the return value type is the type of the first container.
*/

#ifndef GENERAL_OPERATIONS_H
#define GENERAL_OPERATIONS_H

#include <cmath>
#include <stdexcept>
#include <functional>
#include <type_traits>
#include <vector>
#include "General/Meta.h"

/** \namespace op Uniform interface for operations on PODs and containers. */
namespace op {
/** \brief Sum of elements. */
template <typename Container> auto sum(Container const &input) {
  if constexpr (meta::has_begin_v<Container>) {
    std::decay_t<decltype(*input.begin())> output{0};
    for (auto const &val : input)
      output += val;
    return output;
  } else {
    std::decay_t<decltype(Container()[0])> output{};
    for (std::size_t ii = 0; ii < input.size(); ++ii)
      output += input[ii];
    return output;
  }
}

/** \brief Sum of elements. */
inline auto sum(double input) { return input; }

/** \brief Sum of elements. */
inline auto sum(int input) { return input; }

/** \brief Sum of elements. */
inline auto sum(std::size_t input) { return input; }

/** \brief Product of elements. */
template <typename Container> auto prod(Container const &input) {
  if constexpr (meta::has_begin_v<Container>) {
    std::decay_t<decltype(*input.begin())> output{1};
    for (auto const &val : input)
      output *= val;
    return output;
  } else {
    std::decay_t<decltype(Container()[0])> output{1};
    for (std::size_t ii = 0; ii < input.size(); ++ii)
      output *= input[ii];
    return output;
  }
}

/** \brief Product of elements. */
inline auto prod(double input) { return input; }

/** \brief Product of elements. */
inline auto prod(int input) { return input; }

/** \brief Product of elements. */
inline auto prod(std::size_t input) { return input; }

/** \brief Element-wise sum of scalar. */
template <typename Container, typename Scalar, typename Container_out>
void plus_scalar(Container const &input, Scalar cc, Container_out &output) {
  if constexpr (meta::is_convertible_from_plus_v<Container, Container, Scalar>)
    output = input + cc;
  else
    for (std::size_t ii = 0; ii < output.size(); ++ii)
      output[ii] = input[ii] + cc;
}

/** \brief Element-wise sum of scalar. */
template <typename Container, typename Scalar>
auto plus_scalar(Container const &input, Scalar cc) {
  if constexpr (meta::is_convertible_from_plus_v<Container, Container, Scalar>)
    return input + cc;
  if constexpr (meta::is_constructible_from_size_v<Container>) {
    Container output(input.size());
    plus_scalar(input, cc, output);
    return output;
  } else {
    Container output{};
    plus_scalar(input, cc, output);
    return output;
  }
}

/** \brief Element-wise sum of scalar. */
template <typename Container, typename Scalar>
void plus_scalar_inplace(Container &input, Scalar cc) {
  plus_scalar(input, cc, input);
}

/** \brief Element-wise sum. */
template <typename Container_1, typename Container_2, typename Container_out>
void plus(Container_1 const &input_1, Container_2 const &input_2,
          Container_out &output) {
  if constexpr (meta::is_convertible_from_plus_v<Container_out, Container_1,
                                                 Container_2>)
    output = input_1 + input_2;
  else {
    for (std::size_t ii = 0; ii < output.size(); ++ii)
      output[ii] = input_1[ii] + input_2[ii];
  }
}

/** \brief Element-wise sum. */
template <typename Container_1, typename Container_2>
auto plus(Container_1 const &input_1, Container_2 const &input_2) {
  if constexpr (meta::is_convertible_from_plus_v<Container_1, Container_1,
                                                 Container_2>)
    return input_1 + input_2;
  else {
    if constexpr (meta::is_constructible_from_size_v<Container_1>) {
      Container_1 output(input_1.size());
      plus(input_1, input_2, output);
      return output;
    } else {
      Container_1 output{};
      plus(input_1, input_2, output);
      return output;
    }
  }
}

/** \brief Element-wise sum. */
template <typename Container_1, typename Container_2>
void plus_inplace(Container_1 &input_1, Container_2 const &input_2) {
  plus(input_1, input_2, input_1);
}

/** \brief Element-wise subtraction of scalar. */
template <typename Container, typename Scalar, typename Container_out>
void minus_scalar(Container const &input, Scalar cc, Container_out &output) {
  if constexpr (meta::is_convertible_from_minus_v<Container_out, Container,
                                                  Scalar>)
    output = input - cc;
  else
    for (std::size_t ii = 0; ii < output.size(); ++ii)
      output[ii] = input[ii] - cc;
}

/** \brief Element-wise subtraction of scalar. */
template <typename Container, typename Scalar>
auto minus_scalar(Container const &input, Scalar cc) {
  if constexpr (meta::is_convertible_from_minus_v<Container, Container, Scalar>)
    return input - cc;
  else {
    if constexpr (meta::is_constructible_from_size_v<Container>) {
      Container output(input.size());
      minus_scalar(input, cc, output);
      return output;
    } else {
      Container output{};
      minus_scalar(input, cc, output);
      return output;
    }
  }
}

/** \brief Element-wise subtraction of scalar. */
template <typename Container, typename Scalar>
void minus_scalar_inplace(Container &input, Scalar cc) {
  minus_scalar(input, cc, input);
}

/** \brief Element-wise subtraction from scalar. */
template <typename Container, typename Scalar, typename Container_out>
void scalar_minus(Scalar cc, Container const &input, Container_out &output) {
  if constexpr (meta::is_convertible_from_minus_v<Container_out, Scalar,
                                                  Container>)
    output = cc - input;
  else {
    for (std::size_t ii = 0; ii < output.size(); ++ii)
      output[ii] = cc - input[ii];
  }
}

/** \brief Element-wise subtraction from scalar. */
template <typename Container, typename Scalar>
auto scalar_minus(Scalar cc, Container const &input) {
  if constexpr (meta::is_convertible_from_minus_v<Container, Scalar, Container>)
    return cc - input;
  else {
    if constexpr (meta::is_constructible_from_size_v<Container>) {
      Container output(input.size());
      scalar_minus(cc, input, input);
      return output;
    } else {
      Container output{};
      scalar_minus(cc, input, input);
      return output;
    }
  }
}

/** \brief Element-wise subtraction from scalar. */
template <typename Container, typename Scalar>
void scalar_minus_inplace(Scalar cc, Container &input) {
  scalar_minus(cc, input, input);
}

/** \brief Element-wise subtraction. */
template <typename Container_1, typename Container_2, typename Container_out>
void minus(Container_1 const &input_1, Container_2 const &input_2,
           Container_out &output) {
  if constexpr (meta::is_convertible_from_minus_v<Container_out, Container_1,
                                                  Container_2>)
    output = input_1 - input_2;
  else {
    for (std::size_t ii = 0; ii < output.size(); ++ii)
      output[ii] = input_1[ii] - input_2[ii];
  }
}

/** \brief Element-wise subtraction. */
template <typename Container_1, typename Container_2>
auto minus(Container_1 const &input_1, Container_2 const &input_2) {
  if constexpr (meta::is_convertible_from_minus_v<Container_1, Container_1,
                                                  Container_2>)
    return input_1 - input_2;
  else {
    Container_1 output = input_1;
    minus(input_1, input_2, output);
    return output;
  }
}

/** \brief Element-wise subtraction. */
template <typename Container_1, typename Container_2>
void minus_inplace(Container_1 &input_1, Container_2 const &input_2) {
  minus(input_1, input_2, input_1);
}

/** \brief Element-wise multiplication by scalar. */
template <typename Container, typename Scalar, typename Container_out>
void times_scalar(Scalar lambda, Container const &input,
                  Container_out &output) {
  if constexpr (meta::is_convertible_from_multiplies_v<Container_out, Scalar,
                                                       Container>)
    output = lambda * input;
  else
    for (std::size_t ii = 0; ii < input.size(); ++ii)
      output[ii] = lambda * input[ii];
}

/** \brief Element-wise multiplication by scalar. */
template <typename Container, typename Scalar>
auto times_scalar(Scalar lambda, Container const &input) {
  if constexpr (meta::is_convertible_from_multiplies_v<Container, Scalar,
                                                       Container>)
    return lambda * input;
  if constexpr (meta::is_constructible_from_size_v<Container>) {
    Container output(input.size());
    times_scalar(lambda, input, output);
    return output;
  } else {
    Container output{};
    times_scalar(lambda, input, output);
    return output;
  }
}

/** \brief Element-wise multiplication by scalar. */
template <typename Container, typename Scalar>
void times_scalar_inplace(Scalar lambda, Container &input) {
  times_scalar(lambda, input, input);
}

/** \brief Element-wise multiplication. */
template <typename Container_1, typename Container_2, typename Container_out>
void times(Container_1 const &input_1, Container_2 const &input_2,
           Container_out &output) {
  if constexpr (meta::is_convertible_from_multiplies_v<
                    Container_out, Container_1, Container_2>)
    output = input_1 * input_2;
  else
    for (std::size_t ii = 0; ii < output.size(); ++ii)
      output[ii] = input_1[ii] * input_2[ii];
}

/** \brief Element-wise multiplication. */
template <typename Container_1, typename Container_2>
auto times(Container_1 const &input_1, Container_2 const &input_2) {
  if constexpr (meta::is_convertible_from_multiplies_v<Container_1, Container_1,
                                                       Container_2>)
    return input_1 * input_2;
  else {
    if constexpr (meta::is_constructible_from_size_v<Container_1>) {
      Container_1 output(input_1.size());
      times(input_1, input_2, output);
      return output;
    } else {
      Container_1 output{};
      times(input_1, input_2, output);
      return output;
    }
  }
}

/** \brief Element-wise multiplication. */
template <typename Container_1, typename Container_2>
void times_inplace(Container_1 &input_1, Container_2 const &input_2) {
  times(input_1, input_2, input_1);
}

/** \brief Element-wise division by scalar. */
template <typename Container, typename Scalar, typename Container_out>
void div_scalar(Container const &input, Scalar lambda, Container_out &output) {
  if constexpr (meta::is_convertible_from_divides_v<Container_out, Container,
                                                    Scalar>)
    output = input / lambda;
  else
    for (std::size_t ii = 0; ii < input.size(); ++ii)
      output[ii] = input[ii] / lambda;
}

/** \brief Element-wise division by scalar. */
template <typename Container, typename Scalar>
auto div_scalar(Container const &input, Scalar lambda) {
  if constexpr (meta::is_convertible_from_divides_v<Container, Container,
                                                    Scalar>)
    return input / lambda;
  else {
    if constexpr (meta::is_constructible_from_size_v<Container>) {
      Container output(input.size());
      div_scalar(input, lambda, output);
      return output;
    } else {
      Container output{};
      div_scalar(input, lambda, output);
      return output;
    }
  }
}

/** \brief Element-wise division by scalar. */
template <typename Container, typename Scalar>
void div_scalar_inplace(Container &input, Scalar lambda) {
  div_scalar(input, lambda, input);
}

/** \brief Element-wise division. */
template <typename Container_1, typename Container_2, typename Container_out>
void div(Container_1 const &input_1, Container_2 const &input_2,
         Container_out &output) {
  if constexpr (meta::is_convertible_from_divides_v<Container_out, Container_1,
                                                    Container_2>)
    output = input_1 / input_2;
  else {
    for (std::size_t ii = 0; ii < output.size(); ++ii)
      output[ii] = input_1[ii] / input_2[ii];
  }
}

/** \brief Element-wise division. */
template <typename Container_1, typename Container_2>
auto div(Container_1 const &input_1, Container_2 const &input_2) {
  if constexpr (meta::is_convertible_from_divides_v<Container_1, Container_1,
                                                    Container_2>)
    return input_1 / input_2;
  else {
    if constexpr (meta::is_constructible_from_size_v<Container_1>) {
      Container_1 output(input_1.size());
      div(input_1, input_2, output);
      return output;
    } else {
      Container_1 output{};
      div(input_1, input_2, output);
      return output;
    }
  }
}

/** \brief Element-wise division. */
template <typename Container_1, typename Container_2>
void div_inplace(Container_1 &input_1, Container_2 const &input_2) {
  div(input_1, input_2, input_1);
}

/** \brief \c lambda_1*input_1+lambda_2*input2. */
template <typename Container_1, typename Scalar_1, typename Container_2,
          typename Scalar_2, typename Container_out>
void linearop(Scalar_1 lambda_1, Container_1 const &input_1, Scalar_2 lambda_2,
              Container_2 const &input_2, Container_out &output) {
  if constexpr (meta::is_convertible_from_multiplies_v<Container_1, Scalar_1,
                                                       Container_1> &&
                meta::is_convertible_from_multiplies_v<Container_2, Scalar_2,
                                                       Container_2> &&
                meta::is_convertible_from_plus_v<Container_out, Container_1,
                                                 Container_2>)
    output = lambda_1 * input_1 + lambda_2 * input_2;
  else
    for (size_t ii = 0; ii < output.size(); ++ii)
      output[ii] = lambda_1 * input_1[ii] + lambda_2 * input_2[ii];
}

/** \brief \c lambda_1*input_1+lambda_2*input2. */
template <typename Container_1, typename Scalar_1, typename Container_2,
          typename Scalar_2>
auto linearop(Scalar_1 lambda_1, Container_1 const &input_1, Scalar_2 lambda_2,
              Container_2 const &input_2) {
  if constexpr (meta::is_convertible_from_multiplies_v<Container_1, Scalar_1,
                                                       Container_1> &&
                meta::is_convertible_from_multiplies_v<Container_2, Scalar_2,
                                                       Container_2> &&
                meta::is_convertible_from_plus_v<Container_1, Container_1,
                                                 Container_2>)
    return lambda_1 * input_1 + lambda_2 * input_2;
  else {
    if constexpr (meta::is_constructible_from_size_v<Container_1>) {
      Container_1 output(input_1.size());
      linearop(lambda_1, input_1, lambda_2, input_2, output);
      return output;
    } else {
      Container_1 output{};
      linearop(lambda_1, input_1, lambda_2, input_2, output);
      return output;
    }
  }
}

/** \brief \c lambda_1*input_1+lambda_2*input2. */
template <typename Container_1, typename Scalar_1, typename Container_2,
          typename Scalar_2>
void linearop_inplace(Scalar_1 lambda_1, Container_1 &input_1,
                      Scalar_2 lambda_2, Container_2 const &input_2) {
  linearop(lambda_1, input_1, lambda_2, input_2, input_1);
}

/** \brief \c lambda*input_1+input2. */
template <typename Container_1, typename Scalar, typename Container_2,
          typename Container_out>
void linearop(Scalar lambda, Container_1 const &input_1,
              Container_2 const &input_2, Container_out &output) {
  if constexpr (meta::is_convertible_from_multiplies_v<Container_1, Scalar,
                                                       Container_1> &&
                meta::is_convertible_from_plus_v<Container_out, Container_1,
                                                 Container_2>)
    output = lambda * input_1 + input_2;
  else
    for (size_t ii = 0; ii < output.size(); ++ii)
      output[ii] = lambda * input_1[ii] + input_2[ii];
}

/** \brief \c lambda*input_1+input2. */
template <typename Container_1, typename Type, typename Container_2>
auto linearop(Type lambda, Container_1 const &input_1,
              Container_2 const &input_2) {
  if constexpr (meta::has_multiplies_v<Type, Container_1> &&
                meta::has_plus_v<Container_1, Container_2>)
    return lambda * input_1 + input_2;
  else {
    if constexpr (meta::is_constructible_from_size_v<Container_1>) {
      Container_1 output(input_1.size());
      linearop(lambda, input_1, input_2, output);
      return output;
    } else {
      Container_1 output{};
      linearop(lambda, input_1, input_2, output);
      return output;
    }
  }
}

/** \brief \c lambda*input_1+input2. */
template <typename Container_1, typename Type, typename Container_2>
void linearop_inplace(Type lambda, Container_1 &input_1,
                      Container_2 const &input_2) {
  linearop(lambda, input_1, input_2, input_1);
}

/** \brief Element-wise square. */
template <typename Container, typename Container_out>
void square(Container const &input, Container_out &output) {
  if constexpr (meta::is_convertible_from_multiplies_v<Container, Container,
                                                       Container>)
    output = input * input;
  else
    for (std::size_t ii = 0; ii < output.size(); ++ii)
      output[ii] = input[ii] * input[ii];
}

/** \brief Element-wise square. */
template <typename Container> auto square(Container const &input) {
  if constexpr (meta::is_convertible_from_multiplies_v<Container, Container,
                                                       Container>)
    return input * input;
  else {
    if constexpr (meta::is_constructible_from_size_v<Container>) {
      Container output(input.size());
      square(input, output);
      return output;
    } else {
      Container output{};
      square(input, output);
      return output;
    }
  }
}

/** \brief Element-wise square. */
template <typename Container> void square_inplace(Container &input) {
  square(input, input);
}

/** \brief Element-wise square root. */
template <typename Container, typename Container_out>
void sqrt(Container const &input, Container_out &output) {
  if constexpr (meta::is_convertible_from_sqrt_v<Container_out, Container>)
    output = std::sqrt(input);
  else
    for (std::size_t ii = 0; ii < output.size(); ++ii)
      output[ii] = std::sqrt(input[ii]);
}

/** \brief Element-wise square root. */
template <typename Container> auto sqrt(Container const &input) {
  if constexpr (meta::is_convertible_from_sqrt_v<Container, Container>)
    return std::sqrt(input);
  else {
    if constexpr (meta::is_constructible_from_size_v<Container>) {
      Container output(input.size());
      sqrt(input, output);
      return output;
    } else {
      Container output{};
      sqrt(input, output);
      return output;
    }
  }
}

/** \brief Element-wise square root. */
template <typename Container> void sqrt_inplace(Container &input) {
  sqrt(input, input);
}

/** \brief Element-wise mean of two containers. */
template <typename Container_1, typename Container_2, typename Container_out>
void mean(Container_1 const &input_1, Container_2 const &input_2,
          Container_out &output) {
  if constexpr (meta::is_convertible_from_multiplies_v<
                    Container_out, Container_1, Container_2> &&
                meta::is_convertible_from_multiplies_v<Container_out, double, Container_out>)
    output = 0.5 * (input_1 + input_2);
  else {
    for (std::size_t ii = 0; ii < input_1.size(); ++ii)
      output[ii] = 0.5 * (input_1[ii] + input_2[ii]);
  }
}

/** \brief Element-wise mean of two containers. */
template <typename Container_1, typename Container_2>
auto mean(Container_1 const &input_1, Container_2 const &input_2) {
  if constexpr (meta::is_convertible_from_multiplies_v<Container_1, Container_1,
                                                       Container_2> &&
                meta::is_convertible_from_multiplies_v<Container_1, double, Container_1>)
    return 0.5 * (input_1 + input_2);
  else {
    if constexpr (meta::is_constructible_from_size_v<Container_1>) {
      Container_1 output(input_1.size());
      mean(input_1, input_2, output);
      return output;
    } else {
      Container_1 output{};
      mean(input_1, input_2, output);
      return output;
    }
  }
}
/** \brief Element-wise mean of two containers. */
template <typename Container>
void mean_inplace(Container &input_1, Container const &input_2) {
  mean(input_1, input_2, input_1);
}

/** \brief Euclidean norm squared. */
template <typename Container> auto abs_sq(Container const &input) {
  if constexpr (meta::is_convertible_from_abs_v<double, Container>) {
    auto abs = std::abs(input);
    return abs * abs;
  } else if constexpr (meta::has_begin_v<Container>) {
    std::decay_t<decltype(*input.begin())> abs_sq{0};
    for (auto const &val : input)
      abs_sq += val * val;
    return abs_sq;
  } else {
    std::decay_t<decltype(Container()[0])> abs_sq{0};
    for (std::size_t ii = 0; ii < input.size(); ++ii)
      abs_sq += input[ii] * input[ii];
    return abs_sq;
  }
}

/** \brief Euclidean norm. */
template <typename Container> auto abs(Container const &input) {
  if constexpr (meta::is_convertible_from_abs_v<double, Container>)
    return std::abs(input);
  else
    return std::sqrt(abs_sq(input));
}

/** \brief Like \c std::adjacent_difference but without the first element. */
template <class InputIterator, class OutputIterator, class BinaryOperation>
OutputIterator
adjacent_difference(InputIterator first, InputIterator last,
                    OutputIterator result,
                    BinaryOperation binary_op =
                        std::minus<typename InputIterator::value_type>{}) {
  if (first != last) {
    InputIterator prev = first++;
    while (first != last) {
      InputIterator val = first++;
      *result++ = binary_op(*val, *prev);
      prev = val;
    }
  }
  return result;
}

/** \brief Averages of adjacent elements. */
template <typename Container>
void midpoints(Container const &input, Container &output) {
  for (size_t ii = 0; ii < input.size() - 1; ++ii)
    output[ii] = (input[ii + 1] + input[ii]) / 2.;
}

/** \brief Averages of adjacent elements. */
template <typename Container> Container midpoints(Container const &input) {
  if constexpr (meta::is_constructible_from_size_v<Container>) {
    Container output(input.size() - 1);
    midpoints(input, output);
    return output;
  } else {
    Container output{};
    midpoints(input, output);
    return output;
  }
}

/** \brief Differences of adjacent elements. */
template <typename Container>
void diff(Container const &input, Container &output) {
  for (size_t ii = 0; ii < input.size() - 1; ++ii)
    output[ii] = input[ii + 1] - input[ii];
}

/** \brief Differences of adjacent elements. */
template <typename Container> Container diff(Container const &input) {
  Container output(input.size() - 1);
  diff(input, output);
  return output;
}

/** \brief Compute widths of bins having given \p values as midpoints and given
 \p minimum (left edge). \details Widths are computed sequentially as twice the
 distance between the next midpoint and the current left edge. */
template <template <typename...> typename Container = std::vector,
          typename Scalar = double, typename... Args>
auto bin_widths(Container<Scalar, Args...> const &values, Scalar minimum) {
  Container<Scalar, Args...> widths;
  if constexpr (meta::has_reserve_v<Container<Scalar, Args...>>)
    widths.reserve(values.size());
  else if constexpr (meta::is_constructible_from_size_v<
                         Container<Scalar, Args...>>)
    widths = Container<Scalar, Args...>(values.size());
  double left_edge = minimum;
  if constexpr (meta::has_push_back_v<Container<Scalar, Args...>>)
    for (auto const &val : values) {
      widths.push_back(2. * (val - left_edge));
      if (widths.back() < 0.)
        throw std::runtime_error{
            "Bad minimum value, could not determine bin widths"};
      left_edge += widths.back();
    }
  else
    for (std::size_t ii = 0; ii < values.size(); ++ii) {
      widths[ii] = 2. * (values[ii] - left_edge);
      if (widths[ii] < 0.)
        throw std::runtime_error{
            "Bad minimum value, could not determine bin widths"};
      left_edge += widths[ii];
    }

  return widths;
}

/** \brief Compute widths of bins having given \p values as midpoints.
 \details Leftmost bin edge is the minimum value minus half the width of the
 first bin.

 Widths are computed sequentially as twice the distance between the next
 midpoint and the current left edge. */
template <typename Container = std::vector<double>>
auto bin_widths(Container const &values) {
  if (values.size() < 2)
    throw std::runtime_error{
        "Less than two values, could not determine bin widths"};
  double midpoint = (values[1] + values[0]) / 2.;
  return get_bin_widths(values, values[0] - (midpoint - values[0]));
}

/** \brief Dot product. */
template <typename Container>
auto dot(Container const &input_1, Container const &input_2) {
  typename std::decay_t<decltype(input_1[0])> result{};
  for (size_t ii = 0; ii < input_1.size(); ++ii)
    result += input_1[ii] * input_2[ii];
  return result;
}

/** \brief Dot product. */
inline auto dot(double input_1, double input_2) { return input_1 * input_2; }

/** \brief Dot product. */
inline auto dot(int input_1, int input_2) { return input_1 * input_2; }

/** \brief Dot product. */
inline auto dot(std::size_t input_1, std::size_t input_2) { return input_1 * input_2; }

/** \brief Container of containers dotted into container, \c
 * =sum_j(\input_1)_{ij}(input_2)_j. */
template <typename Container_outer, typename Container_inner>
void dot(Container_outer const &input_1, Container_inner const &input_2,
         Container_inner &output) {
  std::size_t counter = 0;
  for (auto const &vec : input_1)
    output[counter++] = dot(vec, input_2);
}

/** \brief Container of containers dotted into container, \c
 * =sum_j(input_1)_{ij}(input_2)_j. */
template <typename Container_outer, typename Container_inner>
auto dot(Container_outer const &input_1, Container_inner const &input_2) {
  if constexpr (meta::has_begin_v<decltype(input_1[0])>) {
    Container_inner output(input_2.size());
    dot(input_1, input_2, output);
    return output;
  } else {
    typename std::decay_t<decltype(input_1[0])> output{};
    for (size_t ii = 0; ii < input_1.size(); ++ii)
      output += input_1[ii] * input_2[ii];
    return output;
  }
}

/** \brief Rotate vector by angle. */
template <typename Container>
void rotate(Container const &input, double theta, Container &output) {
  double cos = std::cos(theta);
  double sin = std::sin(theta);
  output[0] = cos * input[0] - sin * input[1];
  output[1] = sin * input[0] + cos * input[1];
}

/** \brief Rotate vector by angle. */
template <typename Container> auto rotate(Container const &vec, double theta) {
  Container output(vec.size());
  rotate(vec, theta, output);

  return output;
}

/** \brief Total number of elements in a container of containers. */
template <template <class> class Container_outer, typename Container_inner>
std::size_t nr_elements(Container_outer<Container_inner> const &vv) {
  std::size_t nr_elements = 0;
  if constexpr (meta::has_begin_v<Container_outer<Container_inner>>)
    for (auto const &cont : vv)
      nr_elements += cont.size();
  else
    for (std::size_t ii = 0; ii < vv.size(); ++vv)
      nr_elements += vv[ii].size();
  return nr_elements;
}

/** \brief Convolution sum. */
template <typename Container>
auto convolution(Container const &cont_1, Container const &cont_2,
                 std::size_t idx_start, std::size_t idx_end) {
  typename std::decay_t<decltype(cont_1[0])> res{};
  for (std::size_t ii = idx_start; ii < idx_end; ++ii)
    res += cont_1[ii] * cont_2[idx_end - 1 - ii];
  return res;
}

/** \brief Convolution integral using trapezoidal rule. */
template <typename Container>
auto convolution_trap(Container const &cont_1, Container const &cont_2,
                      std::size_t idx_start, std::size_t idx_end) {
  typename std::decay_t<decltype(cont_1[0])> res{};
  res += 0.5 * (cont_1[idx_start] * cont_2[idx_end] +
                cont_1[idx_end] * cont_2[idx_start]);
  for (std::size_t ii = idx_start + 1; ii < idx_end; ++ii)
    res += cont_1[ii] * cont_2[idx_end - ii];
  return res;
}

/** \brief Hamming distance. */
template <typename Container>
auto hamming(Container const &vec, Container const &other_vec) {
  typename std::decay_t<decltype(vec[0])> hamming{};
  auto it = vec.begin();
  auto other_it = other_vec.begin();
  for (; it != vec.end(); it++, other_it++)
    hamming += *it > *other_it ? *it - *other_it : *other_it - *it;
  return hamming;
}

/** \brief Euclidean distance squared. */
template <typename Container>
double dist_sq(Container const &vec, Container const &other_vec) {
  return abs_sq(minus(vec, other_vec));
}

/** \brief Euclidean distance. */
template <typename Container>
double dist(Container const &vec, Container const &other_vec) {
  return abs(minus(vec, other_vec));
}

/** \brief Get component. */
template <std::size_t dd, typename Container>
auto project(Container const &container) {
  return container[dd];
}

/** \brief Get component overload to get number itself from a \c double . */
template <> inline auto project<0, double>(double const &val) { return val; }

/** \brief Get component overload to get number itself from an \c int . */
template <> inline auto project<0, int>(int const &val) { return val; }

/** \brief Get component overload to get number itself from a \c size_t . */
template <> inline auto project<0, std::size_t>(std::size_t const &val) { return val; }

/** \return Factorial of \c nn. */
inline std::size_t factorial(std::size_t nn) {
  if (nn == 0)
    return 1;
  return nn * factorial(nn - 1);
}

/** \return \c nn(nn-1)...(nn-mm). */
inline std::size_t factorial_incomplete(std::size_t nn, std::size_t mm) {
  std::size_t result = 1.;
  ++nn;
  for (std::size_t ii = 0; ii < mm && nn-- > 0; ++ii)
    result *= nn;
  return result;
}

/** \return Sign of \c val. */
template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }
} // namespace op

#endif /* GENERAL_OPERATIONS_H */
