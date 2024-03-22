/**
 \file General/Operations.h
 \author Tomás Aquino
 \date 08/06/2019
 
\brief  Miscelaneous operations on containers.
 
 \note
 - Many methods assume containers with consistent sizes are passed in.
 - Many methods require random access and operator[].
 - Some operations involve casting; the latter are spelled out explicitly for clarity and to avoid warnings.
 - In most cases, the return value type is the type of the first container.
*/



#ifndef GENERAL_OPERATIONS_H
#define GENERAL_OPERATIONS_H

#include <cmath>
#include <functional>
#include <type_traits>
#include <vector>
#include "General/Useful.h"

namespace operation
{
  /** \brief Sum of elements. */
  template <typename Container>
  auto sum(Container const& input)
  {
    typename Container::value_type output{};
    for (auto const& val : input)
      output += val;
    return output;
  }

  /** \brief Product of elements. */
  template <typename Container>
  auto prod(Container const& input)
  {
    typename Container::value_type output{ 1. };
    for (auto const& val : input)
      output *= val;
    return output;
  }

  /** \brief Element-wise sum of scalar. */
  template <typename Container, typename Scalar, typename Container_out>
  void plus_scalar(Container const& input, Scalar cc, Container_out& output)
  {
    if constexpr (useful::has_plus<Container, Scalar>{})
      output = input + cc;
    else
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input[ii])+value_type(cc);
    }
  }

  /** \brief Element-wise sum of scalar. */
  template <typename Container, typename Scalar>
  auto plus_scalar(Container const& input, Scalar cc)
  {
    if constexpr (useful::has_plus<Container, Scalar>{})
      return input + cc;
    else
    {
      Container output(input.size());
      plus_scalar(input, cc, output);
      return output;
    }
  }

  /** \brief Element-wise sum of scalar. */
  template <typename Container, typename Scalar>
  void plus_scalar_InPlace(Container& input, Scalar cc)
  { plus_scalar(input, cc, input); }

  /** \brief Element-wise sum. */
  template <typename Container_1, typename Container_2, typename Container_out>
  void plus(Container_1 const& input_1, Container_2 const& input_2, Container_out& output)
  {
    if constexpr (useful::has_plus<Container_1, Container_2>{})
      output = input_1 + input_2;
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input_1[ii]) + value_type(input_2[ii]);
    }
    else
    {
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = input_1[ii] + input_2[ii];
    }
  }

  /** \brief Element-wise sum. */
  template <typename Container_1, typename Container_2>
  auto plus(Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_plus<Container_1, Container_2>{})
      return input_1 + input_2;
    else
    {
      Container_1 output(input_1.size());
      plus(input_1, input_2, output);
      return output;
    }
  }

  /** \brief Element-wise sum. */
  template <typename Container_1, typename Container_2>
  void plus_InPlace(Container_1& input_1, Container_2 const& input_2)
  { plus(input_1, input_2, input_1); }

  /** \brief Element-wise subtraction of scalar. */
  template <typename Container, typename Scalar, typename Container_out>
  void minus_scalar(Container const& input, Scalar cc, Container_out& output)
  {
    if constexpr (useful::has_minus<Container, Scalar>{})
      output = input - cc;
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input[ii]) - value_type(cc);
    }
    else
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = input[ii] - value_type(cc);
  }

  /** \brief Element-wise subtraction of scalar. */
  template <typename Container, typename Scalar>
  auto minus_scalar(Container const& input, Scalar cc)
  {
    if constexpr (useful::has_minus<Container, Scalar>{})
      return input - cc;
    else
    {
      Container output(input.size());
      minus_scalar(input, cc, output);
      return output;
    }
  }

  /** \brief Element-wise subtraction of scalar. */
  template <typename Container, typename Scalar>
  void minus_scalar_InPlace(Container& input, Scalar cc)
  { minus_scalar(input, cc, input); }
  
  /** \brief Element-wise subtraction from scalar. */
  template <typename Container, typename Scalar, typename Container_out>
  void scalar_minus(Scalar cc, Container const& input, Container_out& output)
  {
    if constexpr (useful::has_minus<Container, Scalar>{})
      output = cc - input;
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(cc)-value_type(input[ii]);
    }
    else
    {
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = cc-input[ii];
    }
  }

  /** \brief Element-wise subtraction from scalar. */
  template <typename Container, typename Scalar>
  auto scalar_minus(Scalar cc, Container const& input)
  {
    if constexpr (useful::has_minus<Container, Scalar>{})
      return cc - input;
    else
    {
      Container output(input.size());
      scalar_minus(cc, input, input);
      return output;
    }
  }

  /** \brief Element-wise subtraction from scalar. */
  template <typename Container, typename Scalar>
  void scalar_minus_InPlace(Scalar cc, Container& input)
  { scalar_minus(cc, input, input); }

  /** \brief Element-wise subtraction. */
  template <typename Container_1, typename Container_2, typename Container_out>
  void minus(Container_1 const& input_1, Container_2 const& input_2, Container_out& output)
  {
    if constexpr (useful::has_minus<Container_1, Container_2>{})
      output = input_1 - input_2;
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input_1[ii]) - value_type(input_2[ii]);
    }
    else
    {
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = input_1[ii] - input_2[ii];
    }
  }

  /** \brief Element-wise subtraction. */
  template <typename Container_1, typename Container_2>
  auto minus(Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_minus<Container_1, Container_2>{})
      return input_1 - input_2;
    else
    {
      Container_1 output = input_1;
      minus(input_1, input_2, output);
      return output;
    }
  }

  /** \brief Element-wise subtraction. */
  template <typename Container_1, typename Container_2>
  void minus_InPlace(Container_1& input_1, Container_2 const& input_2)
  { minus(input_1, input_2, input_1); }

  /** \brief Element-wise multiplication by scalar. */
  template <typename Container, typename Scalar, typename Container_out>
  void times_scalar(Scalar lambda, Container const& input, Container_out& output)
  {
    if constexpr (useful::has_multiplies<Scalar, Container>{})
    {
      if constexpr (std::is_convertible_v<Container_out,
                    std::decay_t<decltype(lambda*input)>>)
        output = lambda*input;
    }
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < input.size(); ++ii)
        output[ii] = value_type(lambda)*value_type(input[ii]);
    }
    else
    {
      for (std::size_t ii = 0; ii < input.size(); ++ii)
        output[ii] = lambda*input[ii];
    }
  }

  /** \brief Element-wise multiplication by scalar. */
  template <typename Container, typename Scalar>
  auto times_scalar(Scalar lambda, Container const& input)
  {
    if constexpr (useful::has_multiplies<Scalar, Container>{})
    {
      if constexpr (std::is_convertible_v<Container,
                    std::decay_t<decltype(lambda*input)>>)
        return lambda*input;
    }
    else
    {
      Container output(input.size());
      times_scalar(lambda, input, output);
      return output;
    }
  }

  /** \brief Element-wise multiplication by scalar. */
  template <typename Container, typename Scalar>
  void times_scalar_InPlace(Scalar lambda, Container& input)
  { times_scalar(lambda, input, input); }

  /** \brief Element-wise multiplication. */
  template <typename Container_1, typename Container_2, typename Container_out>
  void times(Container_1 const& input_1, Container_2 const& input_2, Container_out& output)
  {
    if constexpr (useful::has_multiplies<Container_1, Container_2>{})
    {
      if constexpr (std::is_convertible_v<Container_out,
                    std::decay_t<decltype(input_1*input_2)>>)
        output = input_1*input_2;
    }
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input_1[ii])*value_type(input_2[ii]);
    }
    else
    {
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = input_1[ii]*input_2[ii];
    }
  }

  /** \brief Element-wise multiplication. */
  template <typename Container_1, typename Container_2>
  auto times(Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_multiplies<Container_1, Container_2>{})
    {
      if constexpr (std::is_convertible_v<Container_1,
                    std::decay_t<decltype(input_1*input_2)>>)
        return input_1*input_2;
    }
    else
    {
      Container_1 output(input_1.size());
      times(input_1, input_2, output);

      return output;
    }
  }

  /** \brief Element-wise multiplication. */
  template <typename Container_1, typename Container_2>
  void times_InPlace(Container_1& input_1, Container_2 const& input_2)
  { times(input_1, input_2, input_1); }

  /** \brief Element-wise division by scalar. */
  template <typename Container, typename Scalar, typename Container_out>
  void div_scalar(Container const& input, Scalar lambda, Container_out& output)
  {
    if constexpr (useful::has_divides<Container, Scalar>{})
      output = input/lambda;
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < input.size(); ++ii)
        output[ii] = value_type(input[ii]/lambda);
    }
    else
    {
      for (std::size_t ii = 0; ii < input.size(); ++ii)
        output[ii] = input[ii]/lambda;
    }
  }

  /** \brief Element-wise division by scalar. */
  template <typename Container, typename Scalar>
  auto div_scalar(Container const& input, Scalar lambda)
  {
    if constexpr (useful::has_divides<Container, Scalar>{})
      return input/lambda;
    else
    {
      Container output(input.size());
      div_scalar(input, lambda, output);

      return output;
    }
  }
  
  /** \brief Element-wise division by scalar. */
  template <typename Container, typename Scalar>
  void div_scalar_InPlace(Container& input, Scalar lambda)
  { div_scalar(input, lambda, input); }

  /** \brief Element-wise division. */
  template <typename Container_1, typename Container_2, typename Container_out>
  void div(Container_1 const& input_1, Container_2 const& input_2, Container_out& output)
  {
    if constexpr (useful::has_divides<Container_1, Container_2>{})
      output = input_1/input_2;
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input_1[ii]/input_2[ii]);
    }
    else
    {
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = input_1[ii]/input_2[ii];
    }
  }

  /** \brief Element-wise division. */
  template <typename Container_1, typename Container_2>
  auto div(Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_divides<Container_1, Container_2>{})
      return input_1/input_2;
    else
    {
      Container_1 output(input_1.size());
      div(input_1, input_2, output);

      return output;
    }
  }

  /** \brief Element-wise division. */
  template <typename Container_1, typename Container_2>
  void div_InPlace(Container_1& input_1, Container_2 const& input_2)
  { div(input_1, input_2, input_1); }

  /** \brief \p lambda_1 *\p input_1 + \p lambda_2 * \p input2. */
  template <typename Container_1, typename Type_1,
  typename Container_2, typename Type_2,
  typename Container_out>
  void linearOp
  (Type_1 lambda_1, Container_1 const& input_1,
   Type_2 lambda_2, Container_2 const& input_2,
   Container_out& output)
  {
    if constexpr (useful::has_multiplies<Type_1, Container_1>{}
                  && useful::has_multiplies<Type_2, Container_2>{}
                  && useful::has_plus<Container_1, Container_2>{})
    {
      if constexpr (std::is_convertible_v<Container_out,
                    std::decay_t<decltype(lambda_1*input_1 + lambda_2*input_2)>>)
        output = lambda_1*input_1 + lambda_2*input_2;
    }
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for(size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(lambda_1)*input_1[ii]
        + value_type(lambda_2)*input_2[ii];
    }
    else
    {
      for(size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = lambda_1*input_1[ii]
        + lambda_2*input_2[ii];
    }
  }

  /** \brief \p lambda_1 * \p input_1 + \p lambda_2 * \p input2. */
  template <typename Container_1, typename Type_1,
  typename Container_2, typename Type_2>
  auto linearOp
  (Type_1 lambda_1, Container_1 const& input_1,
   Type_2 lambda_2, Container_2 const& input_2)
  {
    if constexpr (useful::has_multiplies<Type_1, Container_1>{}
                  && useful::has_multiplies<Type_2, Container_2>{}
                  && useful::has_plus<Container_1, Container_2>{})
    {
      if constexpr (std::is_convertible_v<Container_1,
                    std::decay_t<decltype(lambda_1*input_1 + lambda_2*input_2)>>)
        return lambda_1*input_1 + lambda_2*input_2;
    }
    else
    {
      Container_1 output(input_1.size());
      linearOp(lambda_1, input_1, lambda_2, input_2, output);
      return output;
    }
  }

  /** \brief \p lambda_1 * \p input_1 + \p lambda_2 * \p input2. */
  template <typename Container_1, typename Type_1,
  typename Container_2, typename Type_2>
  void linearOp_InPlace
  (Type_1 lambda_1, Container_1& input_1,
   Type_2 lambda_2, Container_2 const& input_2)
  { linearOp(lambda_1, input_1, lambda_2, input_2, input_1); }

  /** \brief \p lambda * \p input_1 + \p input2. */
  template <typename Container_1, typename Type,
  typename Container_2, typename Container_out>
  void linearOp
  (Type lambda, Container_1 const& input_1,
   Container_2 const& input_2, Container_out& output)
  {
    if constexpr (useful::has_multiplies<Type, Container_1>{}
                  && useful::has_plus<Container_1, Container_2>{})
    {
      if constexpr (std::is_convertible_v<Container_out,
                    std::decay_t<decltype(lambda*input_1 + input_2)>>)
        output = lambda*input_1 + input_2;
    }
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for(size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(lambda)*input_1[ii] +
          input_2[ii];
    }
    else
    {
      for(size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = lambda*input_1[ii] +
          input_2[ii];
    }
  }

  /** \brief \p lambda * \p input_1 + \p input2. */
  template <typename Container_1, typename Type, typename Container_2>
  auto linearOp
  (Type lambda, Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_multiplies<Type, Container_1>{}
                  && useful::has_plus<Container_1, Container_2>{})
    {
      if constexpr (std::is_convertible_v<Container_1,
                    std::decay_t<decltype(lambda*input_1 + input_2)>>)
        return lambda*input_1 + input_2;
    }
    else
    {
      Container_1 output(input_1.size());
      linearOp(lambda, input_1, input_2, output);

      return output;
    }
  }

  /** \brief \p lambda * \p input_1 + \p input2. */
  template <typename Container_1, typename Type,
  typename Container_2>
  void linearOp_InPlace
  (Type lambda, Container_1& input_1, Container_2 const& input_2)
  { linearOp(lambda, input_1, input_2, input_1); }

  /** \brief Element-wise square. */
  template <typename Container, typename Container_out>
  void square(Container const& input, Container_out& output)
  {
    if constexpr (useful::has_multiplies<Container, Container>{})
    {
      if constexpr (std::is_convertible_v<Container_out,
                    std::decay_t<decltype(input*input)>>)
        output = input*input;
    }
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = value_type(input[ii])*value_type(input[ii]);
    }
    else
    {
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = input[ii]*input[ii];
    }
  }

  /** \brief Element-wise square. */
  template <typename Container>
  auto square(Container const& input)
  {
    if constexpr (useful::has_multiplies<Container, Container>{})
    {
      if constexpr (std::is_convertible_v<Container,
                    std::decay_t<decltype(input*input)>>)
        input*input;
    }
    else
    {
      Container output(input.size());
      square(input,output);

      return output;
    }
  }

  /** \brief Element-wise square. */
  template <typename Container>
  void square_InPlace(Container& input)
  { square(input, input); }

  /** \brief Element-wise square root. */
  template <typename Container, typename Container_out>
  void sqrt(Container const& input, Container_out& output)
  {
    if constexpr (useful::can_call_sqrt_v<Container>)
      output = std::sqrt(input);
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = std::sqrt(value_type(input[ii]));
    }
    else
    {
      for (std::size_t ii = 0; ii < output.size(); ++ii)
        output[ii] = std::sqrt(input[ii]);
    }
  }

  /** \brief Element-wise square root. */
  template <typename Container>
  auto sqrt(Container const& input)
  {
    if constexpr (useful::can_call_sqrt_v<Container>)
      return std::sqrt(input);
    else
    {
      Container output(input.size());
      sqrt(input, output);

      return output;
    }
  }
  
  /** \brief Element-wise square root. */
  template <typename Container>
  void sqrt_InPlace(Container& input)
  { sqrt(input, input); }

  /** \brief Element-wise mean of two containers. */
  template <typename Container_1, typename Container_2,
  typename Container_out>
  void mean
  (Container_1 const& input_1, Container_2 const& input_2,
   Container_out& output)
  {
    if constexpr (useful::has_multiplies<Container_1, Container_2>{}
                  && useful::has_multiplies<double, Container_1>{})
    {
      if constexpr (std::is_convertible_v<Container_out,
                    std::decay_t<decltype(0.5*(input_1+input_2))>>)
        output = 0.5*(input_1+input_2);
    }
    else if constexpr (useful::has_value_type<Container_out>::value)
    {
      using value_type = typename Container_out::value_type;
      for (std::size_t ii = 0; ii <input_1.size(); ++ii)
        output[ii] = (value_type(input_1[ii]) + value_type(input_2[ii]))/2.;
    }
    else
    {
      for (std::size_t ii = 0; ii <input_1.size(); ++ii)
        output[ii] = (input_1[ii] + input_2[ii])/2.;
    }
  }

  /** \brief Element-wise mean of two containers. */
  template <typename Container_1, typename Container_2>
  auto mean(Container_1 const& input_1, Container_2 const& input_2)
  {
    if constexpr (useful::has_multiplies<Container_1, Container_2>{}
                  && useful::has_multiplies<double, Container_1>{})
    {
      if constexpr (std::is_convertible_v<Container_1,
                    std::decay_t<decltype(0.5*(input_1+input_2))>>)
        return 0.5*(input_1+input_2);
    }
    else
    {
      Container_1 output(input_1.size());
      mean(input_1, input_2, output);
      return output;
    }
  }

  /** \brief Element-wise mean of two containers. */
  template <typename Container>
  void mean_InPlace(Container& input_1, Container const& input_2)
  { mean(input_1, input_2, input_1); }

  /** \brief Euclidean norm squared. */
  template <typename Container>
  auto abs_sq(Container const& input)
  {
    if constexpr (useful::can_call_abs_v<Container>)
    {
      auto abs = std::abs(input);
      return abs*abs;
    }
    else if constexpr (useful::has_value_type<Container>::value)
    {
      typename Container::value_type abs_sq{ 0. };
      for (auto const& val : input)
        abs_sq += val*val;
      return abs_sq;
    }
    else
    {
      typename std::decay_t<decltype(input[0])> abs_sq{ 0. };
      for (auto const& val : input)
        abs_sq += val*val;
      return abs_sq;
    }
  }

  /** \brief Euclidean norm. */
  template <typename Container>
  auto abs(Container const& input)
  {
    if constexpr (useful::can_call_abs_v<Container>)
      return std::abs(input);
    else
      return std::sqrt(abs_sq(input));
  }

  /** \brief Like \c std::adjacent_difference but without the first element. */
  template
  <class InputIterator, class OutputIterator,
  class BinaryOperation>
  OutputIterator adjacent_difference
  (InputIterator first, InputIterator last,
   OutputIterator result,
   BinaryOperation binary_op = std::minus<typename InputIterator::value_type>{})
  {
    if (first != last)
    {
      InputIterator prev = first++;
      while (first != last)
      {
        InputIterator val = first++;
        *result++ = binary_op(*val, *prev);
        prev = val;
      }
    }
    return result;
  }
  
  /** \brief Averages of adjacent elements. */
  template <typename Container>
  void midpoints(Container const& input, Container& output)
  {
    for(size_t ii = 0; ii < input.size()-1; ++ii)
      output[ii] = (input[ii+1] + input[ii])/2.;
  }
  
  /** \brief Averages of adjacent elements. */
  template <typename Container>
  Container midpoints(Container const& input)
  {
    Container output(input.size()-1);
    midpoints(input, output);
    return output;
  }
  
  /** \brief Averages of pairs of elements. */
  template
  <typename... TArgs1, typename... TArgs2,
  template<typename...> typename Container,
  typename Value_type>
  void edge_midpoints
  (Container<std::pair<Value_type, Value_type>, TArgs1...> const& input,
   Container<Value_type, TArgs2...>& output)
  {
    for(size_t ii = 0; ii < input.size(); ++ii)
      output[ii] = (input[ii].first + input[ii].second)/2.;
  }
  
  /** \brief Averages of pairs of elements. */
  template
  <typename... TArgs,
  template<typename...> typename Container,
  typename Value_type>
  auto edge_midpoints
  (Container<std::pair<Value_type, Value_type>, TArgs...> const& input)
  {
    Container<Value_type, TArgs...> output(input.size());
    edge_midpoints(input, output);
    return output;
  }
  
  /** \brief Averages of pairs of elements. */
  auto edge_midpoints
  (std::vector<std::pair<double, double>> const& input)
  {
    std::vector<double> output(input.size());
    edge_midpoints(input, output);
    return output;
  }
  
  /** \brief Differences of pairs of elements. */
  template
  <template<typename> typename Container = std::vector,
  typename Value_type = double>
  void widths
  (Container<std::pair<Value_type, Value_type>> const& input,
   Container<Value_type>& output)
  {
    for(size_t ii = 0; ii < input.size(); ++ii)
      output[ii] = input.second[ii] - input.first[ii];
  }
  
  /** \brief Differences of pairs of elements. */
  template
  <template<typename> typename Container = std::vector,
  typename Value_type = double>
  Container<Value_type> widths
  (Container<std::pair<Value_type, Value_type>> const& input)
  {
    Container<Value_type> output(input.size());
    widths(input, output);
    return output;
  }
  
  /** \brief Differences of adjacent elements. */
  template <typename Container>
  void diff(Container const& input, Container& output)
  {
    for(size_t ii = 0; ii < input.size()-1; ++ii)
      output[ii] = input[ii+1] - input[ii];
  }
  
  /** \brief Differences of adjacent elements. */
  template <typename Container>
  Container diff(Container const& input)
  {
    Container output(input.size()-1);
    diff(input, output);
    return output;
  }

  /** \brief Dot product. */
  template <typename Container>
  auto dot(Container const& input_1, Container const& input_2)
  {
    typename std::decay_t<decltype(input_1[0])> result{};
    for (size_t ii = 0; ii < input_1.size(); ++ii)
      result += input_1[ii] * input_2[ii];
    return result;
  }
  
  /** \brief Dot product. */
  auto dot(double input_1, double input_2)
  {
    return input_1*input_2;
  }

  /** \brief Container of containers dotted into container, = sum_j (\p input_1)_{ij} (\p input_2)_j. */
  template <typename Container_outer, typename Container_inner>
  void dot
  (Container_outer const& input_1, Container_inner const& input_2, Container_inner& output)
  {
    std::size_t counter = 0;
    for (auto const& vec : input_1)
      output[counter++] = dot(vec, input_2);
  }

  /** \brief Container of containers dotted into container, = sum_j (\p input_1)_{ij} (\p input_2)_j. */
  template <typename Container_outer, typename Container_inner>
  auto dot(Container_outer const& input_1, Container_inner const& input_2)
  {
    if constexpr (useful::has_begin<decltype(input_1[0])>::value)
    {
      Container_inner output(input_2.size());
      dot(input_1, input_2, output);
      return output;
    }
    else
    {
      typename std::decay_t<decltype(input_1[0])> output{};
      for (size_t ii = 0; ii < input_1.size(); ++ii)
        output += input_1[ii] * input_2[ii];
      return output;
    }
  }
  
  /** \brief Rotate vector by angle. */
  template <typename Container>
  void rotate
  (Container const& input, double theta, Container& output)
  {
    double cos = std::cos(theta);
    double sin = std::sin(theta);
    output[0] = cos*input[0]-sin*input[1];
    output[1] = sin*input[0]+cos*input[1];
  }
  
  /** \brief Rotate vector by angle. */
  template <typename Container>
  auto rotate(Container const& vec, double theta)
  {
    Container output(vec.size());
    rotate(vec, theta, output);

    return output;
  }

  /** \brief Total number of elements in a container of containers. */
  template< template<class> class Container_outer, typename Container_inner >
  std::size_t nr_elements(Container_outer<Container_inner> const& vv)
  {
    std::size_t nr_elements = 0;
    for( auto const& cont : vv )
      nr_elements += cont.size();
    return nr_elements;
  }

  /** \brief Convolution sum. */
  template <typename Container>
  typename Container::value_type convolution
  (Container const& cont_1, Container const& cont_2,
   std::size_t idx_start, std::size_t idx_end)
  {
    typename std::decay_t<decltype(cont_1[0])> res{};
    for (std::size_t ii = idx_start; ii < idx_end; ++ii)
      res += cont_1[ii]*cont_2[idx_end-1-ii];
    return res;
  }

  /** \brief Convolution integral using trapezoidal rule. */
  template <typename Container>
  typename Container::value_type convolution_trap
  (Container const& cont_1, Container const& cont_2,
   std::size_t idx_start, std::size_t idx_end)
  {
    typename std::decay_t<decltype(cont_1[0])> res{};
    res += 0.5*(cont_1[idx_start]*cont_2[idx_end] + cont_1[idx_end]*cont_2[idx_start]);
    for (std::size_t ii = idx_start + 1; ii < idx_end; ++ii)
      res += cont_1[ii]*cont_2[idx_end-ii];
    return res;
  }

  /** \brief Hamming distance. */
  template <typename Container>
  auto hamming(Container const& vec, Container const& other_vec)
  {
    typename std::decay_t<decltype(vec[0])> hamming{};
    auto it = vec.begin();
    auto other_it = other_vec.begin();
    for (; it != vec.end(); it++, other_it++)
      hamming += *it > *other_it ? *it - *other_it : *other_it - *it;
    return hamming;
  }
  
  /** \brief Euclidean distance squared. */
  template <typename Container>
  double dist_sq(Container const& vec, Container const& other_vec)
  {
    return abs_sq(operation::minus(vec, other_vec));
  }
  
  /** \brief Euclidean distance. */
  template <typename Container>
  double dist(Container const& vec, Container const& other_vec)
  {
    return abs(operation::minus(vec, other_vec));
  }
  
  /** \brief Get component. */
  template <std::size_t dd, typename Container>
  auto project(Container const& container)
  { return container[dd]; }
  
  /** \brief Get component overload to get number itself from a \c double . */
  template <>
  auto project<0, double>(double const& val)
  { return val; }
  
  /** \brief Get component overload to get number itself from an \c int . */
  template <>
  auto project<0, int>(int const& val)
  { return val; }
  
  /** \brief Get component overload to get number itself from a \c size_t . */
  template <>
  auto project<0, std::size_t>(std::size_t const& val)
  { return val; }
  
  /** \return Factorial of \p nn. */
  std::size_t factorial(std::size_t nn)
  {
    if (nn == 0)
      return 1;
    return nn*factorial(nn-1);
  }

  /** \return \p nn (\p nn - 1)...(\p nn - \p mm). */
  std::size_t factorial_incomplete(std::size_t nn, std::size_t mm)
  {
    std::size_t result = 1.;
    ++nn;
    for (std::size_t ii = 0; ii < mm && nn --> 0; ++ii)
      result *= nn;
    return result;
  }
}


#endif /* GENERAL_OPERATIONS_H */
