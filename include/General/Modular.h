/**
   \file General/Modular.h
   \author Tomas Aquino
   \date 08/06/2019
   \brief Modular arithmetic, basis change, and similar operations
*/

#ifndef GENERAL_MODULAR_H
#define GENERAL_MODULAR_H

#include <cmath>
#include <cstdlib>
#include <vector>

namespace modular {
/** \return \p number + 1 mod \p mod. */
inline std::size_t circPlus(std::size_t number, std::size_t mod) {
  return ++number == mod ? 0 : number;
}

/** \return \p number - 1  mod \p mod. */
inline std::size_t circMinus(std::size_t number, std::size_t mod) {
  return number == 0 ? mod - 1 : --number;
}

/** \return \p number + \p cc mod \p mod. */
inline std::size_t circPlus(std::size_t number, std::size_t cc,
                            std::size_t mod) {
  return (number + cc) % mod;
}

/** \return \p number - \p cc mod \p mod. */
inline std::size_t circMinus(std::size_t number, std::size_t cc,
                             std::size_t mod) {
  return number < cc ? (mod - (cc - number) % mod) % mod : number - cc;
}

/**
   \param to_convert Number given as a container of integer digits.
   \param basis Basis of the number given.
   \return Number converted to base ten.
   \note
   - Units digit in <tt>to_convert[0]</tt>.
*/
template <typename Container>
auto convert(Container const &to_convert, std::size_t basis) {
  typename Container::value_type converted = 0;
  typename Container::value_type power = 1;
  for (std::size_t pos = to_convert.size(); pos-- > 0;) {
    converted += to_convert[pos] * power;
    power *= basis;
  }
  return converted;
}

/**
   \param to_convert Base-10 integer.
   \param basis Basis to convert to.
   \param converted Result of conversion, given as a container of integer
   digits.
   \note Units digit in <tt>converted[0]</tt>.
*/
template <typename Container, typename Integer_Type>
void convert(Integer_Type to_convert, Integer_Type basis,
             Container &converted) {
  Integer_Type to_convert_old;
  for (std::size_t pos = converted.size(); pos-- > 0;) {
    to_convert_old = to_convert;
    to_convert /= basis;
    converted[pos] = to_convert_old - to_convert * basis;
  }
}

/**
   \return Base \p basis representation of base-10 integer \p to_convert as a
   string of digits.
   \note Units digit in <tt>to_convert[0]</tt>.
*/
template <typename Integer_Type>
std::vector<Integer_Type> convert(Integer_Type to_convert, Integer_Type basis) {
  std::vector<Integer_Type> converted(std::size_t(std::log2(to_convert) + 1));
  convert(to_convert, basis, converted);
  return converted;
}

/**
   \return Base \p basis representation of base-10 integer \p to_convert as a
   std::ring of digits.
   \details The total number of digits is \p nr_digits, which must be at least
   the number of digits needed to represent the number. If bigger, the result is
   padded with zeros on the left.
   \note Units digit in <tt>to_convert[0]</tt>.
*/
template <typename Integer_Type>
std::vector<std::size_t> convert(Integer_Type to_convert, Integer_Type basis,
                                 std::size_t nr_digits) {
  std::vector<Integer_Type> converted(nr_digits);
  Integer_Type to_convert_old;
  for (std::size_t pos = std::size_t(std::log2(to_convert) + 1); pos-- > 0;) {
    to_convert_old = to_convert;
    to_convert /= basis;
    converted[pos] = to_convert_old - to_convert * basis;
  }
  return converted;
}

/**
   \brief Given 1d index \p idx, convert to position <tt>(x_idx, y_idx,
   ...)</tt> in a regular grid; output in \p position.
   \details \p idx goes through every x for each value of y, every y for each
   value of z, etc, in ascending order.
   \note: \p position container must be passed with the correct size.
*/
template <typename Container, typename IntegerType = std::size_t>
void indexToPosition(IntegerType idx, Container const &nr_points,
                     Container &position) {
  for (size_t dd = 0; dd < nr_points.size(); ++dd) {
    position[dd] = idx % nr_points[dd];
    idx /= nr_points[dd];
  }
}

/**
   \return Position <tt>(x_idx, y_idx, ...)</tt> in a regular grid, given 1d
   index \p idx. \details \p idx goes through every x for each value of y, every
   y for each value of z, etc, in ascending order.
*/
template <typename Container, typename IntegerType = std::size_t>
auto indexToPosition(IntegerType idx, Container const &nr_points) {
  std::vector<IntegerType> position;
  position.reserve(nr_points.size());
  for (size_t dd = 0; dd < nr_points.size(); ++dd) {
    position.push_back(idx % nr_points[dd]);
    idx /= nr_points[dd];
  }
  return position;
}

/**
   \retval idx 1d index corresponding to position <tt>(x_idx, y_idx, ...)</tt>
   in a regular grid.
   \details \p idx goes through every x for each value of y, every y for each
   value of z, etc, in ascending order.
*/
template <typename Container>
auto positionToIndex(Container const &position, Container const &nr_points) {
  typename Container::value_type idx = 0;
  typename Container::value_type offset = 1;
  for (size_t dd = 0; dd < nr_points.size(); ++dd) {
    idx += offset * position[dd];
    offset *= nr_points[dd];
  }
  return idx;
}

/**
   \brief \p num are the digits of a basis \p basis number; increment it by 1.
   \note
   - Units digit in <tt>num[0]</tt>.
   - Maximum value wraps to <tt>00...0</tt>.
*/
template <typename Container, typename IntegerType>
void increment(Container &num, IntegerType basis) {
  bool carry;
  for (size_t dig = 0; dig < num.size(); ++dig) {
    carry = ++num[dig] / basis;
    num[dig] %= basis;
    if (!carry) {
      break;
    }
  }
}

/**
   \brief \p num are the digits of a number such that each digit ranges from 0
   to <tt>basis[dig] - 1</tt>; increment it by 1.
   \note
   - Units digit in <tt>num[0]</tt>.
   - Maximum value wraps to <tt>00...0</tt>.
*/
template <typename Container>
void increment(Container &num, Container const &basis) {
  bool carry;
  for (size_t dig = 0; dig < num.size(); ++dig) {
    carry = ++num[dig] / basis[dig];
    num[dig] %= basis[dig];
    if (!carry) {
      break;
    }
  }
}

/**
   \brief \p num are the digits of a number such that each digit ranges from 0
   to <tt>basis[dig] - 1</tt>; increment it by \p inc.
   \note
   - Units digit in <tt>num[0]</tt>.
   - Maximum value wraps to <tt>00...0</tt>.
*/
template <typename Container>
void increment(Container &num, Container const &basis, std::size_t inc) {
  for (std::size_t ii = 0; ii < inc; ++ii) {
    increment(num, basis);
  }
}
} // namespace modular

#endif /* GENERAL_MODULAR_H */
