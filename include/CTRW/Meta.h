/**
 * @file   Meta.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Thu Mar 23 00:00:00 2023
 *
 * @brief Template metamagic utilities relating to CTRW.
 */

#ifndef CTRW_META_H
#define CTRW_META_H

#include "General/Meta.h"
#include <type_traits>
#include <utility>

namespace meta {
/** @brief Type of \c X::time_step(double). */
template <typename X>
using time_step_setter_double_result_t =
    decltype(std::declval<X>().time_step(std::declval<double>()));

/** @brief Type of \c X::time_step(int). */
template <typename X>
using time_step_setter_int_result_t =
    decltype(std::declval<X>().time_step(std::declval<int>()));

/** @brief Type of \c X::time_step(std::size_t). */
template <typename X>
using time_step_setter_uint_result_t =
    decltype(std::declval<X>().time_step(std::declval<std::size_t>()));

/**
 * @brief Check if \c X has method <tt>void time_step([double, int, or
 *        std::size_t])</tt>.
 */
template <typename X>
inline constexpr bool has_time_step_setter_v = std::disjunction_v<
    std::conjunction<has_member<time_step_setter_double_result_t, X>,
                     is_void<time_step_setter_double_result_t, X>>,
    std::conjunction<has_member<time_step_setter_int_result_t, X>,
                     is_void<time_step_setter_int_result_t, X>>,
    std::conjunction<has_member<time_step_setter_uint_result_t, X>,
                     is_void<time_step_setter_uint_result_t, X>>>;

/** @brief Type of \c X::time_step(). */
template <typename X>
using time_step_getter_result_t = decltype(std::declval<X>().time_step());

/**
 * @brief Check if \c X has method <tt>Type time_step()</tt>, where \c Type is
 *        an arithmetic type.
 */
template <typename X>
inline constexpr bool has_time_step_getter_v =
    std::conjunction_v<has_member<time_step_getter_result_t, X>,
                       is_arithmetic<time_step_getter_result_t, X>>;

/** @brief Type of \c X::time_step. */
template <typename X> using time_step_t = decltype(std::declval<X>().time_step);

/**
 * @brief Check if \c X has member <tt>Type time_step()</tt>, where \c Type is
 *        an arithmetic type.
 */
template <typename X>
inline constexpr bool has_time_step_v =
    std::conjunction_v<has_member<time_step_t, X>,
                       is_arithmetic<time_step_t, X>>;

/** @brief Type of \c X::periodicity. */
template <typename X>
using periodicity_t = decltype(std::declval<X>().periodicity);

/** @brief Check if \c X has member <tt>periodicity</tt>. */
template <typename X>
inline constexpr bool has_periodicity_v = has_member<periodicity_t, X>::value;

/** @brief Type of \c X::cell. */
template <typename X> using cell_t = decltype(std::declval<X>().cell);

/** @brief Check if \c X has member <tt>integral_type cell</tt>. */
template <typename X>
inline constexpr bool has_cell_v =
    std::conjunction_v<has_member<cell_t, X>, is_integral<cell_t, X>>;

/** @brief Type of \c X::mass. */
template <typename X> using mass_t = decltype(std::declval<X>().mass);

/** @brief Check if \c X has member <tt>arithmetic_type mass</tt>. */
template <typename X>
inline constexpr bool has_mass_v =
    std::conjunction_v<has_member<mass_t, X>, is_arithmetic<mass_t, X>>;

/** @brief Type of \c X::tag. */
template <typename X> using tag_t = decltype(std::declval<X>().tag);

/** @brief Check if \c X has member <tt>integral_type tag</tt>. */
template <typename X>
inline constexpr bool has_tag_v =
    std::conjunction_v<has_member<tag_t, X>, is_integral<tag_t, X>>;

/** @brief Type of \c X::time. */
template <typename X> using time_t = decltype(std::declval<X>().time);

/** @brief Check if \c X has member <tt>time</tt>. */
template <typename X>
inline constexpr bool has_time_v = has_member<time_t, X>::value;

/** @brief Type of \c X::boundary_periodic. */
template <typename X>
using boundary_periodic_t = decltype(std::declval<X>().boundary_periodic);

/** @brief Check if \c X has member <tt>boundary_periodic</tt>. */
template <typename X>
inline constexpr bool has_boundary_periodic_v =
    has_member<boundary_periodic_t, X>::value;

/** @brief Helper for meta::Parameters_or_empty_t. */
template <bool condition, typename T> struct Parameters_or_empty {
  using type = Empty;
};

/** @brief Helper for meta::Parameters_or_empty_t. */
template <typename T> struct Parameters_or_empty<true, T> {
  using type = typename T::Parameters;
};

/**
 * @brief Get \c T::Parameters type when \c condition is \c true or meta::Empty
 *        type when \c condition is \c false.
 */
template <bool condition, typename T>
using Parameters_or_empty_t = typename Parameters_or_empty<condition, T>::type;
} // namespace meta

#endif /* CTRW_META_H */
