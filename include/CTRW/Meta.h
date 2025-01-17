/**
   \file CTRW/Meta.h
   \author Tomás Aquino
   \date 23/03/2023
   \brief Template metamagic utilities.
*/

#ifndef CTRW_META_H
#define CTRW_META_H

#include "General/Meta.h"
#include "General/Useful.h"
#include <type_traits>
#include <utility>

namespace meta {
/** \brief Type of \c X::time_step(double). */
template <typename X>
using time_step_setter_double_result_t =
    decltype(std::declval<X>().time_step(std::declval<double>()));
/** \brief Type of \c X::time_step(int). */
template <typename X>
using time_step_setter_int_result_t =
    decltype(std::declval<X>().time_step(std::declval<int>()));
/** \brief Type of \c X::time_step(std::size_t). */
template <typename X>
using time_step_setter_uint_result_t =
    decltype(std::declval<X>().time_step(std::declval<std::size_t>()));
/**
   \brief Check if a class has method <tt>void time_step([double, int, or
   std::size_t])</tt>.
*/
template <typename X>
inline constexpr bool has_time_step_setter_v = std::disjunction_v<
    std::conjunction<has_member<time_step_setter_double_result_t, X>,
                     is_void<time_step_setter_double_result_t, X>>,
    std::conjunction<has_member<time_step_setter_int_result_t, X>,
                     is_void<time_step_setter_int_result_t, X>>,
    std::conjunction<has_member<time_step_setter_uint_result_t, X>,
                     is_void<time_step_setter_uint_result_t, X>>>;

/** \brief Type of \c X::time_step(). */
template <typename X>
using time_step_getter_result_t = decltype(std::declval<X>().time_step());
/**
   \brief Check if a class has method <tt>Type time_step()</tt>, where \c Type
   is an arithmetic type.
*/
template <typename X>
inline constexpr bool has_time_step_getter_v =
    std::conjunction_v<has_member<time_step_getter_result_t, X>,
                       is_arithmetic<time_step_getter_result_t, X>>;

/** \brief Type of \c X::time_step. */
template <typename X> using time_step_t = decltype(std::declval<X>().time_step);
/**
   \brief Check if a class has member <tt>Type time_step()</tt>, where \c Type
   is an arithmetic type.
*/
template <typename X>
inline constexpr bool has_time_step_v =
    std::conjunction_v<has_member<time_step_t, X>,
                       is_arithmetic<time_step_t, X>>;

/** \brief Type of \c X::periodicity. */
template <typename X>
using periodicity_t = decltype(std::declval<X>().periodicity);
/** \brief Check if a class has member <tt>periodicity</tt>. */
template <typename X>
inline constexpr bool has_periodicity_v = has_member<periodicity_t, X>::value;

/** \brief Type of \c X::cell. */
template <typename X> using cell_t = decltype(std::declval<X>().cell);
/** \brief Check if a class has member <tt>integral_type cell</tt>. */
template <typename X>
inline constexpr bool has_cell_v =
    std::conjunction_v<has_member<cell_t, X>, is_integral<cell_t, X>>;

/** \brief Type of \c X::mass. */
template <typename X> using mass_t = decltype(std::declval<X>().mass);
/** \brief Check if a class has member <tt>arithmetic_type mass</tt>. */
template <typename X>
inline constexpr bool has_mass_v =
    std::conjunction_v<has_member<mass_t, X>, is_arithmetic<mass_t, X>>;

/** \brief Type of \c X::tag. */
template <typename X> using tag_t = decltype(std::declval<X>().tag);
/** \brief Check if a class has member <tt>integral_type tag</tt>. */
template <typename X>
inline constexpr bool has_tag_v =
    std::conjunction_v<has_member<tag_t, X>, is_integral<tag_t, X>>;

/** \brief Type of \c X::time. */
template <typename X> using time_t = decltype(std::declval<X>().time);
/** \brief Check if a class has member <tt>time</tt>. */
template <typename X>
inline constexpr bool has_time_v = has_member<time_t, X>::value;

/** \class Parameters_or_empty General/Meta.h "General/Meta.h"
    \brief  Get Parameters type for true or Empty type for false. */
template <bool, typename> struct Parameters_or_empty { using type = Empty; };
template <typename T> struct Parameters_or_empty<true, T> {
  using type = typename T::Parameters;
};
template <bool condition, typename T>
using Parameters_or_empty_t = typename Parameters_or_empty<condition, T>::type;
} // namespace meta

#endif /* CTRW_META_H */
