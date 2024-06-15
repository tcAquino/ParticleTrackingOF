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

/** \namespace useful Template metamagic utilities. */
namespace meta {
/**\brief Type of \c X::time_step(double). */
template <typename X>
using time_step_setter_result_t =
    decltype(std::declval<X>().time_step(std::declval<double>()));
/** \brief Check if a class has method <tt>void time_step(double)</tt>. */
template <typename X>
inline constexpr bool has_time_step_setter_v =
    std::conjunction_v<has_member<time_step_setter_result_t, X>,
                       is_void<time_step_setter_result_t, X>>;

/**\brief Type of \c X::time_step(). */
template <typename X>
using time_step_getter_result_t = decltype(std::declval<X>().time_step());
/** \brief Check if a class has method <tt>Type time_step()</tt>, where \c Type
 * is an arithmetic type. */
template <typename X>
inline constexpr bool has_time_step_getter_v =
    std::conjunction_v<has_member<time_step_getter_result_t, X>,
                       is_arithmetic<time_step_getter_result_t, X>>;

/**\brief Type of \c X::time_step. */
template <typename X> using time_step_t = decltype(std::declval<X>().time_step);
/** \brief Check if a class has member <tt>Type time_step()</tt>, where \c Type
 * is an arithmetic type. */
template <typename X>
inline constexpr bool has_time_step_v =
    std::conjunction_v<has_member<time_step_t, X>,
                       is_arithmetic<time_step_t, X>>;

/**\brief Type of \c X::periodicity. */
template <typename X>
using periodicity_t = decltype(std::declval<X>().periodicity);
/** \brief Check if a class has member <tt>std::vector<int> periodicity</tt>. */
template <typename X>
inline constexpr bool has_periodicity_v =
    std::conjunction_v<has_member<periodicity_t, X>,
                       is_same<periodicity_t, X, std::vector<int>>>;

/**\brief Type of \c X::cell. */
template <typename X>
using cell_t = decltype(std::declval<X>().cell);
/** \brief Check if a class has member <tt>integral_type cell</tt>. */
template <typename X>
inline constexpr bool has_cell_v =
    std::conjunction_v<has_member<cell_t, X>,
                       is_integral<cell_t, X>>;

/** \class Parameters_or_empty General/Meta.h "General/Meta.h"
 \brief  Get Parameters type for true or Empty type for false. */
template <bool, typename> struct Parameters_or_empty {
  using type = useful::Empty;
};
template <typename T> struct Parameters_or_empty<true, T> {
  using type = typename T::Parameters;
};
template <bool condition, typename T>
using Parameters_or_empty_t = typename Parameters_or_empty<condition, T>::type;
} // namespace meta

#endif /* CTRW_META_H */
