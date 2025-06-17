/**
   \file PTOF/Meta.h
   \author Tomás Aquino
   \date 23/03/2023
   \brief Template metamagic utilities.
*/

#ifndef PTOF_META_H
#define PTOF_META_H

#include "General/Meta.h"
#include <label.H>
#include <type_traits>
#include <utility>

namespace meta {
/**\brief Type of \c X::velocity_rescaling_factor. */
template <typename X>
using velocity_rescaling_factor_t =
    decltype(std::declval<X>().velocity_rescaling_factor);
/**
   \brief Check if \c X has member <tt>double velocity_rescaling_factor</tt>.
*/
template <typename X>
inline constexpr bool has_velocity_rescaling_factor_v =
    std::conjunction_v<has_member<velocity_rescaling_factor_t, X>,
                       is_same<velocity_rescaling_factor_t, X, double>>;

/**\brief Type of \c X::face. */
template <typename X> using face_t = decltype(std::declval<X>().face);
/**
   \brief Check if \c X has member <tt>Foam::label face</tt>.
*/
template <typename X>
inline constexpr bool has_face_v =
    std::conjunction_v<has_member<face_t, X>,
                       is_same<face_t, X, Foam::label>>;

/**\brief Type of \c X::absorbed. */
template <typename X> using absorbed_t = decltype(std::declval<X>().absorbed);
/**
   \brief Check if \c X has member <tt>bool absorbed</tt>.
*/
template <typename X>
inline constexpr bool has_absorbed_v =
    std::conjunction_v<has_member<absorbed_t, X>, is_same<absorbed_t, X, bool>>;

/**\brief Type of \c X::adsorbed. */
template <typename X> using adsorbed_t = decltype(std::declval<X>().adsorbed);
/**
   \brief Check if \c X has member <tt>bool adsorbed</tt>.
*/
template <typename X>
inline constexpr bool has_adsorbed_v =
    std::conjunction_v<has_member<adsorbed_t, X>, is_same<adsorbed_t, X, bool>>;
} // namespace meta

#endif /* PTOF_META_H */
