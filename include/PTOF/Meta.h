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

/**\brief Type of \c X::patch_id. */
template <typename X> using patch_id_t = decltype(std::declval<X>().patch_id);
/**
   \brief Check if \c X has member <tt>Foam::label patch_id</tt>.
*/
template <typename X>
inline constexpr bool has_patch_id_v =
    std::conjunction_v<has_member<patch_id_t, X>,
                       is_same<patch_id_t, X, Foam::label>>;
} // namespace meta

#endif /* PTOF_META_H */
