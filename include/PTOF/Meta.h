/**
 \file PTOF/Meta.h
 \author Tomás Aquino
 \date 23/03/2023
 \brief Template metamagic utilities.
*/

#ifndef PTOF_META_H
#define PTOF_META_H

#include <type_traits>
#include <utility>
#include "General/Meta.h"

/** \namespace useful Template metamagic utilities. */
namespace meta
{
  /**\brief Type of \c X::velocity_rescaling_factor. */
  template <typename X>
  using velocity_rescaling_factor_t
  = decltype(std::declval<X>().velocity_rescaling_factor);
  /** \brief has_size General/Useful.h "General/Useful.h"
  \brief Check if a class has member <tt>double velocity_rescaling_factor</tt>. */
  template <typename X>
  inline constexpr bool has_velocity_rescaling_factor_v
  = has_member<velocity_rescaling_factor_t, X>::value
    && std::is_same_v<velocity_rescaling_factor_t<X>, double>;
}

#endif /* PTOF_META_H */
