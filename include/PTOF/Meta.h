/**
   \file PTOF/Meta.h
   \author Tomas Aquino
   \date 23/03/2023
   \brief Template metamagic utilities.
*/

#ifndef PTOF_META_H
#define PTOF_META_H

#include "General/Meta.h"
#include <label.H>
#include <point.H>
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

/**\brief Type of \c X::type. */
template <typename X> using type_t = decltype(std::declval<X>().type);
/**
   \brief Check if \c X has member <tt>Foam::label type</tt>.
*/
template <typename X>
inline constexpr bool has_type_v =
    std::conjunction_v<has_member<type_t, X>, is_same<type_t, X, Foam::label>>;

/**\brief Type of \c X::boundary_face. */
template <typename X> using boundary_face_t = decltype(std::declval<X>().face);
/**
   \brief Check if \c X has member <tt>Foam::label boundary_face</tt>.
*/
template <typename X>
inline constexpr bool has_boundary_face_v =
    std::conjunction_v<has_member<boundary_face_t, X>,
                       is_same<boundary_face_t, X, Foam::label>>;

/**\brief Type of \c X::contact_point. */
template <typename X>
using contact_point_t = decltype(std::declval<X>().contact_point);
/**
   \brief Check if \c X has member <tt>Foam::point contact_point</tt>.
*/
template <typename X>
inline constexpr bool has_contact_point_v =
    std::conjunction_v<has_member<contact_point_t, X>,
                       is_same<contact_point_t, X, Foam::point>>;

/**\brief Type of \c X::reinjections. */
template <typename X>
using reinjections_t = decltype(std::declval<X>().reinjections);
/**
   \brief Check if \c X has member <tt>std::size_t reinjections</tt>.
*/
template <typename X>
inline constexpr bool has_reinjections_v =
    std::conjunction_v<has_member<reinjections_t, X>,
                       is_same<reinjections_t, X, std::size_t>>;

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

/**\brief Type of <tt>X::boundaryField()</tt>. */
template <typename X>
using boundaryField_t = decltype(std::declval<X>().boundaryField());
/**
   \brief Check if \c X has member <tt>boundaryField()</tt>.
*/
template <typename X>
inline constexpr bool has_boundaryField_v =
    has_member<boundaryField_t, X>::value;

/**\brief Type of <tt>X::operator[](Foam::label)</tt>. */
template <typename X>
using square_brackets_of_label_t =
    decltype(std::declval<X>()[std::declval<Foam::label>()]);
/**
   \brief Check if \c X has member <tt>operator[](Foam::label)</tt>.
*/
template <typename X>
inline constexpr bool has_square_brackets_of_label_v =
    has_member<square_brackets_of_label_t, X>::value;

/**\brief Type of <tt>X::boundaryField()</tt>. */
template <typename X>
using boundaryField_of_void_t = decltype(std::declval<X>().boundaryField());
/**
   \brief Check if \c X has member <tt>boundaryField()</tt>.
*/
template <typename X>
inline constexpr bool has_boundaryField_of_void_v =
    has_member<boundaryField_of_void_t, X>::value;
} // namespace meta

#endif /* PTOF_META_H */
