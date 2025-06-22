/**
 \file General/Meta.h
 \author Tomas Aquino
 \date 23/03/2023
 \brief Template metamagic utilities.
*/

#ifndef GENERAL_META_H
#define GENERAL_META_H

#include <bits/utility.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <type_traits>
#include <utility>

/** \namespace meta Template metamagic utilities. */
namespace meta {
// The following has_ methods and associated machinery are based on:
// https://stackoverflow.com/questions/29772601/why-is-sfinae-causing-failure-when-there-are-two-functions-with-different-signat
/**
   \class types General/Meta.h "General/Meta.h"
   \brief Bundle of types.
*/
template <typename...> struct types { using type = types; };

/**
   \class has_member_impl General/Meta.h "General/Meta.h"
   \brief Implementation details for \c has_member.
*/
template <template <typename...> typename X, typename types, typename = void>
struct has_member_impl : std::false_type {};
/**
   \class has_member_impl General/Meta.h "General/Meta.h"
   \brief Implementation details for \c has_member.
*/
template <template <typename...> typename X, typename... Ts>
struct has_member_impl<X, types<Ts...>, std::void_t<X<Ts...>>>
    : std::true_type {};
/** \brief Check if \c X has a member. */
template <template <typename...> typename X, typename... Ts>
using has_member = has_member_impl<X, types<Ts...>>;

/** \class is_same General/Meta.h "General/Meta.h"
 \brief Check if X<Y> is the same type as Z. */
template <template <typename> typename X, typename Y, typename Z>
struct is_same {
  static constexpr bool value = std::is_same_v<X<Y>, Z>;
};
/** \brief Check if X<Y> is of void type. */
template <template <typename> typename X, typename Y>
using is_void = is_same<X, Y, void>;
/**
   \class is_arithmetic General/Meta.h "General/Meta.h"
   \brief Check if X<Y> is of arithmetic type.
*/
template <template <typename> typename X, typename Y> struct is_arithmetic {
  static constexpr bool value = std::is_arithmetic_v<X<Y>>;
};
/**
   \class is_integral General/Meta.h "General/Meta.h"
   \brief Check if X<Y> is of integral type.
*/
template <template <typename> typename X, typename Y> struct is_integral {
  static constexpr bool value = std::is_integral_v<X<Y>>;
};

/**\brief Type of \c X::print(). */
template <typename X> using print_t = decltype(std::declval<X>().print());
/** \brief Check if \c X has method \c print(). */
template <typename X> using has_print = has_member<print_t, X>;
/** \brief Check if \c X has method \c print(). */
template <typename X> inline constexpr bool has_print_v = has_print<X>::value;

/**\brief Type of \c X::begin(). */
template <typename X>
using begin_result_t = decltype(std::declval<X>().begin());
/** \brief Check if \c X has method \c begin(). */
template <typename X> using has_begin = has_member<begin_result_t, X>;
template <typename X>
/** \brief Check if \c X has method \c begin(). */
inline constexpr bool has_begin_v = has_begin<X>::value;

/**\brief Type of \c X::size(). */
template <typename X> using size_result_t = decltype(std::declval<X>().size());
/** \brief Check if \c X has method \c size(). */
template <typename X> using has_size = has_member<size_result_t, X>;
/** \brief Check if \c X has method \c size(). */
template <typename X> inline constexpr bool has_size_v = has_size<X>::value;

/**\brief Type of \c X::reserve(). */
template <typename X>
using reserve_result_t = decltype(std::declval<X>().reserve());
/** \brief Check if \c X has method \c reserve(). */
template <typename X> using has_reserve = has_member<reserve_result_t, X>;
template <typename X>
/** \brief Check if \c X has method \c reserve(). */
inline constexpr bool has_reserve_v = has_reserve<X>::value;

/**\brief Type of \c X::push_back(). */
template <typename X>
using push_back_result_t = decltype(std::declval<X>().reserve());
/** \brief Check if \c X has method \c push_back(). */
template <typename X> using has_push_back = has_member<push_back_result_t, X>;
template <typename X>
/** \brief Check if \c X has method \c push_back(). */
inline constexpr bool has_push_back_v = has_push_back<X>::value;

/**\brief Type of static method \c X::info(std::cout). */
template <typename X> using static_info_result_t = decltype(X::info(std::cout));
/** \brief Check if \c X can call static method \c info(std::cout). */
template <typename X>
using has_static_info = has_member<static_info_result_t, X>;
template <typename X>
/** \brief Check if \c X can call static method \c info(std::cout). */
inline constexpr bool has_static_info_v = has_static_info<X>::value;

/**\brief Type of member \c X::Parameters. */
template <typename X> using parameters_type_t = typename X::Parameters;
/** \brief Check if \c X defines type \c Parameters. */
template <typename X>
using has_parameters_type = has_member<parameters_type_t, X>;
template <typename X>
/** \brief Check if \c X defines type \c Parameters. */
inline constexpr bool has_parameters_type_v = has_parameters_type<X>::value;

// The following can_call methods are adapted from Passer By's answer here:
// https://stackoverflow.com/questions/51404763/c-compile-time-check-that-an-overloaded-function-can-be-called-with-a-certain
template <typename = void, typename... Args>
struct can_call_sqrt : std::false_type {};
template <typename... Args>
/**
   \class can_call_sqrt General/Meta.h "General/Meta.h"
   \brief Check whether class can call abs.
*/
struct can_call_sqrt<std::void_t<decltype(std::sqrt(std::declval<Args>()...))>,
                     Args...> : std::true_type {};
/** \brief To check whether class can call sqrt. */
template <typename... Args>
inline constexpr bool can_call_sqrt_v = can_call_sqrt<void, Args...>::value;

template <typename = void, typename... Args>
struct can_call_abs : std::false_type {};
/**
   \class can_call_abs General/Meta.h "General/Meta.h"
   \brief Check whether class can call abs.
*/
template <typename... Args>
struct can_call_abs<std::void_t<decltype(std::abs(std::declval<Args>()...))>,
                    Args...> : std::true_type {};
/** \brief Check whether class can call abs. */
template <typename... Args>
inline constexpr bool can_call_abs_v = can_call_abs<void, Args...>::value;

// The following has_operator methods and associated machinery are adapted from
// Richard Hodges's answer here:
// https://stackoverflow.com/questions/51404763/c-compile-time-check-that-an-overloaded-function-can-be-called-with-a-certain
// .

/** \brief Check whether classes can be used with binary operator. */
template <typename X, typename Y, typename Op> struct op_valid {
  template <typename U, typename L, typename R>
  static auto test(int)
      -> decltype(std::declval<U>()(std::declval<L>(), std::declval<R>()),
                  void(), std::true_type());

  template <typename U, typename L, typename R>
  static auto test(...) -> std::false_type;

  using type = decltype(test<Op, X, Y>(0));
};
/** \brief Check whether classes can be used with binary operator. */
template <typename X, typename Y, typename Op>
using op_valid_t = typename op_valid<X, Y, Op>::type;

/** \namespace notstd Operators not in std. */
namespace notstd {
/**
   \struct left_shift General/Meta.h "General/Meta.h"
   \brief Left shift operator.
*/
struct left_shift {
  template <typename L, typename R>
  constexpr auto operator()(L &&l, R &&r) const
      noexcept(noexcept(std::forward<L>(l) << std::forward<R>(r)))
          -> decltype(std::forward<L>(l) << std::forward<R>(r)) {
    return std::forward<L>(l) << std::forward<R>(r);
  }
};

/**
   \struct right_shift General/Meta.h "General/Meta.h"
   \brief Right shift operator.
*/
struct right_shift {
  template <typename L, typename R>
  constexpr auto operator()(L &&l, R &&r) const
      noexcept(noexcept(std::forward<L>(l) >> std::forward<R>(r)))
          -> decltype(std::forward<L>(l) >> std::forward<R>(r)) {
    return std::forward<L>(l) >> std::forward<R>(r);
  }
};

/**
   \struct operator_and General/Meta.h "General/Meta.h"
   \brief And operator.
*/
struct operator_and {
  template <typename L, typename R>
  constexpr auto operator()(L &&l, R &&r) const
      noexcept(noexcept(std::forward<L>(l) & std::forward<R>(r)))
          -> decltype(std::forward<L>(l) & std::forward<R>(r)) {
    return std::forward<L>(l) & std::forward<R>(r);
  }
};
} // namespace notstd
/** \brief Check whether classes can be used with equality operator. */
template <typename X, typename Y>
using has_equality = op_valid_t<X, Y, std::equal_to<>>;
/** \brief Check whether classes can be used with equality operator. */
template <typename X, typename Y>
inline constexpr bool has_equality_v = has_equality<X, Y>::value;
/** \brief Check whether classes can be used with inequality operator. */
template <typename X, typename Y>
using has_inequality = op_valid_t<X, Y, std::not_equal_to<>>;
/** \brief Check whether classes can be used with inequality operator. */
template <typename X, typename Y>
inline constexpr bool has_inequality_v = has_inequality<X, Y>::value;
/** \brief Check whether classes can be used with less than operator. */
template <typename X, typename Y>
using has_less_than = op_valid_t<X, Y, std::less<>>;
/** \brief Check whether classes can be used with less than operator. */
template <typename X, typename Y>
inline constexpr bool has_less_than_v = has_less_than<X, Y>::value;
/**
   \brief Check whether classes can be used with less than or equal operator.
*/
template <typename X, typename Y>
using has_less_equal = op_valid_t<X, Y, std::less_equal<>>;
/**
   \brief Check whether classes can be used with less than or equal operator.
*/
template <typename X, typename Y>
inline constexpr bool has_less_equal_v = has_less_equal<X, Y>::value;
/** \brief Check whether classes can be used with greater than operator. */
template <typename X, typename Y>
using has_greater_than = op_valid_t<X, Y, std::greater<>>;
/** \brief Check whether classes can be used with greater than operator. */
template <typename X, typename Y>
inline constexpr bool has_greater_than_v = has_greater_than<X, Y>::value;
/**
   \brief Check whether classes can be used with greater than or equal operator.
*/
template <typename X, typename Y>
using has_greater_equal = op_valid_t<X, Y, std::greater_equal<>>;
/**
   \brief Check whether classes can be used with greater than or equal operator.
*/
template <typename X, typename Y>
inline constexpr bool has_greater_equal_v = has_greater_equal<X, Y>::value;
/** \brief Check whether classes can be used with bit xor operator. */
template <typename X, typename Y>
using has_bit_xor = op_valid_t<X, Y, std::bit_xor<>>;
/** \brief Check whether classes can be used with bit xor operator. */
template <typename X, typename Y>
inline constexpr bool has_bit_xor_v = has_bit_xor<X, Y>::value;
/** \brief Check whether classes can be used with bit or operator. */
template <typename X, typename Y>
using has_bit_or = op_valid_t<X, Y, std::bit_or<>>;
/** \brief Check whether classes can be used with bit or operator. */
template <typename X, typename Y>
inline constexpr bool has_bit_or_v = has_bit_or<X, Y>::value;
/** \brief Check whether classes can be used with left shift operator. */
template <typename X, typename Y>
using has_left_shift = op_valid_t<X, Y, notstd::left_shift>;
/** \brief Check whether classes can be used with left shift operator. */
template <typename X, typename Y>
inline constexpr bool has_left_shift_v = has_left_shift<X, Y>::value;
/** \brief Check whether classes can be used with right shift operator. */
template <typename X, typename Y>
using has_right_shift = op_valid_t<X, Y, notstd::right_shift>;
/** \brief Check whether classes can be used with right shift operator. */
template <typename X, typename Y>
inline constexpr bool has_right_shift_v = has_right_shift<X, Y>::value;
/** \brief Check whether classes can be used with plus operator. */
template <typename X, typename Y>
using has_plus = op_valid_t<X, Y, std::plus<>>;
/** \brief Check whether classes can be used with plus operator. */
template <typename X, typename Y>
inline constexpr bool has_plus_v = has_plus<X, Y>::value;
/** \brief Check whether classes can be used with minus operator. */
template <typename X, typename Y>
using has_minus = op_valid_t<X, Y, std::minus<>>;
/** \brief Check whether classes can be used with minus operator. */
template <typename X, typename Y>
inline constexpr bool has_minus_v = has_minus<X, Y>::value;
/** \brief Check whether classes can be used with multiplication operator. */
template <typename X, typename Y>
using has_multiplies = op_valid_t<X, Y, std::multiplies<>>;
/** \brief Check whether classes can be used with multiplication operator. */
template <typename X, typename Y>
inline constexpr bool has_multiplies_v = has_multiplies<X, Y>::value;
/** \brief Check whether classes can be used with division operator. */
template <typename X, typename Y>
using has_divides = op_valid_t<X, Y, std::divides<>>;
/** \brief Check whether classes can be used with division operator. */
template <typename X, typename Y>
inline constexpr bool has_divides_v = has_divides<X, Y>::value;
/** \brief Check whether classes can be used with and operator. */
template <typename X, typename Y>
using has_and = op_valid_t<X, Y, notstd::operator_and>;
/** \brief Check whether classes can be used with and operator. */
template <typename X, typename Y>
inline constexpr bool has_and_v = has_and<X, Y>::value;

/**
   \class is_constructible_from_size_impl General/Meta.h "General/Meta.h"
   \brief Check if \c X<Y>(size()) is viable.
*/
template <typename X> struct is_constructible_from_size_impl {
  static constexpr bool value =
      std::is_constructible_v<X, decltype(std::declval<X>().size())>;
};
/** \brief Check if \c X<Y>(size()) is viable. */
template <typename X>
inline constexpr bool is_constructible_from_size_v =
    std::conjunction_v<has_size<X>, is_constructible_from_size_impl<X>>;

/**
   \class is_convertible_from_times_impl General/Meta.h "General/Meta.h"
   \brief Check if \c Z=X*Y is viable.
*/
template <typename X, typename Y, typename Z>
struct is_convertible_from_multiplies_impl {
  static constexpr bool value =
      std::is_convertible_v<decltype(std::declval<X>() * std::declval<Y>()), Z>;
};
/** \brief Check if \c Z=X*Y is viable. */
template <typename X, typename Y, typename Z>
inline constexpr bool is_convertible_from_multiplies_v =
    std::conjunction_v<has_multiplies<X, Y>,
                       is_convertible_from_multiplies_impl<X, Y, Z>>;

/**
   \class is_convertible_from_plus_impl General/Meta.h "General/Meta.h"
   \brief Check if \c Z=X+Y is viable.
*/
template <typename X, typename Y, typename Z>
struct is_convertible_from_plus_impl {
  static constexpr bool value =
      std::is_convertible_v<decltype(std::declval<X>() + std::declval<Y>()), Z>;
};
/** \brief Check if \c Z=X+Y is viable. */
template <typename X, typename Y, typename Z>
inline constexpr bool is_convertible_from_plus_v =
    std::conjunction_v<has_plus<X, Y>, is_convertible_from_plus_impl<X, Y, Z>>;

/**
   \class is_convertible_from_minus_impl General/Meta.h "General/Meta.h"
   \brief Check if \c Z=X-Y is viable.
*/
template <typename X, typename Y, typename Z>
struct is_convertible_from_minus_impl {
  static constexpr bool value =
      std::is_convertible_v<decltype(std::declval<X>() - std::declval<Y>()), Z>;
};
/** \brief Check if \c Z=X-Y is viable. */
template <typename X, typename Y, typename Z>
inline constexpr bool is_convertible_from_minus_v =
    std::conjunction_v<has_minus<X, Y>,
                       is_convertible_from_minus_impl<X, Y, Z>>;

/**
   \class is_convertible_from_divides_impl General/Meta.h "General/Meta.h"
   \brief Check if \c Z=X/Y is viable.
*/
template <typename X, typename Y, typename Z>
struct is_convertible_from_divides_impl {
  static constexpr bool value =
      std::is_convertible_v<decltype(std::declval<X>() / std::declval<Y>()), Z>;
};
/** \brief Check if \c Z=X/Y is viable. */
template <typename X, typename Y, typename Z>
inline constexpr bool is_convertible_from_divides_v =
    std::conjunction_v<has_divides<X, Y>,
                       is_convertible_from_divides_impl<X, Y, Z>>;

/**
   \class is_convertible_from_abs_impl General/Meta.h "General/Meta.h"
   \brief Check if \c Y=std::abs(X) is viable.
*/
template <typename X, typename Y> struct is_convertible_from_abs_impl {
  static constexpr bool value =
      std::is_convertible_v<decltype(std::abs(std::declval<X>())), Y>;
};
/** \brief Check if \c Y=std::abs(X) is viable. */
template <typename X, typename Y>
inline constexpr bool is_convertible_from_abs_v =
    std::conjunction_v<can_call_abs<X>, is_convertible_from_abs_impl<X, Y>>;

/**
   \class is_convertible_from_pow_impl General/Meta.h "General/Meta.h"
   \brief Check if \c Z=std::pow(X,Y) is viable.
*/
template <typename X, typename Y, typename Z>
struct is_convertible_from_pow_impl {
  static constexpr bool value = std::is_convertible_v<
      decltype(std::pow(std::declval<X>(), std::declval<Y>())), Z>;
};
/** \brief Check if \c Z=std::pow(X, Y) is viable. */
template <typename X, typename Y, typename Z>
inline constexpr bool is_convertible_from_pow_v =
    std::conjunction_v<can_call_abs<X>, is_convertible_from_pow_impl<X, Y, Z>>;

/**
   \class is_convertible_from_sqrt_impl General/Meta.h "General/Meta.h"
   \brief Check if \c Y=std::sqrt(X) is viable.
*/
template <typename X, typename Y> struct is_convertible_from_sqrt_impl {
  static constexpr bool value =
      std::is_convertible_v<decltype(std::sqrt(std::declval<X>())), Y>;
};
/** \brief Check if \c Y=std::sqrt(X) is viable. */
template <typename X, typename Y>
inline constexpr bool is_convertible_from_sqrt_v =
    std::conjunction_v<can_call_sqrt<X>, is_convertible_from_sqrt_impl<X, Y>>;

/**
   \class indices General/Meta.h "General/Meta.h"
   \brief Indices for template metamagic.
*/
template <std::size_t... Indices> struct indices {
  using next = indices<Indices..., sizeof...(Indices)>;
};

/**
   \class Selector_t General/Meta.h "General/Meta.h"
   \brief Type for selecting function implementations at compile time.
*/
template <typename TT> struct Selector_t { using type = TT; };

/**
   \class Selector General/Meta.h "General/Meta.h"
   \brief Type for selecting function implementations at compile time.
*/
template <typename TT, TT val> struct Selector { using type = TT; };

/** \brief Type-dependent expression that is always false. */
template <typename> constexpr bool dependent_false_v = false;

/** \brief Type-dependent expression that is always true. */
template <typename> constexpr bool dependent_true_v = true;

/** \brief Check if type T is an instance of template type U with arbitrary
 * template parameters. */
template <typename T, template <typename...> typename U>
inline constexpr bool is_instance_of_v = std::false_type{};

/** \brief Check if type T is an instance of type U<Vs...>. */
template <template <typename...> typename U, typename... Vs>
inline constexpr bool is_instance_of_v<U<Vs...>, U> = std::true_type{};

/**
   \class Empty General/Meta.h "General/Meta.h"
   \brief Class holding nothing and doing nothing.
*/
struct Empty {
  template <typename... Args> Empty(Args...) {}
  Empty() {}
};

/**
   \struct DoNothing General/Meta.h "General/Meta.h"
   \brief Functor that does nothing.
*/
struct DoNothing {
  template <typename... Args> void operator()(Args...) const {}
};

/**
   \struct DoFalse General/Meta.h "General/Meta.h"
   \brief Functor that always returns false.
*/
struct DoFalse {
  template <typename... Args> bool operator()(Args...) const { return false; }
};

/**
   \struct DoTrue General/Meta.h "General/Meta.h"
   \brief Functor that always returns true.
*/
struct DoTrue {
  template <typename... Args> bool operator()(Args...) const { return true; }
};

} // namespace meta

#endif /* GENERAL_META_H */
