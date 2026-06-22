/**
 * @file   Useful.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Tue Mar 15 00:00:00 2011
 *
 * @brief Miscelaneous collection of useful objects and algorithms.
 */

#ifndef GENERAL_USEFUL_H
#define GENERAL_USEFUL_H

#include "General/Meta.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

/**
 * @namespace useful Miscelaneous collection of useful objects and algorithms.
 */
namespace useful {

/**
 * @brief Check if sorted \p container contains \p val.
 *
 * @note Container must be sorted.
 */
template <typename T, typename U>
bool contains(const std::vector<T> &container, U const &val) {
  auto it =
      std::lower_bound(container.begin(), container.end(), val,
                       [](T const &elem, U const &val) { return elem < val; });

  return it != container.end() && *it == val;
}

/**
 * @brief Check whether \p container contains \p value.
 *
 * @details Container must be sorted according to \p comp_less, equality is
 *          checked using \p comp_eq.
 */
template <typename T, typename U, typename Comp_less, typename Comp_eq>
bool contains(std::vector<T> const &container, U const &val,
              Comp_less comp_less = std::less<void>{},
              Comp_eq comp_eq = std::equal_to<void>{}) {
  auto it =
      std::lower_bound(container.begin(), container.end(), val, comp_less);

  return it != container.end() && comp_eq(*it, val);
}

/** @brief Check if \p string ends with \p suffix. */
inline bool ends_with(const std::string &string, const std::string &suffix) {
  return string.size() >= suffix.size() &&
         string.substr(string.size() - suffix.size()) == suffix;
}

/** @brief Class holding const object. */
template <typename Object_Type, typename Return_Type = Object_Type>
struct StoreConst {
  using value_type = Return_Type;

  Return_Type operator()() const { return obj; }
  const Object_Type obj;
};

/** @brief Class holding object. */
template <typename Object_Type, typename Return_Type = Object_Type>
struct Store {
  using value_type = Return_Type;

  Return_Type operator()() const { return obj; }
  Object_Type obj;
};

/**
 * @brief Make an object and initialize it given \p parameters.
 *
 * @note Object must implement:
 *
 * - <tt>initialize(Params const&)</tt>.
 */
template <typename Object, typename Params>
Object create(Params const &parameters = meta::Empty{}) {
  Object object;
  object.initialize(parameters);

  return object;
}

/**
 * @brief Functor to make objects of given type, passing arbitrary parameters.
 */
template <typename Object> struct Maker {
  template <typename... Args> Object operator()(Args... args) const {
    return Object{args...};
  }
};

/**
 * @brief Functor to make an object on the heap, passing arbitrary parameters.
 */
template <typename T> struct Creator {
  using value_type = T;
  using pointer = T *;

  template <typename... Args> pointer operator()(Args... args) const {
    return new T{args...};
  }
};

/** @brief Functor to forward an argument. */
template <typename Object> struct Forward {
  Object &&operator()(Object &&object) const {
    return std::forward<Object>(object);
  }
};

/** @brief Functor to forward a const reference. */
template <typename Object> struct Forward_ref {
  Object const &operator()(Object const &object) const { return object; }
};

/** @brief Convert <tt>()</tt> to <tt>[]</tt>. */
template <typename Object> struct SquareBracketsToRoundBrackets {
  SquareBracketsToRoundBrackets(Object &&object)
      : _object{std::forward<Object>(object)} {}

  template <typename TT> auto operator()(TT val) { return _object[val]; }

private:
  Object _object;
};
template <typename Object>
SquareBracketsToRoundBrackets(Object) -> SquareBracketsToRoundBrackets<Object>;

/** @brief Convert <tt>() const</tt> to <tt>[] const</tt>. */
template <typename Object> struct SquareBracketsToRoundBracketsConst {
  SquareBracketsToRoundBracketsConst(Object &&object)
      : _object{std::forward<Object>(object)} {}

  template <typename TT> auto operator()(TT val) const { return _object[val]; }

private:
  Object _object;
};
template <typename Object>
SquareBracketsToRoundBracketsConst(Object)
    -> SquareBracketsToRoundBracketsConst<Object>;

/** @brief Convert <tt>()</tt> to <tt>[]</tt>. */
template <typename Object> struct RoundBracketsToSquareBrackets {
  RoundBracketsToSquareBrackets(Object &&object)
      : _object{std::forward<Object>(object)} {}

  template <typename TT> auto operator[](TT val) { return _object(val); }

private:
  Object _object;
};
template <typename Object>
RoundBracketsToSquareBrackets(Object) -> RoundBracketsToSquareBrackets<Object>;

/** @brief Convert <tt>() const</tt> to <tt>[] const</tt>. */
template <typename Object> struct RoundBracketsToSquareBracketsConst {
  RoundBracketsToSquareBracketsConst(Object &&object)
      : _object{std::forward<Object>(object)} {}

  template <typename TT> auto operator[](TT val) const { return _object(val); }

private:
  Object _object;
};
template <typename Object>
RoundBracketsToSquareBracketsConst(Object)
    -> RoundBracketsToSquareBracketsConst<Object>;

/** @brief Simple hash for a container by combining element hashes. */
template <typename Container> struct hash_container {
  std::size_t operator()(Container const &container) const {
    using value_type = typename Container::value_type;
    std::size_t hash_val = std::hash<value_type>()(container[0]);
    for (size_t ii = 1; ii < container.size() - 1; ++ii) {
      hash_val =
          (hash_val ^ (std::hash<value_type>()(container[ii]) << 1)) >> 1;
    }
    hash_val = hash_val ^
               (std::hash<value_type>()(container[container.size() - 1]) << 1);

    return hash_val;
  }
};

/** @brief Same as \c std::remove_cvref (only available from C++20). */
template <class T> struct remove_cvref {
  using type = std::remove_cv_t<std::remove_reference_t<T>>;
};

/** @brief Same as \c std::remove_cvref_t (only available from C++20). */
template <class T> using remove_cvref_t = typename remove_cvref<T>::type;

/** @brief Combine \p seed with the hash of an object \p v. */
template <typename T> void hash_combine(std::size_t &seed, T const &vv) {
  if constexpr (meta::has_begin_v<T>) {
    std::hash<typename T::value_type> hasher;
    for (auto const &val : vv) {
      seed ^= hasher(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
  } else {
    std::hash<T> hasher;
    seed ^= hasher(vv) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }
}

/** @brief Combine \p seed with the hash of an object \p v. */
template <typename T1, typename T2>
void hash_combine(std::size_t &seed, std::pair<T1, T2> const &vv) {
  hash_combine(seed, vv.first);
  hash_combine(seed, vv.second);
}

/** @brief Hash for \c std::pair. */
template <typename S, typename T> struct Hash_pair {
  std::size_t operator()(std::pair<S, T> const &vv) const {
    std::size_t seed = 0;
    hash_combine(seed, vv);
    return seed;
  }
};

/** @brief Hash for containers. */
template <typename Container> struct Hash_container {
  std::size_t operator()(Container const &vv) const {
    std::size_t seed = 0;
    for (auto const &val : vv) {
      hash_combine(seed, vv);
    }
    return seed;
  }
};

/** @brief Hash for combining types. */
template <typename... T> struct Hash_combine {
  std::size_t operator()(T... args) const {
    std::size_t seed = 0;
    (hash_combine(seed, args), ...);
    return seed;
  }
};

/** @return Sign of \p val. */
template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

/** @brief If \p to_replace is \c NaN, replace it with \p replace_with. */
template <typename T> T &deNaN(T &to_replace, T replace_with = T()) {
  if constexpr (meta::has_begin_v<T>) {
    for (auto const &val : to_replace) {
      denan(val, replace_with);
    }
  } else {
    if (to_replace != to_replace) {
      to_replace = replace_with;
    }
  }
  return to_replace;
}

/**
 * @brief If \p to_chop is smaller than \p tolerance, replace it with \p
 *        replace_with.
 */
template <typename T> T &chop(T &to_chop, T tolerance, T replace_with = T()) {
  if constexpr (meta::has_begin_v<T>) {
    for (auto const &val : to_chop) {
      denan(val, replace_with);
    }
  } else {
    if (std::abs(to_chop - replace_with) < tolerance) {
      to_chop = replace_with;
    }
  }
  return to_chop;
}

/**
 * @brief Copy end element of \p container to index \p idx, then erase end
 *        element.
 */
template <typename Container>
Container &swap_erase(Container &container, std::size_t idx) {
  if (idx != container.size() - 1) {
    container[idx] = container.back();
  }
  container.pop_back();
  return container;
}

/**
 * @brief Swap-erase elements in \p container at indices \p indices.
 *
 * @note \p positions will be sorted, and must implement \c
 *       sort(std::greater<std::size_t>).
 */
template <typename Container, typename List>
Container &swap_erase(Container &container, List &&indices) {
  indices.sort(std::greater<std::size_t>{});
  for (auto const position : indices) {
    swap_erase(container, indices);
  }
  return container;
}

/** @brief Swap-erase elements of \p container satisfying \p criterion. */
template <typename Container, typename Criterion>
Container &swap_erase_if(Container &container, Criterion &&criterion) {
  for (std::size_t ii = container.size(); ii-- > 0;) {
    if (criterion(container[ii])) {
      swap_erase(container, ii);
    }
  }
  return container;
}

/**
 * @brief Like swap_erase(), but call \c delete on the
 *        removed element.
 */
template <typename Container>
Container &swap_delete(Container &container, std::size_t idx) {
  if (idx != container.size() - 1) {
    container[idx] = container.back();
  }
  delete container.back();
  container.pop_back();
  return container;
}

/** @brief Swap-delete elements in \p container at positions \p positions. */
template <typename Container, typename List>
Container &swap_delete(Container &container, List &indices) {
  indices.sort(std::greater<std::size_t>{});
  for (auto const position : indices) {
    swap_delete(container, position);
  }

  return container;
}

/** @brief Swap-delete elements of \p container satisfying \p criterion. */
template <typename Container, typename Criterion>
Container &swap_delete_if(Container &container, Criterion &&criterion) {
  for (std::size_t ii = container.size(); ii-- > 0;) {
    if (criterion(container[ii])) {
      swap_delete(container, ii);
    }
  }
  return container;
}

/** @brief Check if string represents a number. */
bool is_numeric(std::string const &str) {
  double result{};
  auto val = std::istringstream(str);
  val >> result;
  return !val.fail() && val.eof();
}

/**
 * @brief Throw if container \c container does not have \c nr_after_index
 *        elements after and including \c index.
 */
template <typename Container>
void check_size(Container const &container, std::size_t index,
                std::size_t nr_after_index) {
  if (container.size() < index + nr_after_index) {
    throw std::runtime_error{"check_size : Not enough elements"};
  }
}

/**
 * @return Longest prefix suffix array of \c pattern (for KMP).
 */
std::vector<std::size_t> build_lps(std::string const &pattern) {
  std::vector<std::size_t> lps(pattern.size());
  std::size_t len = 0;

  lps[0] = 0;
  for (std::size_t ii = 1; ii < pattern.length();) {
    if (pattern[ii] == pattern[len]) {
      ++len;
      lps[ii++] = len;
    } else {
      if (len != 0) {
        len = lps[len - 1];
      } else {
        lps[ii++] = 0;
      }
    }
  }

  return lps;
}

/**
 * @return Starting indices of \c pattern matches in \c string (using KMP).
 */
std::vector<std::size_t> pattern_matches(std::string const &pattern,
                                         std::string const &string) {
  auto lps = build_lps(pattern);
  std::vector<std::size_t> indices;

  for (std::size_t ii = 0, jj = 0; ii < string.length();) {
    if (string[ii] == pattern[jj]) {
      ++ii;
      ++jj;
      if (jj == pattern.length()) {
        indices.push_back(ii - jj);
        jj = lps[jj - 1];
      }
    } else {
      if (jj != 0)
        jj = lps[jj - 1];
      else
        ++ii;
    }
  }
  return indices;
}

/**
 * @return Number of \c pattern matches in \c string (using KMP).
 */
std::size_t pattern_matches_count(std::string const &pattern,
                                  std::string const &string) {
  auto lps = build_lps(pattern);
  std::size_t nr_matches = 0;

  for (std::size_t ii = 0, jj = 0; ii < string.length();) {
    if (string[ii] == pattern[jj]) {
      ++ii;
      ++jj;
      if (jj == pattern.length()) {
        ++nr_matches;
        jj = lps[jj - 1];
      }
    } else {
      if (jj != 0)
        jj = lps[jj - 1];
      else
        ++ii;
    }
  }
  return nr_matches;
}

/**
 * @brief Left shift operator.
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
 * @brief Right shift operator.
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
 * @brief And operator.
 */
struct operator_and {
  template <typename L, typename R>
  constexpr auto operator()(L &&l, R &&r) const
      noexcept(noexcept(std::forward<L>(l) & std::forward<R>(r)))
          -> decltype(std::forward<L>(l) & std::forward<R>(r)) {
    return std::forward<L>(l) & std::forward<R>(r);
  }
};
} // namespace useful

#endif /* GENERAL_USEFUL_H */
