/**
 \file General/Useful.h
 \author Tomás Aquino
 \date 03/15/2011
 \brief Miscelaneous collection of useful objects and algorithms.
*/

#ifndef GENERAL_USEFUL_H
#define GENERAL_USEFUL_H

#include "General/IO.h"
#include "General/Meta.h"
#include <algorithm>
#include <bits/chrono.h>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/**
   \namespace useful Miscelaneous collection of useful objects and algorithms.
*/
namespace useful {
// Adapted from Howard Hinnant's answer here:
// https://stackoverflow.com/questions/22590821/convert-stdduration-to-human-readable-time
/** \brief Display execution time in human-readable format. */
inline std::ostream &display_duration(std::ostream &stream,
                                      std::chrono::nanoseconds ns) {
  io::StreamScopeFormat guard{stream};
  using days = std::chrono::duration<std::int_least64_t, std::ratio<86400>>;
  stream.fill('0');
  stream << std::right;
  auto dd = std::chrono::duration_cast<days>(ns);
  ns -= dd;
  auto hh = std::chrono::duration_cast<std::chrono::hours>(ns);
  ns -= hh;
  auto mm = std::chrono::duration_cast<std::chrono::minutes>(ns);
  ns -= mm;
  auto ss = std::chrono::duration_cast<std::chrono::seconds>(ns);
  ns -= ss;
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(ns);
  stream << std::setw(2) << dd.count() << "d:" << std::setw(2) << hh.count()
         << "h:" << std::setw(2) << mm.count() << "m:" << std::setw(2)
         << ss.count() << "s:" << std::setw(3) << ms.count() << "ms";

  return stream;
};

/** \brief Display execution time in human-readable format. */
template <typename Clock>
std::ostream &display_duration(std::ostream &stream,
                               std::chrono::time_point<Clock> start_time,
                               std::chrono::time_point<Clock> end_time) {
  return display_duration(stream,
                          std::chrono::duration_cast<std::chrono::nanoseconds>(
                              end_time - start_time));
}

/**
   \brief Check if sorted \p container contains \p val.
   \note: Container must be sorted.
*/
template <typename T, typename U>
bool contains(const std::vector<T> &container, U const &val) {
  auto it =
      std::lower_bound(container.begin(), container.end(), val,
                       [](T const &elem, U const &val) { return elem < val; });

  return it != container.end() && *it == val;
}

/**
   \brief Check whether \p container contains \p value.
   \details Container must be sorted according to \p comp_less, equality is
   checked using \p comp_eq.
*/
template <typename T, typename U, typename Comp_less, typename Comp_eq>
bool contains(std::vector<T> const &container, U const &val,
              Comp_less comp_less = std::less<void>{},
              Comp_eq comp_eq = std::equal_to<void>{}) {
  auto it =
      std::lower_bound(container.begin(), container.end(), val, comp_less);

  return it != container.end() && comp_eq(*it, val);
}

/** \brief Check if \p string ends with \p suffix. */
inline bool ends_with(const std::string &string, const std::string &suffix) {
  return string.size() >= suffix.size() &&
         string.substr(string.size() - suffix.size()) == suffix;
}

/**
   \class StoreConst General/Useful.h "General/Useful.h"
   \brief Class holding const object.
*/
template <typename Object_Type, typename Return_Type = Object_Type>
struct StoreConst {
  using value_type = Return_Type;

  Return_Type operator()() const { return obj; }
  const Object_Type obj;
};

/**
   \class Store General/Useful.h "General/Useful.h"
   \brief Class holding object.
*/
template <typename Object_Type, typename Return_Type = Object_Type>
struct Store {
  using value_type = Return_Type;

  Return_Type operator()() const { return obj; }
  Object_Type obj;
};

/**
   \brief Make an object and initialize it given \p parameters.
   \note Object must implement:
   - initialize(Params const&).
*/
template <typename Object, typename Params>
Object create(Params const &parameters = meta::Empty{}) {
  Object object;
  object.initialize(parameters);

  return object;
}

/**
   \class Maker General/Useful.h "General/Useful.h"
   \brief Functor to make objects of given type, passing arbitrary parameters.
*/
template <typename Object> struct Maker {
  template <typename... Args> Object operator()(Args... args) const {
    return Object{args...};
  }
};

/**
   \class Creator General/Useful.h "General/Useful.h"
   \brief Functor to make an object on the heap, passing arbitrary parameters.
*/
template <typename T> struct Creator {
  using value_type = T;
  using pointer = T *;

  template <typename... Args> pointer operator()(Args... args) const {
    return new T{args...};
  }
};

/**
   \class Forward General/Useful.h "General/Useful.h"
   \brief Functor to forward an argument.
*/
template <typename Object_Type> struct Forward {
  Object_Type operator()(Object_Type const &object) { return object; }
};

/**
   \class Forward_ref General/Useful.h "General/Useful.h"
   \brief Functor to forward a const reference.
*/
template <typename Object_Type> struct Forward_ref {
  Object_Type const &operator()(Object_Type const &object) { return object; }
};

/**
   \class hash_container General/Useful.h "General/Useful.h"
   \brief Simple hash for a container by combining element hashes.
*/
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

/** \brief Combine \p seed with the hash of an object \p v. */
template <typename T> void hash_combine(std::size_t &seed, T const &vv) {
  std::hash<T> hasher;
  seed ^= hasher(vv) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

/**
   \class Hash_pair General/Useful.h "General/Useful.h"
   \brief Hash for std::pair.
*/
template <typename S, typename T> struct Hash_pair {
  std::size_t operator()(std::pair<S, T> const &vv) const {
    std::size_t seed = 0;
    hash_combine(seed, vv.first);
    hash_combine(seed, vv.second);
    return seed;
  }
};

/** \return Sign of \p val. */
template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

/** \brief If \p to_replace is NaN, replace it with \p replace_with. */
template <typename T> T &deNaN(T &to_replace, T replace_with = T()) {
  if constexpr (meta::has_begin_v<T>) {
    for (auto const &val : to_replace)
      denan(val, replace_with);
  } else {
    if (to_replace != to_replace)
      to_replace = replace_with;
  }

  return to_replace;
}

/**
   \brief If \p to_chop is smaller than \p tolerance, replace it with \p
   replace_with.
*/
template <typename T> T &chop(T &to_chop, T tolerance, T replace_with = T()) {
  if constexpr (meta::has_begin_v<T>) {
    for (auto const &val : to_chop)
      denan(val, replace_with);
  } else {
    if (std::abs(to_chop - replace_with) < tolerance)
      to_chop = replace_with;
  }

  return to_chop;
}

/**
   \brief Copy end element of \p container to index \p idx, then erase end
   element.
*/
template <typename Container>
Container &swap_erase(Container &container, std::size_t idx) {
  if (idx != container.size() - 1)
    container[idx] = container.back();
  container.pop_back();

  return container;
}

/**
   \brief Swap-erase elements in \p container at indices \p indices.
   \note \p positions will be sorted, and must implement <tt>
   %sort(std::greater<std::size_t>)</tt>.
*/
template <typename Container, typename List>
Container &swap_erase(Container &container, List &&indices) {
  indices.sort(std::greater<std::size_t>{});
  for (auto const position : indices)
    swap_erase(container, indices);

  return container;
}

/** \brief Swap-erase elements of \p container satisfying \p criterion. */
template <typename Container, typename Criterion>
Container &swap_erase_if(Container &container, Criterion &&criterion) {
  for (std::size_t ii = container.size(); ii-- > 0;)
    if (criterion(container[ii]))
      swap_erase(container, ii);

  return container;
}

/**
   \brief Like swap_erase(), but call \c delete on the
   removed element.
*/
template <typename Container>
Container &swap_delete(Container &container, std::size_t idx) {
  if (idx != container.size() - 1)
    container[idx] = container.back();
  delete container.back();
  container.pop_back();

  return container;
}

/** \brief Swap-delete elements in \p container at positions \p positions. */
template <typename Container, typename List>
Container &swap_delete(Container &container, List &indices) {
  indices.sort(std::greater<std::size_t>{});
  for (auto const position : indices)
    swap_delete(container, position);

  return container;
}

/** \brief Swap-delete elements of \p container satisfying \p criterion. */
template <typename Container, typename Criterion>
Container &swap_delete_if(Container &container, Criterion &&criterion) {
  for (std::size_t ii = container.size(); ii-- > 0;)
    if (criterion(container[ii]))
      swap_delete(container, ii);

  return container;
}
} // namespace useful

#endif /* GENERAL_USEFUL_H */
