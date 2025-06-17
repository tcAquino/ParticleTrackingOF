/**
   \file General/Parallel.h
   \author Tomás Aquino
   \date 14/05/2024
   \brief Definitions and algorithms for serial and parallel implementations.
*/

#ifndef GENERAL_PARALLEL_H
#define GENERAL_PARALLEL_H

#include <cstddef>
#include <vector>

#if __has_include(<omp.h>)
#include <omp.h>
#endif /* __has_include(<omp.h>) */

/**
   \namespace par Definitions and algorithms for serial and parallel
   implementations.
*/
namespace par {
/**
   \struct ParallelOptions General/Parallel.h "General/Parallel.h"
   \brief Options to choose between serial and parallel implementations.
*/
struct ParallelOptions {
  /**
     \struct ParallelOptions::Serial General/Parallel.h "General/Parallel.h"
     \brief Serial implementations.
  */
  struct Serial {};

  /**
     \struct ParallelOptions::Parallel General/Parallel.h "General/Parallel.h"
     \brief Parallel implementations.
  */
  struct Parallel {};
};

/**
   \brief Get number of current parallel thread.
   \details Passing an empty object of type ParallelOptions::Serial
   chooses this implementation (serial).
   \return Current thread number, always \c 0.
*/
inline std::size_t get_thread_num(ParallelOptions::Serial) { return 0; }

/**
   \brief Get total number of parallel threads.
   \details Passing an empty object of ParallelOptions::Serial
   chooses this implementation (serial).
   \return Number of threads, always \c 1.
*/
inline std::size_t get_num_threads(ParallelOptions::Serial) { return 1; }

/**
   \brief Get number of current parallel thread.
   \details Passing an empty object of type ParallelOptions::Parallel
   chooses this implementation (parallel).
   \return Current thread number.
*/
inline std::size_t get_thread_num(ParallelOptions::Parallel) {
  return omp_get_thread_num();
}

/**
   \brief Get total number of parallel threads.
   \details Passing an empty object of type ParallelOptions::Parallel
   chooses this implementation (parallel).
   \return Number of threads.
*/
inline std::size_t get_num_threads(ParallelOptions::Parallel) {
  std::size_t num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
#endif /** _OPENMP */
  num_threads = omp_get_num_threads();
  return num_threads;
}

/**
   \class ThreatedType General/Parallel.h "General/Parallel.h"
   \brief Set of values of type \c Type for separate thread use.
*/
template <typename Type, typename ParallelOption> struct Threaded {
  Threaded(Type val) : _vals(get_num_threads(ParallelOption{}), val) {}

  Threaded(std::vector<Type> &&vals)
      : _vals{std::forward<std::vector<Type>>(vals)} {}

  Threaded() : _vals(get_num_threads(ParallelOption{})) {}

  Type &operator()() { return _vals[get_thread_num(ParallelOption{})]; }

private:
  std::vector<Type> _vals;
};

template <typename Type> struct Threaded<Type, ParallelOptions::Serial> {
  Threaded(Type val)
      : _val{val}
  {}

  Threaded(std::vector<Type> const &vals) : _val{vals[0]} {}

  Threaded() : _val{} {}

  Type &operator()() { return _val; }

private:
  Type _val;
};
} // namespace par

#endif /* GENERAL_PARALLEL_H */
