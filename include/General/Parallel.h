/**
 \file General/Parallel.h
 \author Tomás Aquino
 \date 14/05/2024
 \brief Algorithms for parallel implementations
*/

#ifndef GENERAL_PARALLEL_H
#define GENERAL_PARALLEL_H

#include "General/Meta.h"
#include <omp.h>

namespace useful {
inline std::size_t get_thread_num(meta::ParallelOptions::Parallel) {
  return omp_get_thread_num();
}

inline std::size_t get_num_threads(meta::ParallelOptions::Parallel) {
  std::size_t num_threads = 1;
#pragma omp parallel
  { num_threads = omp_get_num_threads(); }
  return num_threads;
}
} // namespace useful

#endif /* GENERAL_PARALLEL_H */
