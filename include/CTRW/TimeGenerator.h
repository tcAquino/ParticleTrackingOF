/**
 * @file   TimeGenerator.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Sun Mar  8 12:45:37 2026
 *
 * @brief Objects to generate time increments for CTRW.
 *
 * @details A \c TimeGenerator should implement the following functionality:
 *
 * @code{.cpp}
 * class TimeGenerator {
 * public:
 *   template <typename State>
 *   auto operator() (State const & state) {
 *   // Return time increment of scalar type (typically double).
 *   }
 * };
 * @endcode
 */

#ifndef CTRW_TIMEGENERATOR_H
#define CTRW_TIMEGENERATOR_H

#include "General/Meta.h"
#include "General/Parallel.h"
#include "Stochastic/Random.h"
#include <random>

namespace ctrw {
/** @brief Deterministic time step. */
template <typename ParallelOption = par::ParallelOptions::Serial,
          typename val_type = double>
class TimeGenerator_Step {
  par::Threaded<double, ParallelOption> _time_step; /**< Time step. */

public:
  using value_type = val_type;

  TimeGenerator_Step(value_type time_step = 0.) : _time_step{time_step} {}

  void time_step(value_type time_step) { _time_step() = time_step; }

  value_type time_step() const { return _time_step(); }

  template <typename State = meta::Empty>
  value_type operator()(State const & = {}) {
    return _time_step();
  }
};

/** @brief Time step according to a given distribution. */
template <typename Dist, typename ParallelOption = par::ParallelOptions::Serial,
          typename RNGEngine = std::mt19937>
class TimeGenerator_Dist {
  Dist _dist;
  stochastic::RNGThreaded<ParallelOption, RNGEngine> _rng{
      std::random_device{}()};

public:
  using value_type = decltype(_dist(_rng));

  TimeGenerator_Dist(Dist dist) : _dist{dist} {}

  TimeGenerator_Dist(Dist dist, ParallelOption) : TimeGenerator_Dist{dist} {}

  template <typename State = meta::Empty> auto operator()(State const & = {}) {
    return _dist(_rng);
  }
};
} // namespace ctrw

#endif /* CTRW_TIMEGENERATOR_H */
