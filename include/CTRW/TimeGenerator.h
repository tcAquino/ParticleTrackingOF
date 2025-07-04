/**
   \file CTRW/TimeGenerator.h
   \author Tomas Aquino
   \date 10/23/2019

   \brief Objects to generate time increments for particle tracking.

   \details
   A TimeGenerator should implement the following minimal functionality:

   \code{.cpp}
   class TimeGenerator
   {
   template <typename State>
   val_type operator() (State const&)
   {
   // Return time increment, with val_type
   // a scalar type (e.g. double)
   }
   };
   \endcode
*/

#ifndef CTRW_TIMEGENERATOR_H
#define CTRW_TIMEGENERATOR_H

#include "General/Meta.h"
#include "General/Parallel.h"
#include "Stochastic/Random.h"
#include <random>

namespace ctrw {
/**
   \class TimeGenerator_Step CTRW/TimeGenerator.h "CTRW/TimeGenerator.h"
   \brief Deterministic time step.
*/
template <typename val_type = double> class TimeGenerator_Step {
  val_type _time_step;

public:
  using value_type = val_type;

  TimeGenerator_Step(value_type time_step = 0.) : _time_step{time_step} {}

  void time_step(value_type time_step) { _time_step = time_step; }

  value_type time_step() const { return _time_step; }

  template <typename State = meta::Empty>
  value_type operator()(State const & = {}) {
    return _time_step;
  }
};

/**
   \class TimeGenerator_Dist CTRW/TimeGenerator.h "CTRW/TimeGenerator.h"
   \brief Time step according to a given distribution.
*/
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
