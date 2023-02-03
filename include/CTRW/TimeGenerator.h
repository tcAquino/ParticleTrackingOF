/**
 * \file CTRW/TimeGenerator.h
 * \author Tomás Aquino
 * \date 10/23/2019
*/

#ifndef CTRW_TIMEGENERATOR_H
#define CTRW_TIMEGENERATOR_H

#include <random>

namespace ctrw
{
  // A TimeGenerator should implement the following minimal functionality:
  //
  // class TimeGenerator
  // {
  //   template <typename State>
  //   val_type operator() (State const&)
  //   {
  //     // Return time increment, with val_type
  //     // a scalar type (e.g. double)
  //   }
  // };
  
  /** \class TimeGenerator_Step CTRW/TimeGenerator.h "CTRW/TimeGenerator.h"
   *  \brief Deterministic time step. */
  template <typename val_type = double>
  class TimeGenerator_Step
  {
    val_type dt;

  public:
    using value_type = val_type;

    TimeGenerator_Step(val_type dt = 0.)
    : dt{ dt }
    {}

    void time_step(val_type dt)
    { this->dt = dt; }

    val_type time_step() const
    { return dt; }

    template <typename State = useful::Empty>
    val_type operator()(State const& = {})
    { return dt; }
  };
  
  /** \class TimeGenerator_Dist CTRW/TimeGenerator.h "CTRW/TimeGenerator.h"
   * \brief Time step according to a given distribution. */
  template <typename Dist, typename RNG = std::mt19937>
  class TimeGenerator_Dist
  {
    Dist dist;
    RNG rng{ std::random_device{}() };

  public:
    using value_type = decltype(dist(rng));

    TimeGenerator_Dist(Dist dist)
    : dist{ dist }
    {}

    template <typename State = useful::Empty>
    auto operator()(State const& = {})
    { return dist(rng); }
  };
}



#endif /* CTRW_TIMEGENERATOR_H */
