//
//  TimeGenerator.h
//  CTRW
//
//  Created by Tomas Aquino on 10/23/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

#ifndef TimeGenerator_h
#define TimeGenerator_h

#include <random>

namespace ctrw
{
  // A TimeGenerator should implement the following
  // minimum functionality:
  // class TimeGenerator
  // {
  //   template <typename State>
  //   val_type operator() (State const&)
  //   {
  //     // Return time increment, with val_type
  //     // a scalar type (e.g. double)
  //   }
  // }
  
  // Deterministic time step
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
  
  // Time step according to a given distribution
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



#endif /* TimeGenerator_h */
