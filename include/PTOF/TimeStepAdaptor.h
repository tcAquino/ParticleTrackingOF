/**
* \file PTOF/TimeStepAdaptor.h
* \author Tomás Aquino
* \date 09/26/2017
*/

#ifndef PTOF_TIMESTEPADAPTOR_H
#define PTOF_TIMESTEPADAPTOR_H

#include <algorithm>
#include <MinMax.H>
#include <zero.H>

namespace ptof
{
  /** \struct TimeStepAdaptor_CellSize PTOF/TimeStepAdaptor.h "PTOF/TimeStepAdaptor.h"
   * \brief Adaptive time step control based on local cell size.
   * Takes reaction limitations to time step only if close to reactive boundary */
  template
  <typename Geometry,
  typename VelocityField,
  typename SurfaceReaction,
  typename TransportParameters,
  typename SolverParameters>
  class TimeStepAdaptor_CellSize_SurfaceReaction
  {
  public:
    TimeStepAdaptor_CellSize_SurfaceReaction
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     SurfaceReaction& reaction,
     TransportParameters const& params_transport,
     SolverParameters const& params_solvers)
    : geometry{ geometry }
    , velocity_field{ velocity_field }
    , reaction{ reaction }
    , params_transport{ params_transport }
    , params_solvers{ params_solvers }
    , volume_factor{ 1. }
    {
      for (std::size_t dd = 2; dd --> Geometry::dim;)
      {
        Foam::vector direction = Foam::zero{};
        direction[dd] = 1.;
        std::vector<double> max;
        std::vector<double> min;
        for (auto const& patch : geometry.mesh.boundary())
        {
          max.push_back(Foam::max(patch.Cf() & direction));
          min.push_back(Foam::min(patch.Cf() & direction));
        }
        volume_factor *= *std::max_element(max.begin(), max.end()) -
        *std::min_element(min.begin(), min.end());
      }
    }
    
    template <typename State, typename TimeGenerator, typename JumpGenerator>
    void operator()
    (State& state, TimeGenerator& time_generator, JumpGenerator& jump_generator)
    {
      state.cell = geometry.locator(state);
      double cell_side = std::pow(geometry.mesh.V()[state.cell]/volume_factor,
                                  1./Geometry::dim);
      double advection_time = params_solvers.time_step_accuracy_adv*
        cell_side/Foam::mag(velocity_field(state));
      double diffusion_time = params_solvers.time_step_accuracy_diff*
        cell_side*cell_side/(2.*params_transport.diff_coeff);
      
      double reaction_rate = 0.;
      for (auto face : geometry.mesh.cells()[state.cell])
      {
        double rate = reaction.rate(face);
        if (rate > reaction_rate)
          reaction_rate = rate;
      }
      double reaction_time = params_solvers.time_step_accuracy_react*
        params_transport.diff_coeff/
        (reaction_rate*reaction_rate);
      
      timestep = std::min({ advection_time,
        diffusion_time,
        reaction_time });
      
      if constexpr (has_time_step_setter<TimeGenerator>::value)
        time_generator.time_step(timestep);
      if constexpr (has_time_step<JumpGenerator>::value)
        jump_generator.time_step(timestep);
      if constexpr (has_time_step_setter<SurfaceReaction>::value)
        reaction.time_step(timestep);
    }
    
    double time_step()
    { return timestep; }
    
  private:
    Geometry const& geometry;
    VelocityField const& velocity_field;
    SurfaceReaction& reaction;
    TransportParameters const& params_transport;
    SolverParameters const& params_solvers;
    double volume_factor;
    double timestep;
    double diff_coeff;
    
    // Check if type T has member function void time_step(double)
    // Adapted from kispaljr's answer here:
    // https://stackoverflow.com/questions/257288/templated-check-for-the-existence-of-a-class-member-function
    template <typename T> struct has_time_step_setter
    {
        typedef char (&Yes)[1];
        typedef char (&No)[2];

        template<class U>
        static Yes test(U* data,
                        typename std::enable_if<std::is_void<
                          decltype(data->time_step(0.))>::value>::type* = 0);
        static No test(...);
        static const bool value =
          sizeof(Yes) ==
            sizeof(has_time_step_setter::test((typename std::remove_reference<T>::type*)0));
    };
  };
}

#endif /* PTOF_TIMESTEPADAPTOR_H */
