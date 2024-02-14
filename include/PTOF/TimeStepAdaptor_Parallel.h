/**
* \file PTOF/TimeStepAdaptor_Parallel.h
* \author Tomás Aquino
* \date 13/02/2024
*/

#ifndef PTOF_TIMESTEPADAPTOR_PARALLEL_H
#define PTOF_TIMESTEPADAPTOR_PARALLEL_H

#include <algorithm>
#include <fieldTypes.H>
#include <MinMax.H>
#include <omp.h>
#include <zero.H>
#include "PTOF/Field.h"
#include "PTOF/Reaction.h"

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
  typename ReactionParameters,
  typename SolverParameters,
  bool check_if_outside,
  bool warn_if_outside>
  class TimeStepAdaptor_CellSize_SurfaceReaction_Parallel
  {
  public:
    TimeStepAdaptor_CellSize_SurfaceReaction_Parallel
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     SurfaceReaction& reaction,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     FieldOptions::NoCheck)
    : TimeStepAdaptor_CellSize_SurfaceReaction_Parallel{
      geometry, velocity_field, reaction,
      params_transport, params_reaction, params_solvers }
    {}
    
    TimeStepAdaptor_CellSize_SurfaceReaction_Parallel
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     SurfaceReaction& reaction,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     FieldOptions::Check)
    : TimeStepAdaptor_CellSize_SurfaceReaction_Parallel{
      geometry, velocity_field, reaction,
      params_transport, params_reaction, params_solvers }
    {}
    
    TimeStepAdaptor_CellSize_SurfaceReaction_Parallel
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     SurfaceReaction& reaction,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     FieldOptions::Warn)
    : TimeStepAdaptor_CellSize_SurfaceReaction_Parallel{
      geometry, velocity_field, reaction,
      params_transport, params_reaction, params_solvers }
    {}
    
    template <typename State, typename TimeGenerator, typename JumpGenerator>
    void operator()
    (State& state, TimeGenerator& time_generator, JumpGenerator& jump_generator)
    {
      state.cell = geometry.locator[omp_get_thread_num()](state);
      auto cell_id = state.cell;
      if constexpr (check_if_outside)
        if (outside(make_point(state.position), state.cell))
          cell_id = geometry.locator[omp_get_thread_num()].nearest_cell(state);
      double cell_side = std::pow(geometry.mesh.V()[cell_id]/volume_factor,
                                  1./Geometry::dim);
      double local_advection_time = cell_side/Foam::mag(velocity_field(state));
      double local_diffusion_time = cell_side*cell_side/(2.*params_transport.diff_coeff);
      
      double reaction_rate = 0.;
      if constexpr (!std::is_same_v<SurfaceReaction, Reaction_DoNothing>)
        for (auto face : geometry.mesh.cells()[cell_id])
        {
          double rate = reaction.rate(face);
          if (rate > reaction_rate)
            reaction_rate = rate;
        }
      double local_reaction_time = params_transport.diff_coeff/(reaction_rate*reaction_rate);
      
      timestep = std::max(std::min({
          params_solvers.local_time_step_adv*local_advection_time,
          params_solvers.local_time_step_diff*local_diffusion_time,
          params_solvers.local_time_step_react*local_reaction_time }),
                          std::min({
          params_solvers.global_time_step_adv*params_transport.advection_time,
          params_solvers.global_time_step_adv*params_transport.diffusion_time,
          params_solvers.global_time_step_adv*params_reaction.reaction_time }));
      
      if constexpr (has_time_step_setter<TimeGenerator>::value)
        time_generator.time_step(timestep);
      if constexpr (has_time_step_setter<JumpGenerator>::value)
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
    ReactionParameters const& params_reaction;
    SolverParameters const& params_solvers;
    double volume_factor;
    double timestep;
    
    TimeStepAdaptor_CellSize_SurfaceReaction_Parallel
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     SurfaceReaction& reaction,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers)
    : geometry{ geometry }
    , velocity_field{ velocity_field }
    , reaction{ reaction }
    , params_transport{ params_transport }
    , params_reaction{ params_reaction }
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
    
    /** Check if type T has member function void time_step(double)
     Adapted from kispaljr's answer here: https://stackoverflow.com/questions/257288/templated-check-for-the-existence-of-a-class-member-function */
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
    
    /** Bounds checking. */
    bool outside
    (Foam::point const& position, Foam::label cell) const
    {
      if (cell == -1)
      {
        if constexpr (warn_if_outside)
          std::cerr << "Warning: Requested cell side at position "
                    << "("
                    << position[0] << ", "
                    << position[1] << ", "
                    << position[2] << ")"
                    << " outside mesh. Using nearest cell\n";
        return 1;
      }
      return 0;
    }
  };
  template
  <typename Geometry,
  typename VelocityField,
  typename SurfaceReaction,
  typename TransportParameters,
  typename ReactionParameters,
  typename SolverParameters>
  TimeStepAdaptor_CellSize_SurfaceReaction_Parallel
  (Geometry const&, VelocityField const&, SurfaceReaction&,
   TransportParameters const&, ReactionParameters const&, SolverParameters const&,
   FieldOptions::NoCheck) ->
  TimeStepAdaptor_CellSize_SurfaceReaction_Parallel
  <Geometry, VelocityField, SurfaceReaction,
  TransportParameters, ReactionParameters, SolverParameters, 0, 0>;
  template
  <typename Geometry,
  typename VelocityField,
  typename SurfaceReaction,
  typename TransportParameters,
  typename ReactionParameters,
  typename SolverParameters>
  TimeStepAdaptor_CellSize_SurfaceReaction_Parallel
  (Geometry const&, VelocityField const&, SurfaceReaction&,
   TransportParameters const&, ReactionParameters const&, SolverParameters const&,
   FieldOptions::Check) ->
  TimeStepAdaptor_CellSize_SurfaceReaction_Parallel
  <Geometry, VelocityField, SurfaceReaction,
  TransportParameters, ReactionParameters, SolverParameters, 1, 0>;
  template
  <typename Geometry,
  typename VelocityField,
  typename SurfaceReaction,
  typename TransportParameters,
  typename ReactionParameters,
  typename SolverParameters>
  TimeStepAdaptor_CellSize_SurfaceReaction_Parallel
  (Geometry const&, VelocityField const&, SurfaceReaction&,
   TransportParameters const&, ReactionParameters const&, SolverParameters const&,
   FieldOptions::Warn) ->
  TimeStepAdaptor_CellSize_SurfaceReaction_Parallel
  <Geometry, VelocityField, SurfaceReaction,
  TransportParameters, ReactionParameters, SolverParameters, 1, 1>;
}

#endif /* PTOF_TIMESTEPADAPTOR_PARALLEL_H */
