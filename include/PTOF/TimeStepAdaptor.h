/**
 \file PTOF/TimeStepAdaptor.h
 \author Tomás Aquino
 \date 09/26/2017
*/

#ifndef PTOF_TIMESTEPADAPTOR_H
#define PTOF_TIMESTEPADAPTOR_H

#include <algorithm>
#include <fieldTypes.H>
#include <MinMax.H>
#include <zero.H>
#include "PTOF/Field.h"
#include "PTOF/Reaction.h"

namespace ptof
{

  /** \class TimeStepAdaptor_CellSize PTOF/TimeStepAdaptor.h "PTOF/TimeStepAdaptor.h"
    \brief Adaptive time step control based on local cell size.
    \details Takes reaction limitations to time step only if close to reactive boundary */
  template
  <typename Geometry,
  typename VelocityField,
  typename SurfaceReaction,
  typename TransportParameters,
  typename ReactionParameters,
  typename SolverParameters,
  bool check_if_outside,
  bool warn_if_outside>
  class TimeStepAdaptor_CellSize_SurfaceReaction
  {
  public:
    /** Constructor.
     \brief Without checking if particle is in bounds.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field.
     \param reaction Surface reaction handler.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
    */
    TimeStepAdaptor_CellSize_SurfaceReaction
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     SurfaceReaction& reaction,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     CheckOptions::NoCheck)
    : TimeStepAdaptor_CellSize_SurfaceReaction{
        geometry, velocity_field, reaction,
        params_transport, params_reaction, params_solvers }
    {
      static_assert(check_if_outside == false && warn_if_outside == false,
                    "Bad template arguments for no bounds checking.");
    }
    
    /** Constructor.
     \brief With checking if particle is in bounds but no warning if not.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field.
     \param reaction Surface reaction handler.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
    */
    TimeStepAdaptor_CellSize_SurfaceReaction
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     SurfaceReaction& reaction,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     CheckOptions::Check)
    : TimeStepAdaptor_CellSize_SurfaceReaction{
        geometry, velocity_field, reaction,
        params_transport, params_reaction, params_solvers }
    {
      static_assert(check_if_outside == true && warn_if_outside == false,
                    "Bad template argument for bounds checking.");
    }
    
    /** Constructor.
     \brief With warning if particle is out of bounds.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field.
     \param reaction Surface reaction handler.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
    */
    TimeStepAdaptor_CellSize_SurfaceReaction
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     SurfaceReaction& reaction,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     CheckOptions::Warn)
    : TimeStepAdaptor_CellSize_SurfaceReaction{
        geometry, velocity_field, reaction,
        params_transport, params_reaction, params_solvers }
    {
      static_assert(check_if_outside == true && warn_if_outside == true,
      "Bad template argument for bounds warning.");
    }
    
    /**
     \brief Update time step.
     \param state Particle state.
     \param time_generator Time step generator object to be updated.
     \param jump_generator Jump generator object to be updated.
     \return New time step.
    */
    template <typename State, typename TimeGenerator, typename JumpGenerator>
    double operator()
    (State& state, TimeGenerator& time_generator, JumpGenerator& jump_generator) const
    {
      state.cell = _geometry.locator(state);
      auto cell_id = state.cell;
      if constexpr (check_if_outside)
        if (outside(make_point(state.position), state.cell))
        {
          cell_id = _geometry.locator.nearest_cell(state);
        }
      double cell_side = std::pow(_geometry.mesh().V()[cell_id]/_volume_factor,
                                  1./Geometry::dim);
      double local_advection_time = cell_side/Foam::mag(_velocity_field(state));
      double local_diffusion_time = cell_side*cell_side/(2.*_params_transport.diff_coeff);
      
      double reaction_rate = 0.;
      if constexpr (!std::is_same_v<SurfaceReaction, Reaction_DoNothing>)
        for (auto face : _geometry.mesh().cells()[cell_id])
        {
          double rate = _reaction.rate(face);
          if (rate > reaction_rate)
            reaction_rate = rate;
        }
      double local_reaction_time = _params_transport.diff_coeff/(reaction_rate*reaction_rate);
      
      double time_step = std::max(std::min({
        _params_solvers.local_time_step_adv*local_advection_time,
        _params_solvers.local_time_step_diff*local_diffusion_time,
        _params_solvers.local_time_step_react*local_reaction_time }),
                                  std::min({
        _params_solvers.global_time_step_adv*_params_transport.advection_time,
        _params_solvers.global_time_step_adv*_params_transport.diffusion_time,
        _params_solvers.global_time_step_adv*_params_reaction.reaction_time }));
      
      if constexpr (useful::has_time_step_setter<TimeGenerator>::value)
        time_generator.time_step(time_step);
      if constexpr (useful::has_time_step_setter<JumpGenerator>::value)
        jump_generator.time_step(time_step);
      if constexpr (useful::has_time_step_setter<SurfaceReaction>::value)
        _reaction.time_step(time_step);
      
      return time_step;
    }
    
  private:
    Geometry const& _geometry;                     /**< Domain geometry info and utilities. */
    VelocityField const& _velocity_field;          /**< Velocity field. */
    SurfaceReaction& _reaction;                    /**< Surface reaction. */
    TransportParameters const& _params_transport;  /**< Transport parameters. */
    ReactionParameters const& _params_reaction;    /**< Reaction parameters. */
    SolverParameters const& _params_solvers;       /**< Solver parameters. */
    double _volume_factor;                         /**< Volume factor to correct cell volumes in 1D and 2D. */
    
    /** Base constructor.
     \brief Called by different public constructors.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field.
     \param reaction Surface reaction handler.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
    */
    TimeStepAdaptor_CellSize_SurfaceReaction
    (Geometry const& geometry,
     VelocityField const& velocity_field,
     SurfaceReaction& reaction,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers)
    : _geometry{ geometry }
    , _velocity_field{ velocity_field }
    , _reaction{ reaction }
    , _params_transport{ params_transport }
    , _params_reaction{ params_reaction }
    , _params_solvers{ params_solvers }
    , _volume_factor{ 1. }
    {
      for (std::size_t dd = 2; dd --> Geometry::dim;)
      {
        Foam::vector direction = Foam::zero{};
        direction[dd] = 1.;
        std::vector<double> max;
        std::vector<double> min;
        for (auto const& patch : _geometry.mesh().boundary())
        {
          max.push_back(Foam::max(patch.Cf() & direction));
          min.push_back(Foam::min(patch.Cf() & direction));
        }
        _volume_factor *= *std::max_element(max.begin(), max.end()) -
          *std::min_element(min.begin(), min.end());
      }
    }
    
    /** Check if position is out of bounds
    \param position Position to check.
    \param cell Hint of which mesh cell position is in.
    \return true if outside, false otherwise.
    */
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
  TimeStepAdaptor_CellSize_SurfaceReaction
  (Geometry const&, VelocityField const&, SurfaceReaction&,
   TransportParameters const&, ReactionParameters const&, SolverParameters const&,
   CheckOptions::NoCheck) ->
  TimeStepAdaptor_CellSize_SurfaceReaction
  <Geometry, VelocityField, SurfaceReaction,
  TransportParameters, ReactionParameters, SolverParameters, 0, 0>;
  template
  <typename Geometry,
  typename VelocityField,
  typename SurfaceReaction,
  typename TransportParameters,
  typename ReactionParameters,
  typename SolverParameters>
  TimeStepAdaptor_CellSize_SurfaceReaction
  (Geometry const&, VelocityField const&, SurfaceReaction&,
   TransportParameters const&, ReactionParameters const&, SolverParameters const&,
   CheckOptions::Check) ->
  TimeStepAdaptor_CellSize_SurfaceReaction
  <Geometry, VelocityField, SurfaceReaction,
  TransportParameters, ReactionParameters, SolverParameters, 1, 0>;
  template
  <typename Geometry,
  typename VelocityField,
  typename SurfaceReaction,
  typename TransportParameters,
  typename ReactionParameters,
  typename SolverParameters>
  TimeStepAdaptor_CellSize_SurfaceReaction
  (Geometry const&, VelocityField const&, SurfaceReaction&,
   TransportParameters const&, ReactionParameters const&, SolverParameters const&,
   CheckOptions::Warn) ->
  TimeStepAdaptor_CellSize_SurfaceReaction
  <Geometry, VelocityField, SurfaceReaction,
  TransportParameters, ReactionParameters, SolverParameters, 1, 1>;
}

#endif /* PTOF_TIMESTEPADAPTOR_H */
