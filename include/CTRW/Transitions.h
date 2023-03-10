/**
* \file CTRW/Transitions.h
* \author Tomás Aquino
* \date 08/06/2018
*/

#ifndef CTRW_TRANSITIONS_H
#define CTRW_TRANSITIONS_H

#include <limits>
#include <utility>
#include "CTRW/JumpGenerator.h"
#include "CTRW/TimeGenerator.h"
#include "General/Operations.h"
#include "Geometry/Boundary.h"

namespace ctrw
{
  // Transitions objects should implement the following minimal functionality:
  // class Transitions
  // {
  // public:
  //   template <typename State>
  //   void operator() (State& state)
  //   {
  //     // Update state for a single transition
  //   }
  // };
  
  /** \class Transition_Time_Position CTRW/Transitions.h "CTRW/Transitions.h"
   * \brief Update state.position and state.time
   * by summing quantities returned by a
   * JumpGenerator and a TimeGenerator object given state.\n
   * 
   * Apply Boundary condition given new state and old state. */
	template <typename TimeGenerator, typename JumpGenerator,
  typename Boundary = geometry::Boundary_DoNothing>
	class Transitions_Time_Position
	{
	public:
		Transitions_Time_Position
    (TimeGenerator&& time_generator, JumpGenerator&& jump_generator,
     Boundary&& boundary = {})
    : time_generator{ std::forward<TimeGenerator>(time_generator) }
    , jump_generator{ std::forward<JumpGenerator>(jump_generator) }
    , boundary{ std::forward<Boundary>(boundary) }
		{}

		template <typename State>
		void operator() (State& state)
		{
      auto state_old = state;
      operation::plus_InPlace(state.position, jump_generator(state));
      state.time += time_generator(state);
      boundary(state, state_old);
		}
    
    auto const& get_time_generator() const
    { return time_generator; }
    
    auto const& get_jump_generator() const
    { return jump_generator; }
    
    auto const& get_boundary() const
    { return boundary; }

	private:
		TimeGenerator time_generator;
		JumpGenerator jump_generator;
    Boundary boundary;
	};
  template <typename TimeGenerator, typename JumpGenerator,
  typename Boundary>
  Transitions_Time_Position
  (TimeGenerator&&, JumpGenerator&&, Boundary&&) ->
  Transitions_Time_Position<TimeGenerator, JumpGenerator, Boundary>;
  
  /** \class Transitions_AdaptiveTimeStep_Time_Position CTRW/Transitions.h "CTRW/Transitions.h"
   *  \brief
   * Choose time step based on state
   * Update state.position and state.time
   * by summing quantities returned by a
   * JumpGenerator and a TimeGenerator object given state. 
   * 
   * Apply Boundary condition given new state and old state. \n
   * Locate all particles in cells before jumping. \n
   * State must define: 
   * - cell */
  template <typename TimeStepAdaptor, typename TimeGenerator, typename JumpGenerator,
  typename Boundary = geometry::Boundary_DoNothing>
  class Transitions_AdaptiveTimeStep_Time_Position
  {
  public:
    Transitions_AdaptiveTimeStep_Time_Position
    (TimeStepAdaptor&& time_step_adaptor,
     TimeGenerator&& time_generator,
     JumpGenerator&& jump_generator,
     Boundary&& boundary = {})
    : time_step_adaptor{ std::forward<TimeStepAdaptor>(time_step_adaptor) }
    , time_generator{ std::forward<TimeGenerator>(time_generator) }
    , jump_generator{ std::forward<JumpGenerator>(jump_generator) }
    , boundary{ std::forward<Boundary>(boundary) }
    {}

    template <typename State>
    void operator() (State& state)
    {
      auto state_old = state;
      time_step_adaptor(state, time_generator, jump_generator);
      operation::plus_InPlace(state.position, jump_generator(state));
      state.time += time_generator(state);
      boundary(state, state_old);
    }
    
    auto const& get_time_generator() const
    { return time_generator; }
    
    auto const& get_jump_generator() const
    { return jump_generator; }
    
    auto const& get_boundary() const
    { return boundary; }
    
    double time_step() const
    { return time_step_adaptor.time_step(); }

  private:
    TimeStepAdaptor time_step_adaptor;
    TimeGenerator time_generator;
    JumpGenerator jump_generator;
    Boundary boundary;
  };
  template <typename TimeStepAdaptor, typename TimeGenerator, typename JumpGenerator,
  typename Boundary>
  Transitions_AdaptiveTimeStep_Time_Position
  (TimeStepAdaptor&&, TimeGenerator&&, JumpGenerator&&, Boundary&&) ->
  Transitions_AdaptiveTimeStep_Time_Position
  <TimeStepAdaptor, TimeGenerator, JumpGenerator, Boundary>;
  
  /** \class Transitions_Positions CTRW/Transitions.h "CTRW/Transitions.h"
   *  \brief Update state.position
   *  by summing quantity returned by a
   * JumpGenerator object given state.
   *  
   * Apply Boundary condition given new state and old state. */
  template <typename JumpGenerator, typename Boundary = geometry::Boundary_DoNothing>
  class Transitions_Position
  {
  public:
    Transitions_Position
    (JumpGenerator&& jump_generator, Boundary&& boundary = {})
    : jump_generator{ std::forward<JumpGenerator>(jump_generator) }
    , boundary{ std::forward<Boundary>(boundary) }
    {}

    template <typename State>
    void operator()(State& state)
    {
      auto state_old = state;
      operation::plus_InPlace(state.position, jump_generator(state));
      boundary(state, state_old);
    }
    
    template <typename Time = double>
    void time_step(Time time_step)
    {
      jump_generator.time_step(time_step);
    }
    
    auto time_step() const
    {
      return jump_generator.time_step();
    }
    
    auto const& get_jump_generator() const
    { return jump_generator; }
    
    auto const& get_boundary() const
    { return boundary; }

  private:
    JumpGenerator jump_generator;
    Boundary boundary;
  };
  template <typename JumpGenerator, typename Boundary>
  Transitions_Position(JumpGenerator&&, Boundary&&) ->
  Transitions_Position<JumpGenerator, Boundary>;
  
  /** \struct VelocityFromGenerator
   *  \brief Helper class to return velocity from JumpGenerator object
   * given state. */
  template <typename JumpGenerator>
  struct VelocityFromGenerator
  {
    VelocityFromGenerator
    (JumpGenerator const& jump_generator)
    : jump_generator{ jump_generator }
    {}
    
    template <typename State>
    auto operator()(State const& state) const
    { return this->jump_generator.velocity(state); }
    
  private:
    JumpGenerator const& jump_generator;
  };
  
  /** \class Transitions_Position_VelocityStep CTRW/Transitions.h "CTRW/Transitions.h"
   * \brief Update state.position according to a JumpGenerator object given state
   * with time step according to prescribed jump size divided
   * by velocity calculated accorsing to Velocity object given state. 
   * 
   * Apply Boundary condition given new state and old state. */
  template <typename JumpGenerator,
  typename Velocity,
  typename Boundary = geometry::Boundary_DoNothing>
  class Transitions_Position_VelocityStep
  {
  public:
    Transitions_Position_VelocityStep
    (double step_length,
     JumpGenerator&& jump_generator,
     Velocity&& velocity,
     Boundary&& boundary)
    : step_length{ step_length }
    , jump_generator{ std::forward<JumpGenerator>(jump_generator) }
    , velocity_field{ std::forward<Velocity>(velocity_field) }
    , boundary{ std::forward<Boundary>(boundary) }
    {}
    
    Transitions_Position_VelocityStep
    (double step_length,
     JumpGenerator&& jump_generator,
     Boundary&& boundary = {})
    : step_length{ step_length }
    , jump_generator{ std::forward<JumpGenerator>(jump_generator) }
    , velocity_field{ VelocityFromGenerator{ this->jump_generator } }
    , boundary{ std::forward<Boundary>(boundary) }
    {}

    template <typename State>
    void operator() (State& state)
    {
      auto state_old = state;
      jump_generator.time_step(step_length/operation::abs(velocity(state)));
      if (jump_generator.time_step() == std::numeric_limits<double>::infinity())
      {
        state.time = std::numeric_limits<double>::infinity();
        return;
      }
      operation::plus_InPlace(state.position, jump_generator(state));
      state.time += jump_generator.time_step();
      boundary(state, state_old);
    }
    
    template <typename State>
    auto velocity(State const& state) const
    { return velocity_field(state); }
    
    auto time_step() const
    { return jump_generator.time_step(); }
    
    auto const& get_velocity_field() const
    { return velocity_field; }

    auto const& get_jump_generator() const
    { return jump_generator; }
    
    auto const& get_boundary() const
    { return boundary; }

  private:
    double step_length;
    JumpGenerator jump_generator;
    Velocity velocity_field;  /**< Velocity field given state. */
    Boundary boundary;
  };
  template <typename JumpGenerator, typename Velocity,
  typename Boundary>
  Transitions_Position_VelocityStep
  (double, JumpGenerator&&, Velocity&&, Boundary&&) ->
  Transitions_Position_VelocityStep<JumpGenerator, Velocity, Boundary>;
  template <typename JumpGenerator, typename Boundary>
  Transitions_Position_VelocityStep
  (double, JumpGenerator&&, Boundary&&) ->
  Transitions_Position_VelocityStep
  <JumpGenerator, VelocityFromGenerator<JumpGenerator>, Boundary>;

  /** \class Transitions_PTRW_FlowField_Diff CTRW/Transitions.h "CTRW/Transitions.h"
   *  \brief Update state.position according to forward-Euler
   * with prescribed time step,
   * with advection contribution according to flow field
   * along the 0 direction
   * given position and diffusion.
   * 
   * Apply Boundary condition given new state and old state. */
	template <typename FlowField, typename Boundary = geometry::Boundary_DoNothing>
	class Transitions_PTRW_FlowField_Diff
	{
  private:
    auto make_diff_generators
    (std::vector<double> const& diff, double timestep)
    {
      std::vector<JumpGenerator_Diffusion_1d> diff_generators;
      diff_generators.reserve(diff.size());
      for (auto const& diff_val : diff)
        diff_generators.emplace_back(diff_val, timestep);

      return diff_generators;
    }
    
	public:
    Transitions_PTRW_FlowField_Diff
    (std::vector<double> const& diff, double timestep,
     FlowField&& flow_field, Boundary&& boundary = {})
    : time_generator{ timestep }
    , diff_generators{ make_diff_generators(diff, timestep) }
    , flow_field{ std::forward<FlowField>(flow_field) }
    , boundary{ std::forward<Boundary>(boundary) }
		{}
    
		template <typename State>
		void operator()(State& state)
		{
      auto state_old = state;
      std::vector<double> jump(state.position.size());
      for (std::size_t dd = 0; dd < jump.size(); ++dd)
        jump[dd] = diff_generators[dd]();
      jump[0] += velocity(state)*time_generator.time_step();

      for (std::size_t dd = 0; dd < jump.size(); ++dd)
        state.position[dd] += jump[dd];

			boundary(state, state_old);
		}
    
    template <typename State>
    auto velocity(State const& state) const
    { return velocity_field(state.position); }

		void time_step(double timestep)
		{
      time_generator.time_step(timestep);
      for (auto& gen : diff_generators)
        gen.time_step(timestep);
    }

    void diff(std::size_t dd, double diff)
    { return diff_generators[dd].diff(diff); }

    double time_step() const
    { return time_generator.time_step(); }

    double diff(std::size_t dd) const
    { return diff_generators[dd].diff(); }
    
    auto const& get_time_generator() const
    { return time_generator; }
    
    auto const& get_boundary() const
    { return boundary; }
    
    auto const& get_velocity_field() const
    { return flow_field; }

	private:
    TimeGenerator_Step<double> time_generator;
    std::vector<JumpGenerator_Diffusion_1d> diff_generators;
    FlowField flow_field;
		Boundary boundary;
	};
  template <typename FlowField, typename Boundary>
  Transitions_PTRW_FlowField_Diff
  (std::vector<double> const&, double,
   FlowField&&, Boundary&&) ->
  Transitions_PTRW_FlowField_Diff<FlowField, Boundary>;

  /** \class Transitions_PTRW_Diffusion_1d CTRW/Transitions.h "CTRW/Transitions.h"
   *  \brief Update state.position according to forward-Euler
   * for diffusion in one dimension
   * with prescribed time step.
   * 
   * Apply Boundary condition given new state and old state. */
  template <typename Boundary = geometry::Boundary_DoNothing>
  class Transitions_PTRW_Diffusion_1d
  {
  public:
    Transitions_PTRW_Diffusion_1d
    (double diff, double timestep, Boundary&& boundary = {})
    : time_generator{ timestep }
    , jump_generator{ diff, timestep }
    , boundary{ std::forward<Boundary>(boundary) }
    {}

    template <typename State>
    void operator()(State& state)
    {
      auto state_old = state;
      state.position += jump_generator(state);
      boundary(state, state_old);
      state.time += time_step();
    }

    void time_step(double timestep)
    {
      time_generator.time_step(timestep);
      jump_generator.time_step(timestep);
    }

    void diff(double diff)
    { return jump_generator.diff(diff); }

    double time_step() const
    { return time_generator.time_step(); }

    double diff() const
    { return jump_generator.diff(); }
    
    auto const& get_time_generator() const
    { return time_generator; }
    
    auto const& get_jump_generator() const
    { return jump_generator; }
    
    auto const& get_boundary() const
    { return boundary; }

  private:
    TimeGenerator_Step<double> time_generator;
    JumpGenerator_Diffusion_1d jump_generator;
    Boundary boundary;
  };
  template <typename Boundary>
  Transitions_PTRW_Diffusion_1d
  (double, double, Boundary&&) ->
  Transitions_PTRW_Diffusion_1d<Boundary>;
  
  /** \class Transitions_PTRW_Transport_Reaction CTRW/Transitions.h "CTRW/Transitions.h"
   *  \brief Update state according to a transport transition given state
   * followed by a reaction transition given state and time step,
   * with a prescribed time step.
   * 
   * Note: transport rule should enforce any boundary conditions. */
  template <typename Transitions_Transport, typename Reaction>
  class Transitions_PTRW_Transport_Reaction
  {
  public:
    Transitions_PTRW_Transport_Reaction
    (Transitions_Transport&& transitions_transport,
     Reaction&& reaction)
    : transitions_transport{
      std::forward<Transitions_Transport>(transitions_transport) }
    , reaction{ std::forward<Reaction>(reaction) }
    , timestep{ transitions_transport.time_step() }
    {}
    
    Transitions_PTRW_Transport_Reaction
    (Transitions_Transport&& transitions_transport,
     Reaction reaction, double timestep)
    : transitions_transport{
      std::forward<Transitions_Transport>(transitions_transport) }
    , reaction{ std::forward<Reaction>(reaction) }
    , timestep{ timestep }
    {}

    template <typename State>
    void operator()(State& state)
    {
      transitions_transport(state);
      reaction(state, timestep);
    }
    
    void time_step(double timestep)
    {
      transitions_transport.time_step(timestep);
      this->timestep = timestep;
    }
    
    double time_step()
    { return timestep; }
    
    auto const& get_transitions_transport() const
    { return transitions_transport; }
    
    auto const& get_reaction() const
    { return reaction; }

  private:
    Transitions_Transport transitions_transport;
    Reaction reaction;
    double timestep;
  };
  template <typename Transitions_Transport, typename Reaction>
  Transitions_PTRW_Transport_Reaction
  (Transitions_Transport&&, Reaction&&) ->
  Transitions_PTRW_Transport_Reaction<Transitions_Transport, Reaction>;
  template <typename Transitions_Transport, typename Reaction>
  Transitions_PTRW_Transport_Reaction
  (Transitions_Transport&&, Reaction&&, double) ->
  Transitions_PTRW_Transport_Reaction<Transitions_Transport, Reaction>;
  
  /** \class Transitions_CTR_Transport_Reaction CTRW/Transitions.h "CTRW/Transitions.h"
   * \brief Update state according to a transport transition given state
   * followed by a reaction transition given change in state.time.
   * 
   * Note: transport transition should enforce any boundary conditions. */
  template <typename Transitions_Transport, typename Reaction>
  class Transitions_CTRW_Transport_Reaction
  {
  public:
    Transitions_CTRW_Transport_Reaction
    (Transitions_Transport&& transitions_transport,
     Reaction&& reaction)
    : transitions_transport{
      std::forward<Transitions_Transport>(transitions_transport) }
    , reaction{ std::forward<Reaction>(reaction) }
    {}

    template <typename State>
    void operator()(State& state)
    {
      auto time_old = state.time;
      transitions_transport(state);
      reaction(state, state.time-time_old);
    }
    
    auto const& get_transitions_transport() const
    { return transitions_transport; }
    
    auto const& get_reaction() const
    { return reaction; }

  private:
    Transitions_Transport transitions_transport;
    Reaction reaction;
  };
  template <typename Transitions_Transport, typename Reaction>
  Transitions_CTRW_Transport_Reaction
  (Transitions_Transport&&, Reaction&&) ->
  Transitions_CTRW_Transport_Reaction<Transitions_Transport, Reaction>;

  /** \class Transitions_Reaction_Position CTRW/Transitions.h "CTRW/Transitions.h"
   *  \brief Update state according to a reaction transition given a time step
   * returned by a TimeGenerator given state,
   * followed by a change in state.time according to the time step
   * and a change in position returned by a JumpGenerator given state.
   * 
   * Note: JumpGenerator should enforce any boundary conditions. */
	template <typename TimeGenerator,
  typename JumpGenerator, typename Reaction>
	class Transitions_Reaction_Position
	{
	public:
		Transitions_Reaction_Position
		(TimeGenerator&& time_generator,
     JumpGenerator&& jump_generator,
     Reaction&& reaction)
    : time_generator{
      std::forward<TimeGenerator>(time_generator) }
    , jump_generator{
      std::forward<JumpGenerator>(jump_generator) }
    , reaction{ std::forward<Reaction>(reaction) }
		{}

		template <typename State>
		void operator()(State& state)
		{
			double exposure_time = time_generator(state);
			reaction(state, exposure_time);
			state.position += jump_generator(state);
			state.time += exposure_time;
		}
    
    auto const& get_time_generator() const
    { return time_generator; }
    
    auto const& get_jump_generator() const
    { return jump_generator; }
    
    auto const& get_reaction() const
    { return reaction; }

	private:
		TimeGenerator time_generator;
		JumpGenerator jump_generator;
		Reaction reaction;
	};
  template <typename TimeGenerator,
  typename JumpGenerator, typename Reaction>
  Transitions_Reaction_Position
  (TimeGenerator&&, JumpGenerator&&, Reaction&&) ->
  Transitions_Reaction_Position<TimeGenerator, JumpGenerator, Reaction>;
  
  /** \class Transitions_Reaction_Position_Conditional CTRW/Transitions.h "CTRW/Transitions.h"
   * \brief Update state according to a reaction transition given a time step
   * returned by a TimeGenerator given state.
   * 
   * If the particle does not react, this is
   * followed by a change in state.time according to the time step
   * and a change in position returned by a JumpGenerator given state
   * Note: JumpGenerator should enforce any boundary conditions */
  template <typename TimeGenerator,
  typename JumpGenerator, typename Reaction>
  class Transitions_Reaction_Position_Conditional
  {
  public:
    Transitions_Reaction_Position_Conditional
    (TimeGenerator&& time_generator,
     JumpGenerator&& jump_generator,
     Reaction&& reaction)
    : time_generator{
      std::forward<TimeGenerator>(time_generator) }
    , jump_generator{
      std::forward<JumpGenerator>(jump_generator) }
    , reaction{ std::forward<Reaction>(reaction) }
    {}

    template <typename State>
    void operator()(State& state)
    {
      double exposure_time = time_generator(state);
      bool reacted = reaction(state, exposure_time);
      if (!reacted)
      {
        state.position += jump_generator(state);
        state.time += exposure_time;
      }
    }
    
    auto const& get_time_generator() const
    { return time_generator; }
    
    auto const& get_jump_generator() const
    { return jump_generator; }
    
    auto const& get_reaction() const
    { return reaction; }

  private:
    TimeGenerator time_generator;
    JumpGenerator jump_generator;
    Reaction reaction;
  };
  template <typename TimeGenerator,
  typename JumpGenerator, typename Reaction>
  Transitions_Reaction_Position_Conditional
  (TimeGenerator&&, JumpGenerator&&, Reaction&&) ->
  Transitions_Reaction_Position_Conditional<
  TimeGenerator, JumpGenerator, Reaction>;
  
  /** \class Transitions_RunTumble CTRW/Transitions.h "CTRW/Transitions.h"
   * \brief Update state.time according to state.run.   
   * 
   *  If state.run, update state.time by summing contribution returned
   * by TimeGenerator_run object given state, and state.position by
   * summing contribution returned by JumoGenerator object given state. \n
   * If not state.run, update state.time by summing contribution returned
   * by TimeGenerator_tumble object given state, and state.orientation by
   * summing contribution returned by OrientationGenerator object given state
   * Flip state.run */
  template
  <typename TimeGenerator_run, typename TimeGenerator_tumble,
  typename JumpGenerator, typename OrientationGenerator>
  class Transitions_RunTumble
  {
  public:
    Transitions_RunTumble
    (TimeGenerator_run&& time_generator_run,
     TimeGenerator_tumble&& time_generator_tumble,
     JumpGenerator&& jump_generator,
     OrientationGenerator&& orientation_generator)
    : time_generator_run{
      std::forward<TimeGenerator_run>(time_generator_run) }
    , time_generator_tumble{
      std::forward<TimeGenerator_tumble>(time_generator_tumble) }
    , jump_generator{
      std::forward<JumpGenerator>(jump_generator) }
    , orientation_generator{
      std::forward<OrientationGenerator>(orientation_generator) }
    {}

    template <typename State>
    void operator()(State& state)
    {
      if (state.run)
      {
        state.position += jump_generator(state);
        state.time += time_generator_run(state);
      }
      else
      {
        state.orientation += orientation_generator(state);
        state.time += time_generator_tumble(state);
      }
      state.run = !state.run;
    }
    
    auto const& get_time_generator_run() const
    { return time_generator_run; }
    
    auto const& get_time_generator_tumble() const
    { return time_generator_tumble; }
    
    auto const& get_jump_generator() const
    { return jump_generator; }
    
    auto const& get_orientation_generator() const
    { return orientation_generator; }

  private:
    TimeGenerator_run time_generator_run;
    TimeGenerator_tumble time_generator_tumble;
    JumpGenerator jump_generator;
    OrientationGenerator orientation_generator;
  };
  template
  <typename TimeGenerator_run, typename TimeGenerator_tumble,
  typename JumpGenerator, typename OrientationGenerator>
  Transitions_RunTumble
  (TimeGenerator_run&&, TimeGenerator_tumble&&,
   JumpGenerator&&, OrientationGenerator&&) ->
  Transitions_RunTumble<TimeGenerator_run, TimeGenerator_tumble,
  JumpGenerator, OrientationGenerator>;
  
  /** \class Transitions_RunAndTuble_PTRW CTRW/Transitions.h "CTRW/Transitions.h"
   *  \brief Update state according to current state.state.
   * 
   *  Update state according to current state.state:
   *  - 0 (run): Sum to state.position as returned
   *          by JumpGenerator object given state. \n
   *          Apply Boundary condition. \n
   *          If state.state has not been switched to wall-tumble,
   *         apply the StateSwitcher run rule given state
   *         to obtain new state.state.
   *  - 1 (tumble): Apply the StateSwitcher tumble rule given state
   *             to obtain new state.state. \n
   *             If state.state has been switched to run,
   *            sum to state.orientation as returned
   *             by OrientationGenerator object given state.
   *  - 2 (wall-tumble): Apply the StateSwitcher wall_tumble rule given state
   *                  to obtain new state.state. \n
   *                 If state.state has been switched to run,
   *                 sum to state.orientation as returned
   *                 by OrientationGenerator_Wall object given state. */
  template
  <typename Boundary, typename StateSwitcher,
  typename JumpGenerator, typename OrientationGenerator,
  typename OrientationGenerator_Wall = OrientationGenerator>
  class Transitions_RunAndTumble_PTRW
  {
  public:
    Transitions_RunAndTumble_PTRW
    (JumpGenerator&& jump_generator,
     OrientationGenerator&& orientation_generator,
     OrientationGenerator_Wall&& orientation_generator_wall,
     StateSwitcher&& state_switcher, Boundary&& boundary)
    : boundary{ std::forward<Boundary>(boundary) }
    , state_switcher{
      std::forward<StateSwitcher>(state_switcher) }
    , jump_generator{
      std::forward<JumpGenerator>(jump_generator) }
    , orientation_generator{
      std::forward<OrientationGenerator>(orientation_generator) }
    , orientation_generator_wall{
      std::forward<OrientationGenerator_Wall>(orientation_generator_wall) }
    {}
    
    Transitions_RunAndTumble_PTRW
    (JumpGenerator&& jump_generator,
     OrientationGenerator&& orientation_generator,
     StateSwitcher&& state_switcher, Boundary&& boundary)
    : Transitions_RunAndTumble_PTRW{
      std::forward<JumpGenerator>(jump_generator),
      std::forward<OrientationGenerator>(orientation_generator),
      std::forward<OrientationGenerator>(orientation_generator),
      std::forward<StateSwitcher>(state_switcher),
      std::forward<Boundary>(boundary) }
    {}
    
    template <typename State>
    void operator()(State& state)
    {
      auto const& state_old = state;
      switch (state.state)
      {
        case 0:
          operation::plus_InPlace(state.position,
                                  jump_generator(state));
          if(!boundary(state, state_old))
            state.state = state_switcher.run(state);
          break;
        case 1:
          state.state = state_switcher.tumble(state);
          if (state.state == 0)
            operation::plus_InPlace(state.orientation,
                                    orientation_generator(state));
          break;
        case 2:
          state.state = state_switcher.wall_tumble(state);
          if (state.state == 0)
            operation::plus_InPlace(state.orientation,
                                    orientation_generator_wall(state));
          break;
        default:
          break;
      }
    }
    
    auto const& get_boundary() const
    { return boundary; }
    
    auto const& get_state_switcher() const
    { return state_switcher; }
    
    auto const& get_orientation_generator() const
    { return orientation_generator; }
    
    auto const& get_orientation_generator_wall() const
    { return orientation_generator_wall; }
    
  private:
    Boundary boundary;
    StateSwitcher state_switcher;
    JumpGenerator jump_generator;
    OrientationGenerator orientation_generator;
    OrientationGenerator_Wall orientation_generator_wall;
  };
  template
  <typename Boundary, typename StateSwitcher,
  typename JumpGenerator, typename OrientationGenerator,
  typename OrientationGenerator_Wall>
  Transitions_RunAndTumble_PTRW
  (JumpGenerator&&, OrientationGenerator&&, OrientationGenerator_Wall&&,
  StateSwitcher&&, Boundary&&) ->
  Transitions_RunAndTumble_PTRW
  <Boundary, StateSwitcher, JumpGenerator,
  OrientationGenerator, OrientationGenerator_Wall>;
  template
  <typename Boundary, typename StateSwitcher,
  typename JumpGenerator, typename OrientationGenerator>
  Transitions_RunAndTumble_PTRW
  (JumpGenerator&&, OrientationGenerator&&, StateSwitcher&&, Boundary&&) ->
  Transitions_RunAndTumble_PTRW
  <Boundary, StateSwitcher, JumpGenerator,
  OrientationGenerator, OrientationGenerator>;
  
  template <typename Acceleration, typename Boundary>
  class Transitions_Velocity_Acceleration
  {
  public:
    Transitions_Velocity_Acceleration
    (double timestep,
     Acceleration&& acceleration, Boundary&& boundary = {})
    : timestep{ timestep }
    , acceleration{ std::forward<Acceleration>(acceleration) }
    , boundary{ std::forward<Boundary>(boundary) }
    {}

    template <typename State>
    void operator() (State& state)
    {
      auto state_old = state;
      operation::linearOp(time_step(), state.velocity, state.position,
                          state.position);
      operation::linearOp(time_step(), acceleration(state), state.velocity,
                          state.velocity);
      boundary(state, state_old);
    }
    
    template <typename Time = double>
    void time_step(Time timestep)
    { this->timestep = timestep; }
    
    auto time_step() const
    { return timestep; }
    
    auto const& get_acceleration() const
    { return acceleration; }
    
    auto const& get_boundary() const
    { return boundary; }

  private:
    double timestep;
    Acceleration acceleration;
    Boundary boundary;
  };
  template
  <typename Acceleration, typename Boundary>
  Transitions_Velocity_Acceleration
  (double, Acceleration&&, Boundary&&) ->
  Transitions_Velocity_Acceleration<Acceleration, Boundary>;
}

#endif /* CTRW_TRANSITIONS_H */
