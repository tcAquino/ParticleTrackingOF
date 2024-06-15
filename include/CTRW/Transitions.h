/**
 \file CTRW/Transitions.h
 \author Tomás Aquino
 \date 08/06/2018

 \brief Objects to handle state transitions in particle tracking.

Transitions objects should implement the following minimal functionality:

\code{.cpp}
class Transitions
{
public:
   template <typename State>
   void operator() (State& state)
   {
     // Update state for a single transition
   }
};
\endcode
*/

#ifndef CTRW_TRANSITIONS_H
#define CTRW_TRANSITIONS_H

#include "CTRW/JumpGenerator.h"
#include "CTRW/TimeGenerator.h"
#include "General/Operations.h"
#include "General/Useful.h"
#include "Geometry/Boundary.h"
#include <limits>
#include <utility>

namespace ctrw {
/** \class Transition_Time_Position CTRW/Transitions.h "CTRW/Transitions.h"
 * \brief Update position and time according to specified jumps and time
 * increments. \tparam JumpGenerator Object to return jump given state. \tparam
 * TimeGenerator Object to return time increment given state. \tparam Boundary
 * Object to enforce boundary condition on new state given new state and old
 * state. \note State must define:
 * - \c position
 * - \c time */
template <typename TimeGenerator, typename JumpGenerator,
          typename Boundary = geom::Boundary_DoNothing>
class Transitions_Time_Position {
public:
  Transitions_Time_Position(TimeGenerator &&time_generator,
                            JumpGenerator &&jump_generator,
                            Boundary &&boundary = {})
      : _time_generator{std::forward<TimeGenerator>(time_generator)},
        _jump_generator{std::forward<JumpGenerator>(jump_generator)},
        _boundary{std::forward<Boundary>(boundary)} {}

  template <typename State> void operator()(State &state) {
    auto state_old = state;
    op::plus_inplace(state.position, _jump_generator(state));
    state.time += _time_generator(state);
    _boundary(state, state_old);
  }

  auto const &time_generator() const { return _time_generator; }

  auto const &jump_generator() const { return _jump_generator; }

  auto const &boundary() const { return _boundary; }

private:
  TimeGenerator _time_generator;
  JumpGenerator _jump_generator;
  Boundary _boundary;
};
template <typename TimeGenerator, typename JumpGenerator, typename Boundary>
Transitions_Time_Position(TimeGenerator &&, JumpGenerator &&, Boundary &&)
    -> Transitions_Time_Position<TimeGenerator, JumpGenerator, Boundary>;

/**
   \class Transitions_AdaptiveTimeStep_Time_Position CTRW/Transitions.h
"CTRW/Transitions.h"
\brief Update position and time according to specified jumps and time
increments, using state-based adaptive time steps.
\tparam TimeStepAdaptor Object to return and update time step given state, time
generator, and jump generator. \tparam TimeGenerator Object to return time
increment given state. \tparam JumpGenerator Object to return jump given state.
\tparam Boundary Object to enforce boundary condition on new state given new
state and old state. \note State must define:
- \c position
- \c time
*/
template <typename TimeStepAdaptor, typename TimeGenerator,
          typename JumpGenerator, typename Locator = useful::Empty,
          typename Boundary = geom::Boundary_DoNothing>
class Transitions_AdaptiveTimeStep_Time_Position {
public:
  Transitions_AdaptiveTimeStep_Time_Position(
      TimeStepAdaptor &&time_step_adaptor, TimeGenerator &&time_generator,
      JumpGenerator &&jump_generator, Locator &&locator = {},
      Boundary &&boundary = {})
      : _time_step_adaptor{std::forward<TimeStepAdaptor>(time_step_adaptor)},
        _time_generator{std::forward<TimeGenerator>(time_generator)},
        _jump_generator{std::forward<JumpGenerator>(jump_generator)},
        _locator{std::forward<Locator>(locator)}, _boundary{
                                                      std::forward<Boundary>(
                                                          boundary)} {}

  template <typename State> void operator()(State &state) {
    auto state_old = state;
    _time_step_adaptor(state, _time_generator, _jump_generator);
    op::plus_inplace(state.position, _jump_generator(state));
    locate(state);
    state.time += _time_generator(state);
    if (_boundary(state, state_old))
      locate(state);
  }

  auto const &time_generator() const { return _time_generator; }

  auto const &jump_generator() const { return _jump_generator; }

  auto const &boundary() const { return _boundary; }

  auto const &locator() const { return _locator; }

  double time_step() const { return _time_step; }

private:
  TimeStepAdaptor _time_step_adaptor;
  TimeGenerator _time_generator;
  JumpGenerator _jump_generator;
  Locator _locator;
  Boundary _boundary;
  double _time_step;

  template <typename State> void locate(State &state) {
    if constexpr (meta::has_cell_v<State>)
      state.cell = _locator(state);
  }
};
template <typename TimeStepAdaptor, typename TimeGenerator,
          typename JumpGenerator, typename Locator, typename Boundary>
Transitions_AdaptiveTimeStep_Time_Position(TimeStepAdaptor &&, TimeGenerator &&,
                                           JumpGenerator &&, Locator &&,
                                           Boundary &&)
    -> Transitions_AdaptiveTimeStep_Time_Position<
        TimeStepAdaptor, TimeGenerator, JumpGenerator, Locator, Boundary>;

/** \class Transitions_Positions CTRW/Transitions.h "CTRW/Transitions.h"
 \brief Update position according to specified jump increments.
 \tparam JumpGenerator Object to return jump given state.
 \tparam Boundary Object to enforce boundary condition on new state given new
 state and old state. \note State must define:
 - \c position */
template <typename JumpGenerator, typename Boundary = geom::Boundary_DoNothing>
class Transitions_Position {
public:
  Transitions_Position(JumpGenerator &&jump_generator, Boundary &&boundary = {})
      : _jump_generator{std::forward<JumpGenerator>(jump_generator)},
        _boundary{std::forward<Boundary>(boundary)} {}

  template <typename State> void operator()(State &state) {
    auto state_old = state;
    op::plus_inplace(state.position, _jump_generator(state));
    _boundary(state, state_old);
  }

  template <typename Time = double> void time_step(Time time_step) {
    _jump_generator.time_step(time_step);
  }

  auto time_step() const { return _jump_generator.time_step(); }

  auto const &jump_generator() const { return _jump_generator; }

  auto const &boundary() const { return _boundary; }

private:
  JumpGenerator _jump_generator;
  Boundary _boundary;
};
template <typename JumpGenerator, typename Boundary>
Transitions_Position(JumpGenerator &&, Boundary &&)
    -> Transitions_Position<JumpGenerator, Boundary>;

/** \class VelocityFromGenerator
 *  \brief Helper class to return velocity value from JumpGenerator object given
 * state. \tparam JumpGenerator Object to return jump given state. \note \p
 * JumpGenerator must define:
 *  - <tt>velocity(State const&)</tt> */
template <typename JumpGenerator> struct VelocityFromGenerator {
  VelocityFromGenerator(JumpGenerator const &jump_generator)
      : _jump_generator{jump_generator} {}

  template <typename State> auto operator()(State const &state) const {
    return _jump_generator.velocity(state);
  }

private:
  JumpGenerator const &_jump_generator;
};

/** \class Transitions_Position_VelocityStep CTRW/Transitions.h
 "CTRW/Transitions.h" \brief Update position according to velocity field with
 adapative time step to fix jump increments. \details Tiime step is determined
 according to  jump sizes divided by local velocity. Jump increments use forward
 Euler. \tparam JumpGenerator Object to return jump given state. \tparam
 VelocityField Object to return velocity value given state. \tparam Boundary
 Object to enforce boundary condition on new state given new state and old
 state. \note State must define:
 - \c position
 - \c time

 Velocity field may be passed to constructor or velocities may be obtained from
 jump_generator if JumpGenerator defines:
 - <tt>velocity(State const&)</tt> */
template <typename JumpGenerator, typename VelocityField,
          typename Boundary = geom::Boundary_DoNothing>
class Transitions_Position_VelocityStep {
public:
  Transitions_Position_VelocityStep(double step_length,
                                    JumpGenerator &&jump_generator,
                                    VelocityField &&velocity_field,
                                    Boundary &&boundary)
      : _step_length{step_length}, _jump_generator{std::forward<JumpGenerator>(
                                       jump_generator)},
        _velocity_field{std::forward<VelocityField>(velocity_field)},
        _boundary{std::forward<Boundary>(boundary)} {}

  Transitions_Position_VelocityStep(double step_length,
                                    JumpGenerator &&jump_generator,
                                    Boundary &&boundary = {})
      : _step_length{step_length}, _jump_generator{std::forward<JumpGenerator>(
                                       jump_generator)},
        _velocity_field{VelocityFromGenerator{this->jump_generator}},
        _boundary{std::forward<Boundary>(boundary)} {}

  template <typename State> void operator()(State &state) {
    auto state_old = state;
    _jump_generator.time_step(_step_length / op::abs(_velocity(state)));
    if (_jump_generator.time_step() ==
        std::numeric_limits<double>::infinity()) {
      state.time = std::numeric_limits<double>::infinity();
      return;
    }
    op::plus_inplace(state.position, _jump_generator(state));
    state.time += _jump_generator.time_step();
    _boundary(state, state_old);
  }

  template <typename State> auto velocity(State const &state) const {
    return _velocity_field(state);
  }

  auto time_step() const { return _jump_generator.time_step(); }

  auto time_step(double time_step) const {
    return _jump_generator.time_step(time_step);
  }

  auto const &velocity_field() const { return _velocity_field; }

  auto const &jump_generator() const { return _jump_generator; }

  auto const &boundary() const { return _boundary; }

private:
  double _step_length;
  JumpGenerator _jump_generator;
  VelocityField _velocity_field; /**< Velocity field given state. */
  Boundary _boundary;
};
template <typename JumpGenerator, typename VelocityField, typename Boundary>
Transitions_Position_VelocityStep(double, JumpGenerator &&, VelocityField &&,
                                  Boundary &&)
    -> Transitions_Position_VelocityStep<JumpGenerator, VelocityField,
                                         Boundary>;
template <typename JumpGenerator, typename Boundary>
Transitions_Position_VelocityStep(double, JumpGenerator &&, Boundary &&)
    -> Transitions_Position_VelocityStep<
        JumpGenerator, VelocityFromGenerator<JumpGenerator>, Boundary>;

/** \class Transitions_PTRW_FlowField_Diff CTRW/Transitions.h
 "CTRW/Transitions.h" \brief Update position according to prescribed flow field
 and time step. \details Jump increments use forward Euler. \tparam FlowField
 Object to return velocity value given position. \tparam Boundary Object to
 enforce boundary condition on new state given new state and old state. \note
 State must define:
 - \c position
 - \c time */
template <typename FlowField, typename Boundary = geom::Boundary_DoNothing,
          typename RNG = std::mt19937>
class Transitions_PTRW_FlowField_Diff {
private:
  auto make_diff_generators(std::vector<double> const &diff, double time_step) {
    std::vector<JumpGenerator_Diffusion_1d<RNG>> diff_generators;
    diff_generators.reserve(diff.size());
    for (auto const &diff_val : diff)
      diff_generators.emplace_back(diff_val, time_step);

    return diff_generators;
  }

public:
  Transitions_PTRW_FlowField_Diff(std::vector<double> const &diff,
                                  double time_step, FlowField &&flow_field,
                                  Boundary &&boundary = {})
      : _time_generator{time_step}, _diff_generators{make_diff_generators(
                                        diff, time_step)},
        _flow_field{std::forward<FlowField>(flow_field)},
        _boundary{std::forward<Boundary>(boundary)} {}

  template <typename State> void operator()(State &state) {
    auto state_old = state;
    std::vector<double> jump(state.position.size());
    for (std::size_t dd = 0; dd < jump.size(); ++dd)
      jump[dd] = _diff_generators[dd]();
    jump[0] += _velocity(state) * _time_generator.time_step();

    for (std::size_t dd = 0; dd < jump.size(); ++dd)
      state.position[dd] += jump[dd];

    _boundary(state, state_old);
  }

  template <typename State> auto velocity(State const &state) const {
    return _velocity_field(state.position);
  }

  void time_step(double time_step) {
    _time_generator.time_step(time_step);
    for (auto &gen : _diff_generators)
      gen.time_step(time_step);
  }

  void diff(std::size_t dd, double diff) {
    return _diff_generators[dd].diff(diff);
  }

  double time_step() const { return _time_generator.time_step(); }

  double diff(std::size_t dd) const { return _diff_generators[dd].diff(); }

  auto const &time_generator() const { return _time_generator; }

  auto const &boundary() const { return _boundary; }

  auto const &velocity_field() const { return _flow_field; }

private:
  TimeGenerator_Step<double> _time_generator;
  std::vector<JumpGenerator_Diffusion_1d<RNG>> _diff_generators;
  FlowField _flow_field;
  Boundary _boundary;
};
template <typename FlowField, typename Boundary>
Transitions_PTRW_FlowField_Diff(std::vector<double> const &, double,
                                FlowField &&, Boundary &&)
    -> Transitions_PTRW_FlowField_Diff<FlowField, Boundary, std::mt19937>;

/** \class Transitions_PTRW_Diffusion_1d CTRW/Transitions.h "CTRW/Transitions.h"
 \brief Update position according to diffusion in one dimension with prescribed
 time step. \details Jump increments use stochastic forward Euler. \tparam
 Boundary Object to enforce boundary condition on new state given new state and
 old state. \tparam RNG Random number generator. \note State must define:
 - \c position
 - \c time */
template <typename Boundary = geom::Boundary_DoNothing,
          typename RNG = std::mt19937>
class Transitions_PTRW_Diffusion_1d {
public:
  Transitions_PTRW_Diffusion_1d(double diff, double time_step,
                                Boundary &&boundary = {})
      : _time_generator{time_step}, _jump_generator{diff, time_step},
        _boundary{std::forward<Boundary>(boundary)} {}

  template <typename State> void operator()(State &state) {
    auto state_old = state;
    state.position += _jump_generator(state);
    _boundary(state, state_old);
    state.time += time_step();
  }

  void time_step(double time_step) {
    _time_generator.time_step(time_step);
    _jump_generator.time_step(time_step);
  }

  void diff(double diff) { return _jump_generator.diff(diff); }

  double time_step() const { return _time_generator.time_step(); }

  double diff() const { return _jump_generator.diff(); }

  auto const &time_generator() const { return _time_generator; }

  auto const &jump_generator() const { return _jump_generator; }

  auto const &boundary() const { return _boundary; }

private:
  TimeGenerator_Step<double> _time_generator;
  JumpGenerator_Diffusion_1d<RNG> _jump_generator;
  Boundary _boundary;
};
template <typename Boundary>
Transitions_PTRW_Diffusion_1d(double, double, Boundary &&)
    -> Transitions_PTRW_Diffusion_1d<Boundary, std::mt19937>;

/** \class Transitions_PTRW_Transport_Reaction CTRW/Transitions.h
 "CTRW/Transitions.h" \brief Update state according to specified transport
 transitions followed by reaction. \tparam Transitions_Transport Update particle
 state according to transport. \tparam Reaction Update particle state according
 to reaction, given state and time increment. \note
 - \p Transitions_Transport should enforce any boundary conditions.
 - Time step needed for reaction may be specified or extracted from transport
 transitions if \p Transitions_Transport defines:
  -# \c time_step() */
template <typename Transitions_Transport, typename Reaction>
class Transitions_PTRW_Transport_Reaction {
public:
  Transitions_PTRW_Transport_Reaction(
      Transitions_Transport &&transitions_transport, Reaction &&reaction)
      : _transitions_transport{std::forward<Transitions_Transport>(
            transitions_transport)},
        _reaction{std::forward<Reaction>(reaction)},
        _time_step{transitions_transport.time_step()} {}

  Transitions_PTRW_Transport_Reaction(
      Transitions_Transport &&transitions_transport, Reaction reaction,
      double time_step)
      : _transitions_transport{std::forward<Transitions_Transport>(
            transitions_transport)},
        _reaction{std::forward<Reaction>(reaction)}, _time_step{time_step} {}

  template <typename State> void operator()(State &state) {
    _transitions_transport(state);
    _reaction(state, _time_step);
  }

  void time_step(double time_step) {
    _transitions_transport.time_step(time_step);
    _time_step = time_step;
  }

  double time_step() { return _time_step; }

  auto const &get_transitions_transport() const {
    return _transitions_transport;
  }

  auto const &get_reaction() const { return _reaction; }

private:
  Transitions_Transport _transitions_transport;
  Reaction _reaction;
  double _time_step;
};
template <typename Transitions_Transport, typename Reaction>
Transitions_PTRW_Transport_Reaction(Transitions_Transport &&, Reaction &&)
    -> Transitions_PTRW_Transport_Reaction<Transitions_Transport, Reaction>;
template <typename Transitions_Transport, typename Reaction>
Transitions_PTRW_Transport_Reaction(Transitions_Transport &&, Reaction &&,
                                    double)
    -> Transitions_PTRW_Transport_Reaction<Transitions_Transport, Reaction>;

/** \class Transitions_CTRW_Transport_Reaction CTRW/Transitions.h
 "CTRW/Transitions.h" \brief Update state according to specified transport
 transitions followed by reaction. \tparam Transitions_Transport Update particle
 state according to transport. \tparam Reaction Update particle state according
 to reaction, given state and time increment. \note State must define:
 - \c time

 \p Transitions_Transport should enforce any boundary conditions. */
template <typename Transitions_Transport, typename Reaction>
class Transitions_CTRW_Transport_Reaction {
public:
  Transitions_CTRW_Transport_Reaction(
      Transitions_Transport &&transitions_transport, Reaction &&reaction)
      : _transitions_transport{std::forward<Transitions_Transport>(
            transitions_transport)},
        _reaction{std::forward<Reaction>(reaction)} {}

  template <typename State> void operator()(State &state) {
    auto time_old = state.time;
    _transitions_transport(state);
    _reaction(state, state.time - time_old);
  }

  auto const &transitions_transport() const { return _transitions_transport; }

  auto const &reaction() const { return _reaction; }

private:
  Transitions_Transport _transitions_transport;
  Reaction _reaction;
};
template <typename Transitions_Transport, typename Reaction>
Transitions_CTRW_Transport_Reaction(Transitions_Transport &&, Reaction &&)
    -> Transitions_CTRW_Transport_Reaction<Transitions_Transport, Reaction>;

/** \class Transitions_Reaction_Position CTRW/Transitions.h "CTRW/Transitions.h"
 * \brief Update state according to reaction, followed by a position update,
 * followed by a time update \tparam TimeGenerator Object to return time
 * increment given state. \tparam JumpGenerator Object to return jump given
 * state. \tparam Reaction Update particle state according to reaction, given
 * state and time increment. \note State must define:
 * - \c position
 * - \c time
 *
 * \p JumpGenerator should enforce any boundary conditions. */
template <typename TimeGenerator, typename JumpGenerator, typename Reaction>
class Transitions_Reaction_Position {
public:
  Transitions_Reaction_Position(TimeGenerator &&time_generator,
                                JumpGenerator &&jump_generator,
                                Reaction &&reaction)
      : _time_generator{std::forward<TimeGenerator>(time_generator)},
        _jump_generator{std::forward<JumpGenerator>(jump_generator)},
        _reaction{std::forward<Reaction>(reaction)} {}

  template <typename State> void operator()(State &state) {
    double exposure_time = _time_generator(state);
    _reaction(state, exposure_time);
    state.position += _jump_generator(state);
    state.time += exposure_time;
  }

  auto const &time_generator() const { return _time_generator; }

  auto const &jump_generator() const { return _jump_generator; }

  auto const &reaction() const { return _reaction; }

private:
  TimeGenerator _time_generator;
  JumpGenerator _jump_generator;
  Reaction _reaction;
};
template <typename TimeGenerator, typename JumpGenerator, typename Reaction>
Transitions_Reaction_Position(TimeGenerator &&, JumpGenerator &&, Reaction &&)
    -> Transitions_Reaction_Position<TimeGenerator, JumpGenerator, Reaction>;

/** \class Transitions_Reaction_Position_Conditional CTRW/Transitions.h
 * "CTRW/Transitions.h" \brief Update state according to reaction, followed by a
 * position update, followed by a time update \tparam TimeGenerator Object to
 * return time increment given state. \tparam JumpGenerator Object to return
 * jump given state. \tparam Reaction Update particle state according to
 * reaction, given state and time increment. Must return true when applied if
 * reaction occurs and false otherwise. \note State must define:
 * - \c position
 * - \c time
 *
 * \p JumpGenerator should enforce any boundary conditions. */
template <typename TimeGenerator, typename JumpGenerator, typename Reaction>
class Transitions_Reaction_Position_Conditional {
public:
  Transitions_Reaction_Position_Conditional(TimeGenerator &&time_generator,
                                            JumpGenerator &&jump_generator,
                                            Reaction &&reaction)
      : _time_generator{std::forward<TimeGenerator>(time_generator)},
        _jump_generator{std::forward<JumpGenerator>(jump_generator)},
        _reaction{std::forward<Reaction>(reaction)} {}

  template <typename State> void operator()(State &state) {
    double exposure_time = _time_generator(state);
    bool reacted = _reaction(state, exposure_time);
    if (!reacted) {
      state.position += _jump_generator(state);
      state.time += exposure_time;
    }
  }

  auto const &time_generator() const { return _time_generator; }

  auto const &jump_generator() const { return _jump_generator; }

  auto const &reaction() const { return _reaction; }

private:
  TimeGenerator _time_generator;
  JumpGenerator _jump_generator;
  Reaction _reaction;
};
template <typename TimeGenerator, typename JumpGenerator, typename Reaction>
Transitions_Reaction_Position_Conditional(TimeGenerator &&, JumpGenerator &&,
                                          Reaction &&)
    -> Transitions_Reaction_Position_Conditional<TimeGenerator, JumpGenerator,
                                                 Reaction>;

/** \class Transitions_RunTumble CTRW/Transitions.h "CTRW/Transitions.h"
 * \brief Update position and time when in run state and orientation and time
 * when not in run state. Flip run state at each update. returned by a
 * TimeGenerator given state. Position and time are then update only if no
 * reaction occured. \tparam TimeGenerator_run Object to return time increment
 * when in run mode, given state. \tparam TimeGenerator_tumble Object to return
 * time increment when not run mode, given state. \tparam JumpGenerator Object
 * to return jump given state. \tparam OrientationGenerator Object to return
 * orientation increments given state. \note State must define:
 * - \c position
 * - \c time
 * - \c orientation
 * - \c run
 *
 * \p JumpGenerator should enforce any boundary conditions. */
template <typename TimeGenerator_run, typename TimeGenerator_tumble,
          typename JumpGenerator, typename OrientationGenerator>
class Transitions_RunAndTumble {
public:
  Transitions_RunAndTumble(TimeGenerator_run &&time_generator_run,
                           TimeGenerator_tumble &&time_generator_tumble,
                           JumpGenerator &&jump_generator,
                           OrientationGenerator &&orientation_generator)
      : _time_generator_run{std::forward<TimeGenerator_run>(
            time_generator_run)},
        _time_generator_tumble{
            std::forward<TimeGenerator_tumble>(time_generator_tumble)},
        _jump_generator{std::forward<JumpGenerator>(jump_generator)},
        _orientation_generator{
            std::forward<OrientationGenerator>(orientation_generator)} {}

  template <typename State> void operator()(State &state) {
    if (state.run) {
      state.position += _jump_generator(state);
      state.time += _time_generator_run(state);
    } else {
      state.orientation += _orientation_generator(state);
      state.time += _time_generator_tumble(state);
    }
    state.run = !state.run;
  }

  auto const &time_generator_run() const { return _time_generator_run; }

  auto const &time_generator_tumble() const { return _time_generator_tumble; }

  auto const &jump_generator() const { return _jump_generator; }

  auto const &orientation_generator() const { return _orientation_generator; }

private:
  TimeGenerator_run _time_generator_run;
  TimeGenerator_tumble _time_generator_tumble;
  JumpGenerator _jump_generator;
  OrientationGenerator _orientation_generator;
};
template <typename TimeGenerator_run, typename TimeGenerator_tumble,
          typename JumpGenerator, typename OrientationGenerator>
Transitions_RunAndTumble(TimeGenerator_run &&, TimeGenerator_tumble &&,
                         JumpGenerator &&, OrientationGenerator &&)
    -> Transitions_RunAndTumble<TimeGenerator_run, TimeGenerator_tumble,
                                JumpGenerator, OrientationGenerator>;

/** \class Transitions_RunAndTumble_PTRW CTRW/Transitions.h "CTRW/Transitions.h"
 *  \brief Update state according to current internal state.
 * \details
 *  When internal state is:
 *  - 0 (run): Update position as according to \p JumpGenerator object. \n
 *          Apply \p Boundary condition. \n
 *          If internal state has not been switched to wall-tumble,
 *          apply the \p StateSwitcher run rule given state
 *          to obtain new internal state.
 *  - 1 (tumble): Apply the \p StateSwitcher tumble rule given state
 *           to obtain new internal state. \n
 *           If internal state has been switched to run,
 *           update orientation accoiding to \p OrientationGenerator .
 *  - 2 (wall-tumble): Apply \p StateSwitcher wall_tumble rule
 *               to obtain new internal state. \n
 *               If internal state has been switched to run,
 *               update orientation accoiding to \p OrientationGenerator_Wall .
 * \tparam Boundary Object to enforce boundary condition on new state given new
 * state and old state. \tparam StateSwitcher Object to return new internal
 * state given state. \tparam JumpGenerator Object to return jump given state.
 * \tparam OrientationGenerator Object to return orientation increments given
 * state. \tparam OrientationGenerator_Wall Object to return orientation
 * increments given state. \note State must define:
 *  - \c position
 *  - \c orientation
 *  - \c state
 * \c StateSwitcher must define:
 *  - <tt>int run(State const&)</tt>
 *  - <tt>int tumble(State const&)</tt>
 *  - <tt>int wall_tumble(State const&)</tt>
 *
 * \c Boundary object should set wall_tumble state when appropriate. */
template <typename Boundary, typename StateSwitcher, typename JumpGenerator,
          typename OrientationGenerator, typename OrientationGenerator_Wall>
class Transitions_RunAndTumble_PTRW {
public:
  Transitions_RunAndTumble_PTRW(
      JumpGenerator &&jump_generator,
      OrientationGenerator &&orientation_generator,
      OrientationGenerator_Wall &&orientation_generator_wall,
      StateSwitcher &&state_switcher, Boundary &&boundary)
      : _boundary{std::forward<Boundary>(boundary)},
        _state_switcher{std::forward<StateSwitcher>(state_switcher)},
        _jump_generator{std::forward<JumpGenerator>(jump_generator)},
        _orientation_generator{
            std::forward<OrientationGenerator>(orientation_generator)},
        _orientation_generator_wall{std::forward<OrientationGenerator_Wall>(
            orientation_generator_wall)} {}

  Transitions_RunAndTumble_PTRW(JumpGenerator &&jump_generator,
                                OrientationGenerator &&orientation_generator,
                                StateSwitcher &&state_switcher,
                                Boundary &&boundary)
      : Transitions_RunAndTumble_PTRW{
            std::forward<JumpGenerator>(jump_generator),
            std::forward<OrientationGenerator>(orientation_generator),
            std::forward<OrientationGenerator>(orientation_generator),
            std::forward<StateSwitcher>(state_switcher),
            std::forward<Boundary>(boundary)} {}

  template <typename State> void operator()(State &state) {
    auto const &state_old = state;
    switch (state.state) {
    case 0:
      op::plus_inplace(state.position, _jump_generator(state));
      if (!boundary(state, state_old))
        state.state = _state_switcher.run(state);
      break;
    case 1:
      state.state = _state_switcher.tumble(state);
      if (state.state == 0)
        op::plus_inplace(state.orientation, _orientation_generator(state));
      break;
    case 2:
      state.state = _state_switcher.wall_tumble(state);
      if (state.state == 0)
        op::plus_inplace(state.orientation, _orientation_generator_wall(state));
      break;
    default:
      break;
    }
  }

  auto const &boundary() const { return _boundary; }

  auto const &state_switcher() const { return _state_switcher; }

  auto const &orientation_generator() const { return _orientation_generator; }

  auto const &orientation_generator_wall() const {
    return _orientation_generator_wall;
  }

private:
  Boundary _boundary;
  StateSwitcher _state_switcher;
  JumpGenerator _jump_generator;
  OrientationGenerator _orientation_generator;
  OrientationGenerator_Wall _orientation_generator_wall;
};
template <typename Boundary, typename StateSwitcher, typename JumpGenerator,
          typename OrientationGenerator, typename OrientationGenerator_Wall>
Transitions_RunAndTumble_PTRW(JumpGenerator &&, OrientationGenerator &&,
                              OrientationGenerator_Wall &&, StateSwitcher &&,
                              Boundary &&)
    -> Transitions_RunAndTumble_PTRW<Boundary, StateSwitcher, JumpGenerator,
                                     OrientationGenerator,
                                     OrientationGenerator_Wall>;
template <typename Boundary, typename StateSwitcher, typename JumpGenerator,
          typename OrientationGenerator>
Transitions_RunAndTumble_PTRW(JumpGenerator &&, OrientationGenerator &&,
                              StateSwitcher &&, Boundary &&)
    -> Transitions_RunAndTumble_PTRW<Boundary, StateSwitcher, JumpGenerator,
                                     OrientationGenerator,
                                     OrientationGenerator>;

/** \class Transitions_Velocity_Acceleration CTRW/Transitions.h
 * "CTRW/Transitions.h" \brief Update position and velocity according to imposed
 * acceleration. \tparam Acceleration Object to return acceleration given state.
 * \tparam Boundary Object to enforce boundary condition on new state given new
 * state and old state. \note State must define:
 * - \c position
 * - \c velocity */
template <typename Acceleration, typename Boundary = useful::Empty>
class Transitions_Velocity_Acceleration {
public:
  Transitions_Velocity_Acceleration(double time_step,
                                    Acceleration &&acceleration,
                                    Boundary &&boundary)
      : _time_step{time_step}, _acceleration{std::forward<Acceleration>(
                                   acceleration)},
        _boundary{std::forward<Boundary>(boundary)} {}

  Transitions_Velocity_Acceleration(double time_step,
                                    Acceleration &&acceleration)
      : _time_step{time_step}, _acceleration{std::forward<Acceleration>(
                                   acceleration)},
        _boundary{useful::Empty{}} {}

  template <typename State> void operator()(State &state) {
    auto state_old = state;
    op::linearop(time_step(), state.velocity, state.position, state.position);
    op::linearop(time_step(), acceleration(state), state.velocity,
                 state.velocity);
    boundary(state, state_old);
  }

  template <typename Time = double> void time_step(Time time_step) {
    _time_step = time_step;
  }

  auto time_step() const { return _time_step; }

  auto const &acceleration() const { return _acceleration; }

  auto const &boundary() const { return _boundary; }

private:
  double _time_step;
  Acceleration _acceleration;
  Boundary _boundary;
};
template <typename Acceleration, typename Boundary>
Transitions_Velocity_Acceleration(double, Acceleration &&, Boundary &&)
    -> Transitions_Velocity_Acceleration<Acceleration, Boundary>;
template <typename Acceleration>
Transitions_Velocity_Acceleration(double, Acceleration &&)
    -> Transitions_Velocity_Acceleration<Acceleration, useful::Empty>;
} // namespace ctrw

#endif /* CTRW_TRANSITIONS_H */
