/**
 \file CTRW/Particle.h
 \author Tomás Aquino
 \date 26/09/2017
*/

#ifndef CTRW_PARTICLE_H
#define CTRW_PARTICLE_H

namespace ctrw {
/** \class Particle CTRW/Particle.h "CTRW/Particle.h"
  \brief Particle class for ptrw/ctrw/tdrw methods.
  \details
  A particle has an old state and a new state. \n
  Transitions are handled by an external class that modifies the state. \n
  Upon a transition, the old state becomes the previous new state
  and the new state is obtained by applying the transition to the previous new
 state. \tparam State_t Object holding particle state information.
*/
template <typename State_t> class Particle {
public:
  using State = State_t; /**< Particle state. */

  /** Constructor.
   \param state Particle state, copied to set both old and new state.
   */
  Particle(State const &state) : _state_new{state}, _state_old{state} {}

  /** Constructor
  \param state_new New particle state.
  \param state_old Old particle state.
  */
  Particle(State state_new, State state_old)
      : _state_new{state_new}, _state_old{state_old} {}

  /** \param transition Update particle states according to this object that
   * modifies a state.*/
  template <typename Transition> void transition(Transition &&transition) {
    _state_old = _state_new;
    transition(_state_new);
  }

  /** \param state Set new and old state to a copy of this state. */
  void set(State const &state) {
    _state_old = state;
    _state_new = state;
  }

  /** \param transformation Transformation to apply to both new and old state
   * individually.*/
  template <typename Transformation>
  void transform_both(Transformation &&transformation) {
    transformation(_state_new);
    transformation(_state_old);
  }

  /** \param transformation Transformation to apply to both states
   * simultaneously.*/
  template <typename Transformation>
  void transform(Transformation &&transformation) {
    transformation(_state_new, _state_old);
  }

  /** \return New state. */
  State const &state_new() const { return _state_new; }

  /** \return Old state. */
  State const &state_old() const { return _state_old; }

private:
  State_t _state_new; /**< New state. */
  State_t _state_old; /**< Old state. */
};
} // namespace ctrw

#endif /* CTRW_PARTICLE_H */
