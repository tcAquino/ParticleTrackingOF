/**
 * @file   Particle.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Tue Sep 26 00:00:00 2017
 *
 * @brief Particle for CTRW dynamics.
 */

#ifndef CTRW_PARTICLE_H
#define CTRW_PARTICLE_H

namespace ctrw {
/**
 * @brief Particle class for CTRW.
 *
 * @details
 * - A particle has an old state and a new state.
 *
 * - Transitions are handled by an external class that modifies the state.
 *
 * - Upon a transition, the old state becomes the previous new state and the new
 *   state is obtained by applying the transition to the previous new state.
 *
 * @tparam State_t Object holding particle state information.
 */
template <typename State_t> class Particle {
public:
  using State = State_t; /**< Particle state. */

  /**
   * @brief Constructor.
   *
   * @param state Particle state, copied to set both old and new state.
   */
  Particle(State const &state) : _state_new{state}, _state_old{state} {}

  /**
   * @brief Constructor
   *
   * @param state_new New particle state.
   *
   * @param state_old Old particle state.
   */
  Particle(State state_new, State state_old)
      : _state_new{state_new}, _state_old{state_old} {}

  /**
   * @brief Update particle states.
   *
   * @details Sets the old state to the new state and updates the new state.
   *
   * @param transition Object that modifies a state.
   */
  template <typename Transition> void transition(Transition &&transition) {
    _state_old = _state_new;
    transition(_state_new);
  }

  /** @brief Set new and old state to \p state. */
  void set(State const &state) {
    _state_old = state;
    _state_new = state;
  }

  /**
   * @param transformation Transformation to apply to both new and old state
   *                       individually.
   */
  template <typename Transformation>
  void transform_both(Transformation &&transformation) {
    transformation(_state_new);
    transformation(_state_old);
  }

  /**
   * @param transformation Transformation that takes and may modify both new and
   *                       old state.
   */
  template <typename Transformation>
  void transform(Transformation &&transformation) {
    transformation(_state_new, _state_old);
  }

  /** @return New state. */
  State const &state_new() const { return _state_new; }

  /** @return Old state. */
  State const &state_old() const { return _state_old; }

private:
  State_t _state_new; /**< New state. */
  State_t _state_old; /**< Old state. */
};
} // namespace ctrw

#endif /* CTRW_PARTICLE_H */
