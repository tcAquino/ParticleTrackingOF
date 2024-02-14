/**
* \file CTRW/Particle.h
* \author Tomás Aquino
* \date 26/09/2017
*/

#ifndef CTRW_PARTICLE_H
#define CTRW_PARTICLE_H

namespace ctrw
{
  /** \class Particle CTRW/Particle.h "CTRW/Particle.h"
   * \brief Particle class for ptrw/ctrw/tdrw methods.
   * 
   * A particle has an old state and a new state. \n
   * Transitions are handled by an external class that modifies the state. \n
   * Upon a transition, the old state becomes the previous new state
   * and the new state is obtained by applying the transition to the previous new state. */
  template <typename State_t>
  class Particle
  {
  public:

    using State = State_t;

    Particle(State const& state)
    : new_state{ state }
    , old_state{ state }
    {}
    
    Particle(State state_new, State state_old)
    : new_state{ state_new }
    , old_state{ state_old }
    {}

    /** Transition according to a given rule that modifies the state.*/
    template <typename Transition>
    void transition(Transition&& transition)
    {
      old_state = new_state;
      transition(new_state);
    }
    
    /** Set both states to given state. */
    void set(State const& state)
    {
      old_state = state;
      new_state = state;
    }
    
    /** Apply same transformation to both states. */
    template <typename Transformation>
    void transform_both(Transformation&& transformation)
    {
      transformation(old_state);
      transformation(new_state);
    }
    
    /** Apply arbitrary transformation to both states. */
    template <typename Transformation>
    void transform(Transformation&& transformation)
    {
      transformation(new_state, old_state);
    }

    State const& state_new() const
    { return new_state; }

    State const& state_old() const
    { return old_state; }

  private:
    State_t new_state;
    State_t old_state;
  };
}


#endif /* CTRW_PARTICLE_H */
