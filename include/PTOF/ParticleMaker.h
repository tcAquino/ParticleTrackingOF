/**
   \file PTOF/ParticleMaker.h
   \author Tomas Aquino
   \date 06/03/2024
   \brief Objects to make particles.
*/

#ifndef PTOF_PARTICLEMAKER_H
#define PTOF_PARTICLEMAKER_H

#include "General/Meta.h"
#include "PTOF/State.h"
#include <label.H>
#include <utility>

namespace ptof {
/**
   \class ParticleMaker_Generic PTOF/ParticleMaker.h "PTOF/ParticleMaker.h"
   \brief Functor to make particles given position.
*/
template <typename Particle_t, typename Locator> struct ParticleMaker_Generic {
  using Particle = Particle_t;
  using State = typename Particle::State;
  using Position = typename State::Position;
  using Time = typename State::Time;
  using Mass = typename State::Mass;
  using Tag = typename State::Tag;
  using Info = typename State::Info;

  ParticleMaker_Generic(meta::Selector_t<Particle>, Locator &&locator,
                        Time time = {}, Mass mass = {}, Tag tag = {},
                        Info info = {})
      : locator{std::forward<Locator>(locator)}, time{time}, mass{mass},
        tag{tag}, info{info} {}

  template <typename Position>
  Particle operator()(Position const &position, Foam::label hint = -1) {
    return {make_state(State::make_position(position), info, locator, hint,
                       time, mass, tag++)};
  }

  Particle operator()(Position const &position, Foam::label hint = -1) {
    return {make_state(position, info, locator, hint, time, mass, tag++)};
  }

  Locator locator; /**< Locator to locate particle states in mesh.*/
  Time time;       /**< Particle state time.                      */
  Mass mass;       /**< Particle state mass.                      */
  Tag tag;         /**< Particle state tag.                       */
  Info info;       /**< Particle state info.                       */
};
template <typename Particle, typename Locator>
ParticleMaker_Generic(meta::Selector_t<Particle>, Locator &&locator,
                      typename Particle::State::Time,
                      typename Particle::State::Mass,
                      typename Particle::State::Tag)
    -> ParticleMaker_Generic<Particle, Locator>;
template <typename Particle, typename Locator>
ParticleMaker_Generic(meta::Selector_t<Particle>, Locator &&locator,
                      typename Particle::State::Time,
                      typename Particle::State::Mass)
    -> ParticleMaker_Generic<Particle, Locator>;
template <typename Particle, typename Locator>
ParticleMaker_Generic(meta::Selector_t<Particle>, Locator &&locator,
                      typename Particle::State::Time)
    -> ParticleMaker_Generic<Particle, Locator>;
template <typename Particle, typename Locator>
ParticleMaker_Generic(meta::Selector_t<Particle>, Locator &&locator)
    -> ParticleMaker_Generic<Particle, Locator>;

/** \class ParticleMaker_Periodic PTOF/ParticleMaker.h "PTOF/ParticleMaker.h"
 *  \brief Functor to make particles given position, for states with periodicity
 * information. */
template <typename Particle_t, typename Locator, typename Boundary>
struct ParticleMaker_Periodic {
  using Particle = Particle_t;
  using State = typename Particle::State;
  using Position = typename State::Position;
  using Time = typename State::Time;
  using Mass = typename State::Mass;
  using Tag = typename State::Tag;
  using Info = typename State::Info;

  ParticleMaker_Periodic(meta::Selector_t<Particle>, Locator &&locator,
                         Boundary &&boundary, Time time = {}, Mass mass = {},
                         Tag tag = {}, Info info = {})
      : locator{std::forward<Locator>(locator)},
        boundary{std::forward<Boundary>(boundary)}, time{time}, mass{mass},
        tag{tag}, info{info} {}

  Particle operator()(Position const &position, Foam::label hint = -1) {
    typename Particle::State state{
        make_state(position, {}, info, locator, hint, time, mass, tag++)};
    boundary(state);
    return state;
  }

  template <typename Position>
  Particle operator()(Position const &position, Foam::label hint = -1) {
    return (*this)(State::make_position(position), {}, info, locator, hint,
                   time, mass, tag++);
  }

  Locator locator;   /**< Locator to locate particle states in mesh.  */
  Boundary boundary; /**< Boundary to enforce periodicity given state */
  Time time;         /**< Particle state time.                        */
  Mass mass;         /**< Particle state mass.                        */
  Tag tag;           /**< Particle state tag.                         */
  Info info;         /**< Particle state info.                         */
};
template <typename Particle, typename Locator, typename Boundary>
ParticleMaker_Periodic(meta::Selector_t<Particle>, Locator &&locator,
                       Boundary &&boundary, typename Particle::State::Time,
                       typename Particle::State::Mass,
                       typename Particle::State::Tag,
                       typename Particle::State::Info)
    -> ParticleMaker_Periodic<Particle, Locator, Boundary>;
template <typename Particle, typename Locator, typename Boundary>
ParticleMaker_Periodic(meta::Selector_t<Particle>, Locator &&locator,
                       Boundary &&boundary, typename Particle::State::Time,
                       typename Particle::State::Mass,
                       typename Particle::State::Tag)
    -> ParticleMaker_Periodic<Particle, Locator, Boundary>;
template <typename Particle, typename Locator, typename Boundary>
ParticleMaker_Periodic(meta::Selector_t<Particle>, Locator &&locator,
                       Boundary &&boundary, typename Particle::State::Time,
                       typename Particle::State::Mass)
    -> ParticleMaker_Periodic<Particle, Locator, Boundary>;
template <typename Particle, typename Locator, typename Boundary>
ParticleMaker_Periodic(meta::Selector_t<Particle>, Locator &&locator,
                       Boundary &&boundary, typename Particle::State::Time)
    -> ParticleMaker_Periodic<Particle, Locator, Boundary>;
template <typename Particle, typename Locator, typename Boundary>
ParticleMaker_Periodic(meta::Selector_t<Particle>, Locator &&locator,
                       Boundary &&boundary)
    -> ParticleMaker_Periodic<Particle, Locator, Boundary>;
} // namespace ptof

#endif /* PTOF_PARTICLEMAKER_H */
