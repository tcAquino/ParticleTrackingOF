/**
 \file PTOF/ParticleMaker.h
 \author Tomás Aquino
 \date 06/03/2024
*/

#ifndef PTOF_PARTICLEMAKER_H
#define PTOF_PARTICLEMAKER_H

#include <utility>
#include "General/Useful.h"

namespace ptof
{
  /** \class ParticleMaker PTOF/ParticleMaker.h "PTOF/ParticleMaker.h"
   * \brief Functor to make particles given position. */
  template <typename Particle, typename Locator>
  struct ParticleMaker
  {
    using State = typename Particle::State;
    using Time = typename State::Time;
    using Mass = typename State::Mass;
    using Tag = typename State::Tag;
    using Info = typename State::Info;
    
    ParticleMaker
    (meta::Selector_t<Particle>,
     Locator&& locator,
     Time time = {},
     Mass mass = {},
     Tag tag = {},
     Info info = {})
    : locator{ std::forward<Locator>(locator) }
    , time{ time }
    , mass{ mass }
    , tag{ tag }
    {}
    
    template <typename Position>
    Particle operator()
    (Position const& position)
    {
      return { make_state(State::make_position(position),
                          info, locator, time, mass, tag) };
    }
    
    Locator locator; /**< Locator to locate particle states in mesh.*/
    Time time;       /**< Particle state time.                      */
    Mass mass;       /**< Particle state mass.                      */
    Tag tag;         /**< Particle state tag.                       */
    Info info;       /**< Particle state info.                       */
  };
  template <typename Particle, typename Locator>
  ParticleMaker
  (meta::Selector_t<Particle>,
   Locator&& locator,
   typename Particle::State::Time,
   typename Particle::State::Mass,
   typename Particle::State::Tag) ->
  ParticleMaker<Particle, Locator>;
  template <typename Particle, typename Locator>
  ParticleMaker
  (meta::Selector_t<Particle>,
   Locator&& locator,
   typename Particle::State::Time,
   typename Particle::State::Mass) ->
  ParticleMaker<Particle, Locator>;
  template <typename Particle, typename Locator>
  ParticleMaker
  (meta::Selector_t<Particle>,
   Locator&& locator,
   typename Particle::State::Time) ->
  ParticleMaker<Particle, Locator>;
  template <typename Particle, typename Locator>
  ParticleMaker
  (meta::Selector_t<Particle>,
   Locator&& locator) ->
  ParticleMaker<Particle, Locator>;
  
  /** \class ParticleMaker_Periodic PTOF/ParticleMaker.h "PTOF/ParticleMaker.h"
   *  \brief Functor to make particles given position, for states with periodicity information. */
  template
  <typename Particle, typename Locator, typename Boundary>
  struct ParticleMaker_Periodic
  {
    using State = typename Particle::State;
    using Time = typename State::Time;
    using Mass = typename State::Mass;
    using Tag = typename State::Tag;
    using Info = typename State::Info;
    
    ParticleMaker_Periodic
    (meta::Selector_t<Particle>,
     Locator&& locator,
     Boundary&& boundary,
     Time time = {},
     Mass mass = {},
     Tag tag = {},
     Info info = {})
    : locator{ std::forward<Locator>(locator) }
    , boundary{ std::forward<Boundary>(boundary) }
    , time{ time }
    , mass{ mass }
    , tag{ tag }
    , info{ info }
    {}
    
    template <typename Position>
    Particle operator()
    (Position const& position)
    {
      auto state_position = State::make_position(position);
      typename Particle::State state{
        make_state(state_position,
                   std::vector<int>(state_position.size(), 0),
                   info, locator, time, mass, tag) };
      
      boundary(state);
      return state;
    }
    
    Locator locator;   /**< Locator to locate particle states in mesh.  */
    Boundary boundary; /**< Boundary to enforce periodicity given state */
    Time time;         /**< Particle state time.                        */
    Mass mass;         /**< Particle state mass.                        */
    Tag tag;           /**< Particle state tag.                         */
    Info info;         /**< Particle state info.                         */
  };
  template
  <typename Particle, typename Locator, typename Boundary>
  ParticleMaker_Periodic
  (meta::Selector_t<Particle>,
   Locator&& locator,
   Boundary&& boundary,
   typename Particle::State::Time,
   typename Particle::State::Mass,
   typename Particle::State::Tag) ->
  ParticleMaker_Periodic<Particle, Locator, Boundary>;
  template
  <typename Particle, typename Locator, typename Boundary>
  ParticleMaker_Periodic
  (meta::Selector_t<Particle>,
   Locator&& locator,
   Boundary&& boundary,
   typename Particle::State::Time,
   typename Particle::State::Mass) ->
  ParticleMaker_Periodic<Particle, Locator, Boundary>;
  template
  <typename Particle, typename Locator, typename Boundary>
  ParticleMaker_Periodic
  (meta::Selector_t<Particle>,
   Locator&& locator,
   Boundary&& boundary,
   typename Particle::State::Time) ->
  ParticleMaker_Periodic<Particle, Locator, Boundary>;
  template
  <typename Particle, typename Locator, typename Boundary>
  ParticleMaker_Periodic
  (meta::Selector_t<Particle>,
   Locator&& locator,
   Boundary&& boundary) ->
  ParticleMaker_Periodic<Particle, Locator, Boundary>;
}

#endif /* PTOF_PARTICLEMAKER_H */
