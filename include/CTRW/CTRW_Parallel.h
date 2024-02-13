/**
* \file CTRW/CTRW_Parallel.h
* \author Tomás Aquino
* \date 14/12/2023
*/

#ifndef CTRW_PARALLEL_h
#define CTRW_PARALLEL_h

#include <functional>
#include <list>
#include <omp.h>
#include <vector>
#include "Particle.h"
#include "General/Useful.h"

namespace ctrw
{
  /** \class CTRW_Parallel CTRW/CTRW_Parallel.h "CTRW/CTRW_Parallel.h"
  * \brief Handles a set of particles characterized by a state, for parallel computations.
  *
  * Particles undergo transitions handled by a separate class
  * passed to the evolution methods. \n
  * State_t characterizes each particle's state. \n
  *
  * Notes: \n
  * - Order of particles in container is not preserved by removal methods.
  * - Iterators and references can be invalidated when removing or adding particles. */
  template <typename State_t>
  class CTRW_Parallel
  {
  public:
    using State = State_t;                   /**< Particle state.                           */
    using Particle = ctrw::Particle<State>;  /**< Particle type.                            */
    using Container = std::vector<Particle>; /**< Set of particles.                         */
    
    /** \struct Tag CTRW/CTRW.h "CTRW/CTRW.h"
    * To select constructors that tag particles.*/
    struct Tag{};

    /** Construct empty. */
    CTRW_Parallel()
    {}

    /** Construct given particles. */
    CTRW_Parallel(Container particles)
    : particle_container{ particles }
    {}
    
    /** Construct given particles and tag them in ascending order.
    * State must define: tag. */
    CTRW_Parallel(Container particles, Tag)
    : CTRW_Parallel{ particles }
    { retag(); }
    
    /** Construct given particles and reserve space for given maximum number. */
    CTRW_Parallel(Container particles, std::size_t max_nr_particles)
    {
      particle_container.reserve(max_nr_particles);
      particle_container = particles;
    }
    
    /** Construct given particles and reserve space for given maximum number
    * Tag them in ascending order.
    * State must define:
    * - tag */
    CTRW_Parallel(Container particles, std::size_t max_nr_particles, Tag)
    : CTRW_Parallel(particles, max_nr_particles)
    { retag(); }

    /** Handle nr_particles particles.
    * Each particle's state is produced by StateMaker. */
    template <typename StateMaker>
    CTRW_Parallel(std::size_t nr_particles, StateMaker&& state_maker)
    { push_back(nr_particles, state_maker); }

    /** Add nr_particles particles.
    * Each particle's state is produced by StateMaker. */
    template <typename StateMaker>
    void push_back(std::size_t nr_particles, StateMaker&& state_maker)
    {
      particle_container.reserve(nr_particles);
      std::generate_n(std::back_inserter(particle_container), nr_particles, state_maker);
    }

    /** Add copy of particle. */
    void push_back(Particle const& particle)
    { particle_container.push_back(particle); }

    /** Add copies of particles in a container. */
    void push_back(Container const& particles)
    {
      for (auto const& particle : particles)
        particle_container.push_back(particle);
    }
    
    /** Reserve space for nr_particles particles. */
    void reserve(std::size_t nr_particles)
    { particle_container.reserve(nr_particles); }
    
    /** Tag particles in order of container. */
    void retag()
    {
      std::size_t pp = 0;
      for (auto& particle : particle_container)
      {
        particle.transform_both([pp](State& state){ state.tag = pp; });
        ++pp;
      }
    }
    
    /** Set particle's new and old state to given state. */
    void set(std::size_t part, State const& state)
    { particle_container[part].set(state); }
    
    /** Apply transformation to particle. */
    template <typename Transformation>
    void transform(Transformation&& transformation, std::size_t part)
    { particle_container[part].transform(transformation); }
    
    /** Apply transformation to all particles. */
    template <typename Transformation>
    void transform(Transformation&& transformation)
    {
      for (auto& particle : particle_container)
        particle.transform(transformation);
    }
    
    /** Apply same transformation to both states of particle. */
    template <typename Transformation>
    void transform_both(Transformation&& transformation, std::size_t part)
    { particle_container[part].transform_both(transformation); }
    
    /** Apply same transformation to both states of all particles. */
    template <typename Transformation>
    void transform_both(Transformation&& transformation)
    {
      for (auto& particle : particle_container)
        particle.transform_both(transformation);
    }

    /** Remove particles satisfying a criterion. */
    template<typename Criterion>
    void remove(Criterion&& criterion)
    {
      for (std::size_t part = size(); part --> 0;)
        remove(part, criterion);
    }

    /** Remove particles at given positions in container. */
    template <typename IntegerType>
    void remove(std::list<IntegerType>& list)
    {
      list.sort(std::greater<IntegerType>{});
      for (auto const& part : list)
        remove(part);
    }

    /** Remove all particles. */
    void clear()
    { particle_container.clear(); }

    /** Particles make transitions until their time is >= time_to. */
    template<typename Transitions>
    void evolve_time
    (double time_to, Transitions&& transitions)
    {
      evolve
      ([time_to](Particle const& part)
       { return part.state_new().time < time_to; }, transitions);
    }

    /** Particles make transitions until their time is >= time_to. */
    template<typename Transitions>
    void evolve_space
    (double length_to, Transitions&& transitions)
    {
      evolve
      ([length_to](Particle const& part)
      { return part.state_new().position < length_to; }, transitions);
    }

    /** Particles make transitions until a given (particle-based) criterion is met. */
    template
    <template<typename...> typename Container, typename Transitions, typename Criterion,
    typename... Ts>
    void evolve
    (Criterion&& criterion, Container<Transitions, Ts...>& transitions)
    {
      #pragma omp parallel for
      for (std::size_t part = 0; part < particle_container.size(); ++part)
        while(criterion(particle_container[part]))
          particle_container[part].transition(transitions[omp_get_thread_num()]);
    }

    /** Each particle makes one transition. */
    template<template<typename...> typename Container, typename Transitions, typename... Ts>
    void step(Container<Transitions, Ts...>& transitions)
    {
      #pragma omp parallel for
      for (std::size_t part = 0; part < particle_container.size(); ++part)
        particle_container[part].transition(transitions[omp_get_thread_num()]);
    }
    
    /** Each particle for which criterion is satisfied makes one transition. */
    template
    <template<typename...> typename Container, typename Transitions, typename Criterion,
    typename... Ts>
    std::size_t step(Criterion&& criterion, Container<Transitions, Ts...>& transitions)
    {
      
      std::size_t nr_stepped = 0;
      #pragma omp parallel for reduction (+:nr_stepped)
      for (std::size_t part = 0; part < particle_container.size(); ++part)
        if (criterion(particle_container[part]))
        {
          particle_container[part].transition(transitions[omp_get_thread_num()]);
          ++nr_stepped;
        }
      return nr_stepped;
    }

    /** Get reference to particle container. */
    Container const& particles() const
    { return particle_container; }

    /** Get reference to particle. */
    Particle const& particles(std::size_t part) const
    { return particle_container[part]; }

    /** Get number of particles. */
    std::size_t size() const
    { return particle_container.size(); }

    /** Iterator to begining of particle container. */
    auto cbegin() const
    { return particle_container.cbegin(); }

    /** Iterator to end of particle container. */
    auto cend() const
    { return particle_container.cend(); }
    
    /** Iterator to begining of particle container. */
    auto begin() const
    { return particle_container.cbegin(); }

    /** Iterator to end of particle container. */
    auto end() const
    { return particle_container.cend(); }
    
  private:
    Container particle_container; /**< Current particles. */

    /** Remove particle by position in container if criterion is met. */
    template <typename Criterion>
    void remove(std::size_t part, Criterion&& criterion)
    { if (criterion(particle_container[part])) remove(part); }
    
    /** Remove particle by position in container. */
    void remove(std::size_t part)
    { useful::swap_erase(particle_container, part); }
  };
}

#endif /* CTRW_PARALLEL_h */
