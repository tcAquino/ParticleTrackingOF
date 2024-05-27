/**
 \file CTRW/CTRW_Parallel.h
 \author Tomás Aquino
 \date 14/12/2023
*/

#ifndef CTRW_PARALLEL_h
#define CTRW_PARALLEL_h

#include "General/Useful.h"
#include "Particle.h"
#include <functional>
#include <list>
#include <omp.h>
#include <vector>

namespace ctrw {
/** \class CTRW_Parallel CTRW/CTRW_Parallel.h "CTRW/CTRW_Parallel.h"
 \brief Handles a set of particles characterized by a state, for parallel
computations. \details Particles undergo transitions handled by a separate class
passed to the evolution methods. \tparam State_t Characterizes a particle's
state. \note
 - Order of particles in container is not preserved by removal methods.
 - Iterators and references can be invalidated when removing or adding
particles.
 - Particle evolution methods expect a container of Transitions objects, one per
parallel thread.
 */
template <typename State_t> class CTRW_Parallel {
public:
  using State = State_t; /**< Particle state.                           */
  using Particle = ctrw::Particle<State>;  /**< Particle type.  */
  using Container = std::vector<Particle>; /**< Set of particles. */

  /** \struct Tag CTRW/CTRW.h "CTRW/CTRW.h"
  \brief To select constructors that tag particles. */
  struct Tag {};

  /** Constructor.
   \brief No particles.
  */
  CTRW_Parallel() {}

  /** Constructor.
   \param particles Container with particles.
  */
  CTRW_Parallel(Container particles) : _particles{particles} {}

  /** Constructor.
   \brief Tag particles in ascending order.
   \note State must define:
   - \c tag
   \param particles Container with particles.
  */
  CTRW_Parallel(Container particles, Tag) : CTRW_Parallel{particles} {
    retag();
  }

  /** Constructor.
   \brief Reserve space for maximum number of particles.
   \param particles Container with particles.
   \param max_nr_particles Number of particles to reserve space for.
  */
  CTRW_Parallel(Container particles, std::size_t max_nr_particles) {
    _particles.reserve(max_nr_particles);
    _particles = particles;
  }

  /** Constructor.
   \brief Tag particles in ascending order reserve space for maximum number of
   particles. \note State must define:
   - \c tag
   \param particles Container with particles.
   \param max_nr_particles Number of particles to reserve space for.
  */
  CTRW_Parallel(Container particles, std::size_t max_nr_particles, Tag)
      : CTRW_Parallel(particles, max_nr_particles) {
    retag();
  }

  /** Constructor.
   \param nr_particles Number of particles
   \param state_maker Object to initialize particle states.
  */
  template <typename StateMaker>
  CTRW_Parallel(std::size_t nr_particles, StateMaker &&state_maker) {
    push_back(nr_particles, state_maker);
  }

  /** \brief Add particles.
   \param nr_particles Number of particles to add.
   \param state_maker Object to initialize particle states.
  */
  template <typename StateMaker>
  void push_back(std::size_t nr_particles, StateMaker &&state_maker) {
    _particles.reserve(nr_particles);
    std::generate_n(std::back_inserter(_particles), nr_particles, state_maker);
  }

  /** \brief Add a particle.
   \param particle Particle to add a copy of.
  */
  void push_back(Particle const &particle) { _particles.push_back(particle); }

  /** \brief Add particles.
   \param particles Container of particles to add copies of.
  */
  void push_back(Container const &particles) {
    for (auto const &particle : particles)
      _particles.push_back(particle);
  }

  /** \brief Reserve space for particles.
   \param nr_particles Number of particles to reserve space for.
   */
  void reserve(std::size_t nr_particles) { _particles.reserve(nr_particles); }

  /** \brief Tag particles in order of container. */
  void retag() {
    std::size_t pp = 0;
    for (auto &particle : _particles) {
      particle.transform_both([pp](State &state) { state.tag = pp; });
      ++pp;
    }
  }

  /** \brief Set particle's new and old state.
   \param part Index of particle to change state of.
   \param state New state.
  */
  void set(std::size_t part, State const &state) {
    _particles[part].set(state);
  }

  /** \brief Apply transformation to a particle's state(s).
   \param transformation Transformation to apply.
   \param part Index of particle to apply transformation to.
   */
  template <typename Transformation>
  void transform(Transformation &&transformation, std::size_t part) {
    _particles[part].transform(transformation);
  }

  /** \brief Apply transformation to all particles' state(s).
  \param transformation Transformation to apply.
  */
  template <typename Transformation>
  void transform(Transformation &&transformation) {
    for (auto &particle : _particles)
      particle.transform(transformation);
  }

  /** \brief Apply same transformation to both states of a particle.
  \param transformation Transformation to apply.
  \param part Index of particle to apply transformation to.
  */
  template <typename Transformation>
  void transform_both(Transformation &&transformation, std::size_t part) {
    _particles[part].transform_both(transformation);
  }

  /** \brief Apply same transformation to both states of all particles.
  \param transformation Transformation to apply.
  */
  template <typename Transformation>
  void transform_both(Transformation &&transformation) {
    for (auto &particle : _particles)
      particle.transform_both(transformation);
  }

  /** \brief Remove particles satisfying a criterion.
   \param criterion Criterion to apply to each particle.
  */
  template <typename Criterion> void remove(Criterion &&criterion) {
    for (std::size_t part = size(); part-- > 0;)
      remove(part, criterion);
  }

  /** \brief Remove particles at given indices.
   \param list List of particle indices to remove.
  */
  template <typename IntegerType> void remove(std::list<IntegerType> &list) {
    list.sort(std::greater<IntegerType>{});
    for (auto const &part : list)
      remove(part);
  }

  /** \brief Remove all particles. */
  void clear() { _particles.clear(); }

  /** \brief Evolve particle states in time.
   \note State must define:
   - \c time
   \param time_to Particles make transitions until their time is >= time_to.
   \param transitions Handle changes to particle states in each step.
  */
  template <typename TransitionsContainer>
  void evolve_time(double time_to, TransitionsContainer &&transitions) {
    evolve(
        [time_to](Particle const &part) {
          return part.state_new().time < time_to;
        },
        transitions);
  }

  /** \brief Evolve particle states in space in 1D.
   \note
   State must define:
   - \c position [scalar]
   \param length_to Particles make transitions until their position is >=
   length_to. \param transitions Handle changes to particle states in each step.
  */
  template <typename TransitionsContainer>
  void evolve_space(double length_to, TransitionsContainer &&transitions) {
    evolve(
        [length_to](Particle const &part) {
          return part.state_new().position < length_to;
        },
        transitions);
  }

  /** \brief Evolve each particle until a given criterion is met.
   \param criterion Criterion to apply to each particle.
   \param transitions Handle changes to particle states in each step.
  */
  template <typename Criterion, typename TransitionsContainer>
  void evolve(Criterion &&criterion, TransitionsContainer &&transitions) {
#pragma omp parallel for
    for (std::size_t part = 0; part < _particles.size(); ++part)
      while (criterion(_particles[part]))
        _particles[part].transition(transitions[omp_get_thread_num()]);
  }

  /** \brief Evolve each particle one step.
   \param transitions Handle changes to particle states in each step.
  */
  template <typename TransitionsContainer>
  void step(TransitionsContainer &&transitions) {
#pragma omp parallel for
    for (std::size_t part = 0; part < _particles.size(); ++part)
      _particles[part].transition(transitions[omp_get_thread_num()]);
  }

  /** \brief Evolve each particle that satisfies a criterion one step.
   \param criterion Criterion to apply to each particle.
   \param transitions Handle changes to particle states in each step.
  */
  template <typename Criterion, typename TransitionsContainer>
  std::size_t step(Criterion &&criterion, TransitionsContainer &&transitions) {

    std::size_t nr_stepped = 0;
#pragma omp parallel for reduction(+ : nr_stepped)
    for (std::size_t part = 0; part < _particles.size(); ++part)
      if (criterion(_particles[part])) {
        _particles[part].transition(transitions[omp_get_thread_num()]);
        ++nr_stepped;
      }
    return nr_stepped;
  }

  /** \return Particle container.  */
  Container const &particles() const { return _particles; }

  /**
   \param part Index of particle.
   \return Particle at given index.
  */
  Particle const &particles(std::size_t part) const { return _particles[part]; }

  /** \return Number of particles. */
  std::size_t size() const { return _particles.size(); }

  /** \return Constant iterator to begining of particle container. */
  auto cbegin() const { return _particles.cbegin(); }

  /** \return Constant iterator to end of particle container. */
  auto cend() const { return _particles.cend(); }

  /** \return Constant iterator to begining of particle container. */
  auto begin() const { return _particles.cbegin(); }

  /** \return Constant iterator to end of particle container. */
  auto end() const { return _particles.cend(); }

private:
  Container _particles; /**< Current particles. */

  /** \brief Remove particle if criterion is met.
   \param part Index of particle to remove.
   \param criterion Criterion to apply to each particle.
   */
  template <typename Criterion>
  void remove(std::size_t part, Criterion &&criterion) {
    if (criterion(_particles[part]))
      remove(part);
  }

  /** \brief Remove particle.
   \param part Index of particle to remove.
   */
  void remove(std::size_t part) { useful::swap_erase(_particles, part); }
};
} // namespace ctrw

#endif /* CTRW_PARALLEL_h */
