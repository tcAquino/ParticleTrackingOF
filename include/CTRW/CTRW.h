/**
 * @file   CTRW.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Tue Sep 26 00:00:00 2017
 *
 * @brief CTRW particle dynamics
 */

#ifndef CTRW_CTRW_H
#define CTRW_CTRW_H

#include "CTRW/Particle.h"
#include "General/Parallel.h"
#include "General/Useful.h"
#include <functional>
#include <list>
#include <type_traits>
#include <vector>

/** @namespace ctrw Objects and methods related to CTRW. */
namespace ctrw {
/**
 * @brief Handles a set of particles characterized by a state.
 *
 * @details Particles undergo transitions handled by a separate class passed to
 *          the evolution methods.
 *
 * @tparam State_t Characterizes a particle's state.
 *
 * @tparam ParallelOption_t Choose serial or parallel implementations.
 *
 * @note
 *
 * - Order of particles in container is not preserved by removal methods.
 *
 * - Iterators and references can be invalidated when removing or adding
 *   particles.
 */
template <typename State_t,
          typename ParallelOption_t = par::ParallelOptions::Serial>
class CTRW {
public:
  using State = State_t; /**< Particle state. */
  using Particle = ctrw::Particle<State>;  /**< Particle type. */
  using Container = std::vector<Particle>; /**< Set of particles. */
  using ParallelOption =
      ParallelOption_t; /** To choose serial or parallel implementations. */

  /** @brief To select constructors that tag particles. */
  struct Tag {};

  /**
   * @brief Constructor.
   *
   * @details No particles.
   */
  CTRW() {}

  /**
   * @brief Constructor.
   *
   * @param particles Container with particles.
   */
  CTRW(Container particles) : _particles{particles} {}

  /**
   * @brief Constructor.
   *
   * @details Tag particles in ascending order.
   *
   * @param particles Container with particles.
   *
   * @note State must define:
   *
   * - \c tag
   */
  CTRW(Container particles, Tag) : CTRW{particles} { retag(); }

  /**
   * @brief Constructor.
   *
   * @details Reserve space for maximum number of particles.
   *
   * @param particles Container with particles.
   *
   * @param max_nr_particles Number of particles to reserve space for.
   */
  CTRW(Container particles, std::size_t max_nr_particles) {
    _particles.reserve(max_nr_particles);
    _particles = particles;
  }

  /**
   * @brief Constructor.
   *
   * @details Tag particles in ascending order reserve space for maximum number
   *          of particles.
   *
   * @param particles Container with particles.
   *
   * @param max_nr_particles Number of particles to reserve space for.
   *
   * @note State must define:
   *
   * - \c tag
   */
  CTRW(Container particles, std::size_t max_nr_particles, Tag)
      : CTRW(particles, max_nr_particles) {
    retag();
  }

  /**
   * @brief Constructor.
   * @param nr_particles Number of particles.
   *
   * @param state_maker Object to initialize particle states.
   */
  template <typename StateMaker>
  CTRW(std::size_t nr_particles, StateMaker &&state_maker) {
    push_back(nr_particles, state_maker);
  }

  /**
   * @brief Add particles.
   *
   * @param nr_particles Number of particles to add.
   *
   * @param state_maker Object to initialize particle states.
   */
  template <typename StateMaker>
  void push_back(std::size_t nr_particles, StateMaker &&state_maker) {
    _particles.reserve(nr_particles);
    std::generate_n(std::back_inserter(_particles), nr_particles, state_maker);
  }

  /**
   * @brief Add a particle.
   *
   * @param particle Particle to add a copy of.
   */
  void push_back(Particle const &particle) { _particles.push_back(particle); }

  /**
   * @brief Add particles.
   *
   * @param particles Container of particles to add copies of.
   */
  void push_back(Container const &particles) {
    for (auto const &particle : particles)
      _particles.push_back(particle);
  }

  /**
   * @brief Reserve space for particles.
   *
   * @param nr_particles Number of particles to reserve space for.
   */
  void reserve(std::size_t nr_particles) { _particles.reserve(nr_particles); }

  /** @brief Tag particles in order of container. */
  void retag() {
    std::size_t pp = 0;
    for (auto &particle : _particles) {
      particle.transform_both([pp](State &state) { state.tag = pp; });
      ++pp;
    }
  }

  /**
   * @brief Set particle's new and old state.
   *
   * @param part Index of particle to change state of.
   *
   * @param state New state.
   */
  void set(std::size_t part, State const &state) {
    _particles[part].set(state);
  }

  /**
   * @brief Apply transformation to a particle's states.
   *
   * @param transformation Transformation to apply.
   *
   * @param part Index of particle to apply transformation to.
   */
  template <typename Transformation>
  void transform(Transformation &&transformation, std::size_t part) {
    _particles[part].transform(transformation);
  }

  /**
   * @brief Apply transformation to all particles' states.
   *
   * @param transformation Transformation to apply.
   */
  template <typename Transformation>
  void transform(Transformation &&transformation) {
    for (auto &particle : _particles) {
      particle.transform(transformation);
    }
  }

  /**
   * @brief Apply same transformation to both states of a particle.
   *
   * @param transformation Transformation to apply.
   *
   * @param part Index of particle to apply transformation to.
   */
  template <typename Transformation>
  void transform_both(Transformation &&transformation, std::size_t part) {
    _particles[part].transform_both(transformation);
  }

  /**
   * @brief Apply same transformation to both states of all particles.
   *
   * @param transformation Transformation to apply.
   */
  template <typename Transformation>
  void transform_both(Transformation &&transformation) {
    for (auto &particle : _particles) {
      particle.transform_both(transformation);
    }
  }

  /**
   * @brief Remove particles satisfying a criterion.
   *
   * @param criterion Criterion to apply to each particle.
   */
  template <typename Criterion> void remove(Criterion &&criterion) {
    for (std::size_t part = size(); part-- > 0;) {
      remove(part, criterion);
    }
  }

  /**
   * @brief Remove particles at given indices.
   *
   * @param list List of particle indices to remove.
   */
  template <typename IntegerType> void remove(std::list<IntegerType> &list) {
    list.sort(std::greater<IntegerType>{});
    for (auto const &part : list) {
      remove(part);
    }
  }

  /** @brief Remove all particles. */
  void clear() { _particles.clear(); }

  /**
   * @brief Evolve particle states in time.
   *
   * @note State must define:
   *
   * - \c time
   *
   * @param time_to Particles make transitions until their time is <tt>>=
   *                time_to.</tt>
   *
   * @param transitions Handle changes to particle states in each step.
   */
  template <typename Transitions>
  void evolve_time(double time_to, Transitions &&transitions) {
    evolve(
        [time_to](Particle const &part) {
          return part.state_new().time < time_to;
        },
        transitions);
  }

  /**
   * @brief Evolve particle states in space in 1D.
   *
   * @param length_to Particles make transitions until their position is <tt>>=
   *                  length_to</tt>.
   *
   * @param transitions Handle changes to particle states in each
   *                    step.
   *
   * @note State must define:
   *
   * - \c position [scalar].
   */
  template <typename Transitions>
  void evolve_space(double length_to, Transitions &&transitions) {
    evolve(
        [length_to](Particle const &part) {
          return part.state_new().position < length_to;
        },
        transitions);
  }

  /**
   * @brief Evolve each particle until a given criterion is met.
   *
   * @param criterion Criterion to apply to each particle.
   *
   * @param transitions Handle changes to particle states in each step.
   */
  template <typename Transitions, typename Criterion>
  void evolve(Criterion &&criterion, Transitions &&transitions) {
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Parallel>) {
#ifdef _OPENMP
#pragma omp parallel for
#endif /** _OPENMP */
      for (std::size_t part = 0; part < size(); ++part) {
        while (criterion(_particles[part]))
          _particles[part].transition(transitions);
      }
    } else {
      for (auto &part : _particles) {
        while (criterion(part)) {
          part.transition(transitions);
        }
      }
    }
  }

  /**
   * @brief Evolve each particle one step.
   *
   * @param transitions Handle changes to particle states in each step.
   *
   * @return Number of particles stepped (all particles).
   */
  template <typename Transitions> std::size_t step(Transitions &&transitions) {
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Parallel>) {
#ifdef _OPENMP
#pragma omp parallel for
#endif /** _OPENMP */
      for (std::size_t part = 0; part < size(); ++part) {
        _particles[part].transition(transitions);
      }
    } else {
      for (auto &part : _particles) {
        part.transition(transitions);
      }
    }
    return size();
  }

  /**
   * @brief Evolve each particle that satisfies a criterion one step.
   *
   * @param criterion Criterion to apply to each particle.
   *
   * @param transitions Handle changes to particle states in each step.
   */
  template <typename Criterion, typename Transitions>
  std::size_t step(Criterion &&criterion, Transitions &&transitions) {
    std::size_t nr_stepped = 0;
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Parallel>) {
#ifdef _OPENMP
#pragma omp parallel for reduction(+ : nr_stepped)
#endif /** _OPENMP */
      for (std::size_t part = 0; part < size(); ++part) {
        if (criterion(_particles[part])) {
          _particles[part].transition(transitions);
          ++nr_stepped;
        }
      }
    } else {
      for (auto &part : _particles) {
        if (criterion(part)) {
          part.transition(transitions);
          ++nr_stepped;
        }
      }
    }
    return nr_stepped;
  }

  /**
   * @brief Evolve specified particles one step.
   *
   * @param tags Tags of particles to evolve.
   *
   * @param transitions Handle changes to particle states in each step.
   *
   * @return Number of particles stepped (size of \c tags).
   */
  template <typename Container, typename Transitions>
  std::size_t step_specified(Container const &tags, Transitions &&transitions) {
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Parallel>) {
#ifdef _OPENMP
#pragma omp parallel for
#endif /** _OPENMP */
      for (std::size_t part = 0; part < tags.size(); ++part) {
        _particles[tags[part]].transition(transitions);
      }
    } else {
      for (auto tag : tags) {
        _particles[tag].transition(transitions);
      }
    }
    return tags.size();
  }

  /** @return Particle container. */
  Container const &particles() const { return _particles; }

  /**
   * @param part Index of particle.
   *
   * @return Particle at given index.
   */
  Particle const &particles(std::size_t part) const { return _particles[part]; }

  /** @return Number of particles. */
  std::size_t size() const { return _particles.size(); }

  /** @return Constant iterator to begining of particle container. */
  auto cbegin() const { return _particles.cbegin(); }

  /** @return Constant iterator to end of particle container. */
  auto cend() const { return _particles.cend(); }

  /** @return Constant iterator to begining of particle container. */
  auto begin() const { return _particles.cbegin(); }

  /** @return Constant iterator to end of particle container. */
  auto end() const { return _particles.cend(); }

private:
  Container _particles; /**< Current particles. */

  /**
   * @brief Remove particle if criterion is met.
   *
   * @param part Index of particle to remove.
   *
   * @param criterion Criterion to apply to each particle.
   */
  template <typename Criterion>
  void remove(std::size_t part, Criterion &&criterion) {
    if (criterion(_particles[part])) {
      remove(part);
    }
  }

  /**
   * @brief Remove particle.
   *
   * @param part Index of particle to remove.
   */
  void remove(std::size_t part) { useful::swap_erase(_particles, part); }
};
} // namespace ctrw

#endif /* CTRW_CTRW_H */
