/**
   \file PTOF/Steppers.h
   \author Tomas Aquino
   \date 09/03/2022
   \brief Objects and utilities to handle time stepping along advective and
   diffusive trajectories.
*/

#ifndef PTOF_STEPPERS_H
#define PTOF_STEPPERS_H

#include "CTRW/JumpGenerator.h"
#include "CTRW/Meta.h"
#include "General/Meta.h"
#include "General/Parallel.h"
#include "PTOF/TimeGenerator.h"
#include "PTOF/Useful.h"
#include <algorithm>
#include <cstddef>
#include <type_traits>
#include <utility>

namespace ptof {
/**
   \struct CTRWStepper PTOF/Steppers.h
   "PTOF/Steppers.h"
   \brief Types to choose options of CTRW stepping mode.
*/
struct CTRWSteppers {
  /**
     \struct Asynchronous PTOF/Steppers.h
     "PTOF/Steppers.h"
     \brief Step particles independtly
  */
  struct Asynchronous {
    template <typename CTRW, typename Transitions, typename SolverParameters,
              typename Time, typename Update = meta::DoNothing>
    static void evolve(CTRW &ctrw, Transitions &&transitions,
                       SolverParameters const &parameters_solvers,
                       Time new_time, Time old_time, Update &&update = {}) {
      ctrw.evolve(
          [new_time](typename CTRW::Particle const &part) {
            return part.state_new().time < new_time &&
                   !part.state_new().info.absorbed;
          },
          transitions);
      update(ctrw, new_time, old_time);
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output << io::line() << "CTRW Stepper\n"
             << io::line() << method << "\n"
             << io::line();
      return output;
    }

    inline static std::string method =
        "Asynchronous stepping"; /**< Name or description of stepping method. */
  };

  /**
     \struct ParticleStep PTOF/Steppers.h
     "PTOF/Steppers.h"
     \brief Evolve CTRW particles with synchronization over particle steps.
  */
  struct ParticleStep {
    template <typename CTRW, typename Transitions, typename SolverParameters,
              typename Time, typename Update = meta::DoNothing>
    static void evolve(CTRW &ctrw, Transitions &&transitions,
                       SolverParameters const &parameters_solvers,
                       Time new_time, Time old_time, Update &&update = {}) {
      std::size_t nr_stepped = 1;
      while (nr_stepped != 0) {
        nr_stepped = ctrw.step(
            [new_time](typename CTRW::Particle const &part) {
              return part.state_new().time < new_time &&
                     !part.state_new().info.absorbed;
            },
            transitions);
        update(ctrw, new_time, old_time);
      }
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output << io::line() << "CTRW Stepper\n"
             << io::line() << method << "\n"
             << io::line();
      return output;
    }

    inline static std::string method =
        "Particle step synchronization"; /**< Name or description of stepping
                                          * method.
                                          */
  };

  /**
     \struct ParticleTime PTOF/Steppers.h
     "PTOF/Steppers.h"
     \brief Evolve CTRW particles by stepping earliest-time particles first and
     synchronizing after at most one particle per thread has stepped.
  */
  struct ParticleTime {
    template <typename CTRW, typename Transitions, typename SolverParameters,
              typename Time, typename Update = meta::DoNothing>
    static void evolve(CTRW &ctrw, Transitions &&transitions,
                       SolverParameters const &parameters_solvers,
                       Time new_time, Time old_time, Update &&update = {}) {
      std::size_t nr_stepped = 1;
      while (nr_stepped != 0) {
        nr_stepped = ctrw.step_specified(
            find_earliest_time_tags(ctrw, new_time, old_time), transitions);
        update(ctrw, new_time, old_time);
      }
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output << io::line() << "CTRW Stepper\n"
             << io::line() << method << "\n"
             << io::line();
      return output;
    }

    inline static std::string method =
        "Particle step synchronization"; /**< Name or description of stepping
                                            method. */

  private:
    template <typename CTRW>
    static std::vector<std::size_t> find_earliest_time_tags(CTRW const &ctrw,
                                                            double new_time,
                                                            double old_time) {
      std::vector<
          std::pair<typename CTRW::State::Tag, typename CTRW::State::Time>>
          tags_times;
      tags_times.reserve(ctrw.size());
      for (auto const &part : ctrw) {
        if (part.state_new().time < new_time &&
            !part.state_new().info.absorbed) {
          tags_times.push_back({part.state_new().tag, part.state_old().time});
        }
      }

      std::size_t nr_tags =
          std::min(tags_times.size(),
                   par::get_num_threads(typename CTRW::ParallelOption{}));
      if (nr_tags == 0) {
        return {};
      }

      std::vector<std::size_t> tags;
      tags.reserve(nr_tags);
      std::nth_element(tags_times.begin(), tags_times.begin() + nr_tags - 1,
                       tags_times.end());

      for (auto const &tag_time : tags_times) {
        tags.push_back(tag_time.first);
      }

      return tags;
    }
  };

  /**
     \struct TimeStep PTOF/Steppers.h
     "PTOF/Steppers.h"
     \brief Evolve CTRW particles with synchronization over prescribed time
     step.
  */
  struct TimeStep {
    template <typename CTRW, typename Transitions, typename SolverParameters,
              typename Time, typename Update = meta::DoNothing>
    static void evolve(CTRW &ctrw, Transitions &&transitions,
                       SolverParameters const &parameters_solvers,
                       Time new_time, Time old_time, Update &&update = {}) {
      double current_time = old_time;
      while (current_time < new_time) {
        current_time += parameters_solvers.ctrw_time_step;
        ctrw.evolve(
            [current_time](typename CTRW::Particle const &part) {
              return part.state_new().time < current_time &&
                     !part.state_new().info.absorbed;
            },
            transitions);
        update(ctrw, new_time, old_time);
      }
    }
  };

  /**
       \brief Output generic information about object.
       \param output Output stream.
    */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "CTRW Stepper\n"
           << io::line() << method << "\n"
           << io::line();
    return output;
  }

  inline static std::string method =
      "Time step synchronization"; /**< Name or description of stepping method.
                                    */
};

struct Steppers {
  struct Euler {};
  struct RK2 {};
  struct RK4 {};
  struct Heun {};
};

/**
   \struct Steppers_Advection_Euler_Diffusion_Euler PTOF/Steppers.h
   "PTOF/Steppers.h"
   \brief Time steppers for forward Euler advection and stochastic forward
   Euler diffusion.
*/
struct Steppers_Advection_Euler_Diffusion_Euler {
  Steppers_Advection_Euler_Diffusion_Euler() = delete;

  /**
     \brief Adaptive time step generator.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
     \return Adaptive time step TimeGenerator.
  */
  template <typename Geometry, typename VelocityField, typename Boundary,
            typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  static auto makeTimeGenerator(Geometry const &geometry,
                                VelocityField const &velocity_field,
                                Boundary const &boundary,
                                TransportParameters const &params_transport,
                                ReactionParameters const &params_reaction,
                                SolverParameters const &params_solvers) {
    return TimeGenerator_Adaptive(geometry, velocity_field,
                                  boundary.surface_reaction, params_transport,
                                  params_reaction, params_solvers);
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters.
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Advection--diffusion JumpGenerator.
     \note
     - \p params_transport must define:
         - \c diff_coeff
     - If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
   */
  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto
  makeJumpGenerator(VelocityField &&velocity_field, Boundary &&boundary,
                    TransportParameters const &params_transport,
                    SolverParameters const &params_solvers, std::size_t dim) {
    return ctrw::JumpGenerator_Add{
        makeJumpGenerator_Advection<ParallelOption>(
            std::forward<VelocityField>(velocity_field),
            std::forward<Boundary>(boundary), params_transport, params_solvers,
            dim),
        makeJumpGenerator_Diffusion<ParallelOption>(
            std::forward<Boundary>(boundary), params_transport, params_solvers,
            dim)};
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure advection JumpGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto makeJumpGenerator_Advection(
      VelocityField &&velocity_field, Boundary &&boundary,
      TransportParameters const &params_transport,
      SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>) {
      return ctrw::JumpGenerator_Velocity{
          std::forward<VelocityField>(velocity_field), params_solvers.time_step,
          ParallelOption{}};
    }
    return ctrw::JumpGenerator_Velocity{
        std::forward<VelocityField>(velocity_field), 0., ParallelOption{}};
  }

  /**
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure diffusion JumpGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename ParallelOption, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto makeJumpGenerator_Diffusion(
      Boundary &&boundary, TransportParameters const &params_transport,
      SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>) {
      return ctrw::JumpGenerator_Diffusion{params_transport.diff_coeff,
                                           params_solvers.time_step, dim,
                                           ParallelOption{}};
    }
    return ctrw::JumpGenerator_Diffusion{params_transport.diff_coeff, 0., dim,
                                         ParallelOption{}};
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Steppers\n"
           << io::line() << "Advection: " << advection_method << "\n"
           << "Diffusion: " << diffusion_method << "\n"
           << io::line();
    return output;
  }

  inline static std::string advection_method =
      "Euler"; /**< Name or description of advection method. */
  inline static std::string diffusion_method =
      "Stochastic Euler"; /**< Name or discription of diffusion method. */
};

/**
   \struct Steppers_Advection_RK2_Diffusion_Euler PTOF/Steppers.h
   "PTOF/Steppers.h"
   \brief Time steppers for RK2 advection and stochastic forward Euler
   diffusion.
*/
struct Steppers_Advection_RK2_Diffusion_Euler {
  Steppers_Advection_RK2_Diffusion_Euler() = delete;

  /**
     \brief Adaptive time step generator.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
     \return Adaptive time step TimeGenerator.
  */
  template <typename Geometry, typename VelocityField, typename Boundary,
            typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  static auto makeTimeGenerator(Geometry const &geometry,
                                VelocityField const &velocity_field,
                                Boundary const &boundary,
                                TransportParameters const &params_transport,
                                ReactionParameters const &params_reaction,
                                SolverParameters const &params_solvers) {
    return TimeGenerator_Adaptive(geometry, velocity_field,
                                  boundary.surface_reaction, params_transport,
                                  params_reaction, params_solvers);
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters.
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Advection--diffusion JumpGenerator.
     \note
     - \p params_transport must define:
         - \c diff_coeff
     - If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
   */
  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto
  makeJumpGenerator(VelocityField &&velocity_field, Boundary &&boundary,
                    TransportParameters const &params_transport,
                    SolverParameters const &params_solvers, std::size_t dim) {
    return ctrw::JumpGenerator_Add{
        makeJumpGenerator_Advection<ParallelOption>(
            std::forward<VelocityField>(velocity_field),
            std::forward<Boundary>(boundary), params_transport, params_solvers,
            dim),
        makeJumpGenerator_Diffusion<ParallelOption>(
            std::forward<Boundary>(boundary), params_transport, params_solvers,
            dim)};
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure advection JumpGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto makeJumpGenerator_Advection(
      VelocityField &&velocity_field, Boundary &&boundary,
      TransportParameters const &params_transport,
      SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>) {
      return ctrw::JumpGenerator_Velocity_RK2{
          std::forward<VelocityField>(velocity_field), params_solvers.time_step,
          std::forward<Boundary>(boundary), ParallelOption{}};
    }
    return ctrw::JumpGenerator_Velocity_RK2{
        std::forward<VelocityField>(velocity_field), 0.,
        std::forward<Boundary>(boundary), ParallelOption{}};
  }

  /**
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure diffusion JumpGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename ParallelOption, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto makeJumpGenerator_Diffusion(
      Boundary &&boundary, TransportParameters const &params_transport,
      SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>) {
      return ctrw::JumpGenerator_Diffusion{params_transport.diff_coeff,
                                           params_solvers.time_step, dim,
                                           ParallelOption{}};
    }
    return ctrw::JumpGenerator_Diffusion{params_transport.diff_coeff, 0., dim,
                                         ParallelOption{}};
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Steppers\n"
           << io::line() << "Advection: " << advection_method << "\n"
           << "Diffusion: " << diffusion_method << "\n"
           << io::line();
    return output;
  }

  inline static std::string advection_method =
      "RK2"; /**< Name or description of advection method. */
  inline static std::string diffusion_method =
      "Stochastic Euler"; /**< Name or description of diffusion method. */
};

/**
   \struct Steppers_Advection_RK4_Diffusion_Euler PTOF/Steppers.h
   "PTOF/Steppers.h"
   \brief Time steppers for RK4 advection and stochastic forward Euler
   diffusion.
*/
struct Steppers_Advection_RK4_Diffusion_Euler {
  Steppers_Advection_RK4_Diffusion_Euler() = delete;

  /**
     \brief Adaptive time step generator.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
     \return Adaptive time step TimeGenerator.
  */
  template <typename Geometry, typename VelocityField, typename Boundary,
            typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  static auto makeTimeGenerator(Geometry const &geometry,
                                VelocityField const &velocity_field,
                                Boundary const &boundary,
                                TransportParameters const &params_transport,
                                ReactionParameters const &params_reaction,
                                SolverParameters const &params_solvers) {
    return TimeGenerator_Adaptive(geometry, velocity_field,
                                  boundary.surface_reaction, params_transport,
                                  params_reaction, params_solvers);
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters.
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Advection--diffusion JumpGenerator.
     \note
     - \p params_transport must define:
         - \c diff_coeff
     - If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
   */
  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto
  makeJumpGenerator(VelocityField &&velocity_field, Boundary &&boundary,
                    TransportParameters const &params_transport,
                    SolverParameters const &params_solvers, std::size_t dim) {
    return ctrw::JumpGenerator_Add{
        makeJumpGenerator_Advection<ParallelOption>(
            std::forward<VelocityField>(velocity_field),
            std::forward<Boundary>(boundary), params_transport, params_solvers,
            dim),
        makeJumpGenerator_Diffusion<ParallelOption>(
            std::forward<Boundary>(boundary), params_transport, params_solvers,
            dim)};
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure advection JumpGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto makeJumpGenerator_Advection(
      VelocityField &&velocity_field, Boundary &&boundary,
      TransportParameters const &params_transport,
      SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>) {
      return ctrw::JumpGenerator_Velocity_RK4{
          std::forward<VelocityField>(velocity_field), params_solvers.time_step,
          std::forward<Boundary>(boundary), ParallelOption{}};
    }
    return ctrw::JumpGenerator_Velocity_RK4{
        std::forward<VelocityField>(velocity_field), 0.,
        std::forward<Boundary>(boundary), ParallelOption{}};
  }

  /**
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure diffusion JumpGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename ParallelOption, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto makeJumpGenerator_Diffusion(
      Boundary &&boundary, TransportParameters const &params_transport,
      SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>) {
      return ctrw::JumpGenerator_Diffusion{params_transport.diff_coeff,
                                           params_solvers.time_step, dim,
                                           ParallelOption{}};
    }
    return ctrw::JumpGenerator_Diffusion{params_transport.diff_coeff, 0., dim,
                                         ParallelOption{}};
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Steppers\n"
           << io::line() << "Advection: " << advection_method << "\n"
           << "Diffusion: " << diffusion_method << "\n"
           << io::line();
    return output;
  }

  inline static std::string advection_method =
      "RK4"; /**< Name or description of advection method. */
  inline static std::string diffusion_method =
      "Stochastic Euler"; /**< Name or description of diffusion method. */
};

/**
   \struct Steppers_Advection_Heun_Diffusion_Euler PTOF/Steppers.h
   "PTOF/Steppers.h"
   \brief Time steppers for Heun advection and stochastic forward Euler
   diffusion.
*/
struct Steppers_Advection_Heun_Diffusion_Euler {
  Steppers_Advection_Heun_Diffusion_Euler() = delete;

  /**
     \brief Adaptive time step generator.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
     \return Adaptive time step TimeGenerator.
  */
  template <typename Geometry, typename VelocityField, typename Boundary,
            typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  static auto makeTimeGenerator(Geometry const &geometry,
                                VelocityField const &velocity_field,
                                Boundary const &boundary,
                                TransportParameters const &params_transport,
                                ReactionParameters const &params_reaction,
                                SolverParameters const &params_solvers) {
    return TimeGenerator_Adaptive(geometry, velocity_field,
                                  boundary.surface_reaction, params_transport,
                                  params_reaction, params_solvers);
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters.
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Advection--diffusion JumpGenerator.
     \note
     - \p params_transport must define:
         - \c diff_coeff
     - If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
   */
  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto
  makeJumpGenerator(VelocityField &&velocity_field, Boundary &&boundary,
                    TransportParameters const &params_transport,
                    SolverParameters const &params_solvers, std::size_t dim) {
    return ctrw::JumpGenerator_Add{
        makeJumpGenerator_Advection<ParallelOption>(
            std::forward<VelocityField>(velocity_field),
            std::forward<Boundary>(boundary), params_transport, params_solvers,
            dim),
        makeJumpGenerator_Diffusion<ParallelOption>(
            std::forward<Boundary>(boundary), params_transport, params_solvers,
            dim)};
  }

  /**
     \param velocity_field Velocity field as a function of state.
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure advection JumpGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto makeJumpGenerator_Advection(
      VelocityField &&velocity_field, Boundary &&boundary,
      TransportParameters const &params_transport,
      SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>) {
      return ctrw::JumpGenerator_Velocity_Heun{
          std::forward<VelocityField>(velocity_field), params_solvers.time_step,
          std::forward<Boundary>(boundary), ParallelOption{}};
    }
    return ctrw::JumpGenerator_Velocity_Heun{
        std::forward<VelocityField>(velocity_field), 0.,
        std::forward<Boundary>(boundary), ParallelOption{}};
  }

  /**
     \param boundary Boundary condition enforcer.
     \param params_transport Transport parameters (unused).
     \param params_solvers Solver parameters.
     \param dim Spatial dimension.
     \return Pure diffusion JumpGenerator.
     \note If \p params_solvers does not define \c time_step, the time step is
     initially set to zero.
  */
  template <typename ParallelOption, typename Boundary,
            typename TransportParameters, typename SolverParameters>
  static auto makeJumpGenerator_Diffusion(
      Boundary &&boundary, TransportParameters const &params_transport,
      SolverParameters const &params_solvers, std::size_t dim) {
    if constexpr (meta::has_time_step_v<SolverParameters>) {
      return ctrw::JumpGenerator_Diffusion{params_transport.diff_coeff,
                                           params_solvers.time_step, dim,
                                           ParallelOption{}};
    }
    return ctrw::JumpGenerator_Diffusion<ParallelOption>{
        params_transport.diff_coeff, 0., dim, ParallelOption{}};
  }

  /**
       \brief Output generic information about object.
       \param output Output stream.
    */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Steppers\n"
           << io::line() << "Advection: " << advection_method << "\n"
           << "Diffusion: " << diffusion_method << "\n"
           << io::line();
    return output;
  }

  inline static std::string advection_method =
      "Heun"; /**< Name or description of advection method. */
  inline static std::string diffusion_method =
      "Stochastic Euler"; /**< Name or description of diffusion method. */
};
} // namespace ptof

#endif /* PTOF_STEPPERS_H */
