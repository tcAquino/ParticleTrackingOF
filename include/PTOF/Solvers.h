/**
   \file PTOF/Solvers.h
   \author Tomás Aquino
   \date 19/01/2025
   \brief Solver parameters and utilities.
*/

#ifndef PTOF_SOLVERS_H
#define PTOF_SOLVERS_H

#include "General/IO.h"
#include "General/Meta.h"
#include "PTOF/Directories.h"
#include "PTOF/Steppers.h"
#include "PTOF/TimeUnits.h"
#include <ostream>
#include <string>
#include <type_traits>
#include <utility>

namespace ptof {
template <typename Stepper_Advection, typename Stepper_Diffusion,
          typename Stepper_CTRW>
struct SolverParameters_Generic {
  static constexpr bool advection =
      !std::is_same_v<Stepper_Advection, meta::Empty>;
  static constexpr bool diffusion =
      !std::is_same_v<Stepper_Diffusion, meta::Empty>;

  std::size_t nr_particles;
  double local_time_step_adv = std::numeric_limits<double>::infinity();
  double local_time_step_diff = std::numeric_limits<double>::infinity();
  double local_time_step_react = std::numeric_limits<double>::infinity();
  double global_time_step_adv = std::numeric_limits<double>::infinity();
  double global_time_step_diff = std::numeric_limits<double>::infinity();
  double global_time_step_react = std::numeric_limits<double>::infinity();
  double ctrw_time_step = 0.;

  template <typename Geometry, typename TransportParameters,
            typename ReactionParameters>
  SolverParameters_Generic(
      Directories const &directories, std::string const &parameter_set_name,
      Geometry const &geometry, TransportParameters const &params_transport,
      ReactionParameters const &params_reaction = meta::Empty{}) {
    std::string filename = directories.dir_parameters + "/parameters_solvers_" +
                           parameter_set_name + ".dat";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse number of particles", nr_particles);

    split_line = io::split_line(input);
    param_index = 0;
    if constexpr (advection && diffusion) {
      if (split_line.size() >= 3) {
        io::read(split_line, param_index,
                 in_file + "Could not parse local advective, diffusive, and "
                           "reactive time step accuracy",
                 local_time_step_adv, local_time_step_diff,
                 local_time_step_react);
      } else {
        io::read(split_line, param_index,
                 in_file + "Could not parse local advective and diffusive time "
                           "step accuracy",
                 local_time_step_adv, local_time_step_diff);
      }
    } else if constexpr (advection) {
      if (split_line.size() >= 2) {
        io::read(split_line, param_index,
                 in_file + "Could not parse local advective and reactive "
                           "time step accuracy",
                 local_time_step_adv, local_time_step_react);
      } else if constexpr (diffusion) {
        io::read(split_line, param_index,
                 in_file + "Could not parse local diffusive time step accuracy",
                 local_time_step_diff);
      }
    }

    split_line = io::split_line(input);
    param_index = 0;
    if constexpr (advection && diffusion) {
      if (split_line.size() >= 3) {
        io::read(split_line, param_index,
                 in_file + "Could not parse global advective, diffusive, and "
                           "reactive time step accuracy",
                 global_time_step_adv, global_time_step_diff,
                 global_time_step_react);
      } else {
        io::read(split_line, param_index,
                 in_file + "Could not parse global advective and diffusive "
                           "time step accuracy",
                 global_time_step_adv, global_time_step_diff);
      }
    } else if constexpr (advection) {
      if (split_line.size() >= 2) {
        io::read(split_line, param_index,
                 in_file + "Could not parse global advective and reactive time "
                           "step accuracy",
                 global_time_step_adv, global_time_step_adv);
      } else if constexpr (diffusion) {
        io::read(split_line, param_index,
                 in_file +
                     "Could not parse global diffusive time step accuracy",
                 global_time_step_diff, global_time_step_react);
      }
    }

    if constexpr (std::is_same_v<Stepper_CTRW, CTRWStepper::TimeStep>) {
      split_line = io::split_line(input);
      param_index = 0;
      auto time_units = io::read<std::string>(
          split_line, param_index, in_file + "Could not parse time units");
      if (!TimeUnits{}.contains(time_units))
        throw std::runtime_error{in_file + "Not supported"};
      double time_unit_factor =
          ptof::time_unit_factor(time_units, params_transport, params_reaction);
      io::read(split_line, param_index,
               in_file + "Could not parse Time step for synchronizing CTRW",
               ctrw_time_step);
      ctrw_time_step *= time_unit_factor;
    }
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output
        << io::line() << "Solver parameters\n"
        << io::line()
        << "- Number of Lagrangian particles in each injection step\n"
           "- Local time step accuracy:\n"
           "  (Note:\n"
           "    - Minimum between processes is used\n"
           "    - Pass inf to deactivate specific constraints; Pass zeros\n"
           "      deactivate all local constraints\n"
           "    - At least one local or one global constraint must be active,\n"
           "      otherwise the time step will be zero or infinity\n"
           "  - Pass on same line:\n"
           "    - Time step accuracy with respect to local advection time\n"
           "    - Time step accuracy with respect to local diffusion time\n"
           "    - Time step accuracy with respect to local reaction time\n"
           "      (optional)\n"
           "- Global time step accuracy:\n"
           "  (Note:\n"
           "    - Minimum between processes is used\n"
           "    - Initial values (e.g., of flow) are used\n"
           "    - Maximum between local and global is used\n"
           "    - Inf deactivates specific constraints; at least one 0\n"
           "      deactivates all global constraints\n"
           "    - At least one local or one global constraint must be active\n"
           "  - Pass on same line:\n"
           "    - Time step accuracy with respect to global advection\n"
           "      time\n"
           "    - Time step accuracy with respect to global diffusion\n"
           "      time\n"
           "    - Time step accuracy with respect to global reaction time\n"
           "      (optional)\n";
    if constexpr (std::is_same_v<Stepper_CTRW, CTRWStepper::TimeStep>) {
      output << "- Time units for CTRW synchronization time step:\n"
                "  - diffusion\n"
                "    - Diffusion time units\n"
                "  - advection\n"
                "    - Advection time units\n"
                "  - reaction\n"
                "    - Reaction time units\n"
                "  - arbitrary\n"
                "    - Arbitary units (no rescaling)\n"
                "  - Pass on same line:\n"
                "    - Time step for synchronizing CTRW\n";
    }
    output << io::line();
    return output;
  }
};

template <typename Stepper_Advection_t, typename Stepper_Diffusion_t,
          typename Stepper_CTRW_t>
struct Solvers_Generic {
  Solvers_Generic() = delete;

  using Stepper_Advection = Stepper_Advection_t;
  using Stepper_Diffusion = Stepper_Diffusion_t;
  using Stepper_CTRW = Stepper_CTRW_t;
  using Parameters = SolverParameters_Generic<Stepper_Advection,
                                              Stepper_Diffusion, Stepper_CTRW>;
  static_assert(Parameters::advection || Parameters::diffusion,
                "No advection or diffusion must be present");
  static_assert(!(Parameters::diffusion &&
                  !std::is_same_v<Stepper_Diffusion, Stepper::Euler>),
                "Currently only no diffusion or stochastic forward Euler "
                "stepping are supported for diffusion");
  static_assert(
      !(Parameters::advection &&
        !std::is_same_v<Stepper_Advection, Stepper::Euler> &&
        !std::is_same_v<Stepper_Advection, Stepper::RK2> &&
        !std::is_same_v<Stepper_Advection, Stepper::RK4> &&
        !std::is_same_v<Stepper_Advection, Stepper::Heun>),
      "Currently only no advection, forward Euler, RK2, RK4, or Heun stepping "
      "are supported for advection");
  using Steppers = std::conditional_t<
      std::is_same_v<Stepper_Advection, Stepper::Euler>,
      Steppers_Advection_Euler_Diffusion_Euler,
      std::conditional_t<
          std::is_same_v<Stepper_Advection, Stepper::RK2>,
          Steppers_Advection_RK2_Diffusion_Euler,
          std::conditional_t<std::is_same_v<Stepper_Advection, Stepper::RK4>,
                             Steppers_Advection_RK4_Diffusion_Euler,
                             Steppers_Advection_Heun_Diffusion_Euler>>>;

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Solvers\n" << io::line();
    if constexpr (Parameters::advection) {
      output << "Advection: " << Steppers::advection_method << "\n";
    }
    if constexpr (Parameters::diffusion) {
      output << "Diffusion: " << Steppers::diffusion_method << "\n";
    }
    output << "CTRW: " << Stepper_CTRW::method << "\n";
    output << io::line();
    return output;
  }

  template <typename CTRW, typename Transitions, typename Time,
            typename Update = meta::DoNothing>
  static void evolve(CTRW &ctrw, Transitions &&transitions,
                     Parameters const &parameters_solvers, Time new_time,
                     Time old_time, Update &&update = {}) {
    return Stepper_CTRW::evolve(ctrw, transitions, parameters_solvers, new_time,
                                old_time, update);
  }

  static auto makeTimeGenerator(Parameters const &params_solvers) {
    return Steppers::makeTimeGenerator(params_solvers);
  }

  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters>
  static auto
  makeJumpGenerator(VelocityField &&velocity_field, Boundary &&boundary,
                    TransportParameters const &params_transport,
                    Parameters const &params_solvers, std::size_t dim) {
    if constexpr (Parameters::advection && Parameters::diffusion)
      return Steppers::template makeJumpGenerator<ParallelOption>(
          std::forward<VelocityField>(velocity_field),
          std::forward<Boundary>(boundary), params_transport, params_solvers,
          dim);
    else if constexpr (Parameters::advection)
      return Steppers::makeJumpGenerator_Advection(
          std::forward<VelocityField>(velocity_field),
          std::forward<Boundary>(boundary), params_transport, params_solvers,
          dim);
    else if constexpr (Parameters::diffusion)
      return Steppers::template makeJumpGenerator_Diffusion<ParallelOption>(
          std::forward<Boundary>(boundary), params_transport, params_solvers,
          dim);
  }
};
} // namespace ptof

#endif /* PTOF_SOLVERS_H */
