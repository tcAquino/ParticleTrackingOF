/**
   \file PTOF/SolverParameters.h
   \author Tomas Aquino
   \date 19/01/2025
   \brief Solver parameters and utilities.
*/

#ifndef PTOF_SOLVERPARAMETERS_H
#define PTOF_SOLVERPARAMETERS_H

#include "General/IO.h"
#include "General/Meta.h"
#include "PTOF/Directories.h"
#include "PTOF/Steppers.h"
#include "PTOF/TimeUnitsList.h"
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

    auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse number of particles", nr_particles);

    split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
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
      } else {
        io::read(split_line, param_index,
                 in_file + "Could not parse local advective time step accuracy",
                 local_time_step_adv);
      }
    } else {
      if (split_line.size() >= 2) {
        io::read(split_line, param_index,
                 in_file + "Could not parse local diffusive and reactive "
                           "time step accuracy",
                 local_time_step_diff, local_time_step_react);
      } else {
        io::read(split_line, param_index,
                 in_file + "Could not parse local diffusive time step accuracy",
                 local_time_step_diff);
      }
    }

    split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
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
                 global_time_step_adv, global_time_step_react);
      } else {
        io::read(split_line, param_index,
                 in_file +
                     "Could not parse global diffusive time step accuracy",
                 global_time_step_adv);
      }
    } else {
      if (split_line.size() >= 2) {
        io::read(split_line, param_index,
                 in_file + "Could not parse global diffusive and reactive time "
                           "step accuracy",
                 global_time_step_diff, global_time_step_react);
      } else {
        io::read(split_line, param_index,
                 in_file +
                     "Could not parse global diffusive time step accuracy",
                 global_time_step_diff);
      }
    }

    if constexpr (std::is_same_v<Stepper_CTRW, CTRWSteppers::TimeStep>) {
      split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
      param_index = 0;
      auto time_units = io::read<std::string>(
          split_line, param_index, in_file + "Could not parse time units");
      if (!TimeUnitsList::contains(time_units)) {
        throw std::runtime_error{in_file + "Not supported"};
      }
      double time_unit_factor =
          ptof::time_unit_factor(time_units, params_transport, params_reaction);
      io::read(split_line, param_index,
               in_file + "Could not parse time step for synchronizing CTRW",
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
           "    - Pass inf to deactivate specific constraints; pass at least\n"
           "      one 0 to deactivate all global constraints\n"
           "    - At least one local or one global constraint must be active)\n"
           "  - Pass on same line:\n"
           "    - Time step accuracy with respect to local advection time\n"
           "      (pass only if there is advection)\n"
           "    - Time step accuracy with respect to local diffusion time\n"
           "      (pass only if there is diffusion)\n"
           "    - Time step accuracy with respect to local reaction time\n"
           "      (optional)\n"
           "- Global time step accuracy:\n"
           "  (Note:\n"
           "    - Minimum between processes is used\n"
           "    - Initial values (e.g., of flow) are used\n"
           "    - Maximum between local and global is used\n"
           "    - Pass inf to deactivate specific constraints; pass at least\n"
           "      one 0 to deactivate all global constraints\n"
           "    - At least one local or one global constraint must be active)\n"
           "  - Pass on same line:\n"
           "    - Time step accuracy with respect to global advection time\n"
           "      (pass only if there is advection)\n"
           "    - Time step accuracy with respect to global diffusion time\n"
           "      (pass only if there is diffusion)\n"
           "    - Time step accuracy with respect to global reaction time\n"
           "      (optional)\n";
    if constexpr (std::is_same_v<Stepper_CTRW, CTRWSteppers::TimeStep>) {
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
} // namespace ptof

#endif /* PTOF_SOLVERPARAMETERS */
