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
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace ptof {
template <typename Stepper_Advection, typename Stepper_Diffusion,
          bool reaction_v, typename Stepper_CTRW>
struct SolverParameters_Generic {
  static constexpr bool advection =
      !std::is_same_v<Stepper_Advection, meta::Empty>;
  static constexpr bool diffusion =
      !std::is_same_v<Stepper_Diffusion, meta::Empty>;
  static constexpr bool reaction = reaction_v;

  std::size_t nr_particles;
  bool min_global_local;
  double local_step_factor_adv = std::numeric_limits<double>::infinity();
  double local_step_factor_diff = std::numeric_limits<double>::infinity();
  double local_step_factor_surf_react = std::numeric_limits<double>::infinity();
  double local_step_factor_temporal_adv =
      std::numeric_limits<double>::infinity();
  double global_step_factor_adv = std::numeric_limits<double>::infinity();
  double global_step_factor_diff = std::numeric_limits<double>::infinity();
  double global_step_factor_surf_react =
      std::numeric_limits<double>::infinity();
  double ctrw_time_step = 0.;

  template <typename Geometry, typename TransportParameters,
            typename ReactionParameters>
  SolverParameters_Generic(Directories const &directories,
                           std::string const &parameter_set_name,
                           Geometry const &geometry,
                           TransportParameters const &params_transport,
                           ReactionParameters const &params_reaction) {
    std::string filename = directories.dir_parameters + "/solvers_" +
                           parameter_set_name + ".param";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse number of particles", nr_particles);

    std::string min_or_max_global_local;
    split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse whether to use minimum or maximum of "
                       "local and global time step constraints",
             min_or_max_global_local);
    std::string for_min_or_max_global_local =
        std::string{"Minimum or maximum of local and global time step "
                    "constraints option "} +
        min_or_max_global_local + " : ";
    if (min_or_max_global_local != "min" && min_or_max_global_local != "max") {
      throw std::runtime_error{in_file + for_min_or_max_global_local +
                               "Expected min or max"};
    }
    min_global_local = min_or_max_global_local == "min" ? true : false;

    split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    param_index = 0;
    if constexpr (advection) {
      io::read_or_default(
          split_line, param_index, std::numeric_limits<double>::infinity(),
          in_file + "Could not parse local advective step accuracy",
          local_step_factor_adv);
    }
    if constexpr (diffusion) {
      io::read_or_default(
          split_line, param_index, std::numeric_limits<double>::infinity(),
          in_file + "Could not parse local diffusive step accuracy factor",
          local_step_factor_diff);
    }
    if constexpr (reaction) {
      io::read_or_default(
          split_line, param_index, std::numeric_limits<double>::infinity(),
          in_file + "Could not parse local reactive step accuracy factor",
          local_step_factor_surf_react);
    }
    if constexpr (advection) {
      io::read_or_default(
          split_line, param_index, std::numeric_limits<double>::infinity(),
          in_file +
              "Could not parse local temporal advective step accuracy factor",
          local_step_factor_temporal_adv);
    }

    split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    param_index = 0;
    if constexpr (advection) {
      io::read_or_default(
          split_line, param_index, std::numeric_limits<double>::infinity(),
          in_file + "Could not parse global advective step accuracy factor",
          global_step_factor_adv);
    }
    if constexpr (diffusion) {
      io::read_or_default(
          split_line, param_index, std::numeric_limits<double>::infinity(),
          in_file + "Could not parse global diffusive step accuracy factor",
          global_step_factor_diff);
    }
    if constexpr (reaction) {
      io::read_or_default(
          split_line, param_index, std::numeric_limits<double>::infinity(),
          in_file + "Could not parse global reactive step accuracy factor",
          global_step_factor_surf_react);
    }

    check_constraints(in_file);

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

  void check_constraints(std::string const &in_file) {
    if (!(local_step_factor_adv >= 0. && local_step_factor_diff >= 0. &&
          local_step_factor_surf_react >= 0. &&
          local_step_factor_temporal_adv >= 0. &&
          global_step_factor_adv >= 0. && global_step_factor_diff >= 0. &&
          global_step_factor_surf_react >= 0.)) {
      throw std::runtime_error{in_file +
                               "Time step constraints should be non-negative"};
    }

    bool local_constraints_are_zero =
        (local_step_factor_adv == 0. || local_step_factor_diff == 0. ||
         local_step_factor_surf_react == 0.);
    bool local_transport_constraints_are_inf =
        (local_step_factor_adv == std::numeric_limits<double>::infinity() &&
         local_step_factor_diff == std::numeric_limits<double>::infinity() &&
         std::numeric_limits<double>::infinity());

    bool global_constraints_are_zero =
        (global_step_factor_adv == 0. || global_step_factor_diff == 0. ||
         global_step_factor_surf_react == 0.);
    bool global_transport_constraints_are_inf =
        (global_step_factor_adv == std::numeric_limits<double>::infinity() &&
         global_step_factor_diff == std::numeric_limits<double>::infinity() &&
         std::numeric_limits<double>::infinity());

    if (local_constraints_are_zero && global_constraints_are_zero) {
      throw std::runtime_error{in_file +
                               "Local and global transport-related time "
                               "step constraints are zero"};
    }
    if (local_transport_constraints_are_inf &&
        global_transport_constraints_are_inf) {
      throw std::runtime_error{in_file +
                               "All local and global transport-related time "
                               "step constraints are infinite"};
    }
    if (min_global_local &&
        (local_constraints_are_zero || global_constraints_are_zero)) {
      throw std::runtime_error{
          in_file + "Either local or global transport-related time "
                    "step constraints are zero and minimum is required"};
    }
    if (!min_global_local && (local_transport_constraints_are_inf ||
                              global_transport_constraints_are_inf)) {
      throw std::runtime_error{
          in_file + "Either all local or all global transport-related time "
                    "step constraints are infinite and maximum is required"};
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
           "- How to combine local and global time step constraints\n"
           "  - min\n"
           "    - Use minimum of local and global constraints\n"
           "  - max\n"
           "    - Use maximum of local and global constraints\n"
           "- Local step accuracy:\n"
           "  (Note:\n"
           "    - Pass inf to deactivate specific constraints; pass at least\n"
           "      one 0 to deactivate all local constraints with maximum)\n"
           "  - Pass on same line:\n"
           "    - Local step accuracy factor with respect to advection [inf]\n"
           "    - Local step accuracy factor with respect to diffusion [inf]\n"
           "    - Local step accuracy factor with respect to reaction [inf]\n"
           "    - Local step accuracy with respect to temporal advection\n"
           "      variability [inf]\n"
           "- Global step accuracy:\n"
           "  (Note:\n"
           "    - Initial values, e.g., of flow are used\n"
           "    - Pass inf to deactivate specific constraints; pass at least\n"
           "      one 0 to deactivate all global constraints with maximum)\n"
           "  - Pass on same line:\n"
           "    - Global step accuracy factor with respect to advection [inf]\n"
           "    - Global step accuracy factor with respect to diffusion [inf]\n"
           "    - Global step accuracy factor with respect to reaction [inf]\n";
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
