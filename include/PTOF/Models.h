/**
 \file PTOF/Models.h
 \author Tomás Aquino
 \date 22/02/2022
 */

#ifndef PTOF_MODELS_H
#define PTOF_MODELS_H

#include "CTRW/CTRW.h"
#include "CTRW/Meta.h"
#include "General/Meta.h"
#include "General/Serial.h"
#include "General/Useful.h"
#include "PTOF/Advection.h"
#include "PTOF/Geometry.h"
#include "PTOF/Info.h"
#include "PTOF/InitialCondition.h"
#include "PTOF/Locator.h"
#include "PTOF/Output.h"
#include "PTOF/ParticleMaker.h"
#include "PTOF/Reaction.h"
#include "PTOF/State.h"
#include "PTOF/Steppers.h"
#include "PTOF/Transitions.h"
#include "PTOF/Useful.h"
#include <cmath>
#include <cstddef>
#include <limits>
#include <meshSearch.H>
#include <stdexcept>
#include <string>
#include <type_traits>

/** \namespace ptof Objects and methods for ParticleTrackingOF. */
namespace ptof {
/** \namespace ptof::model_advection_diffusion_2d
 Definitions for 2D advective--diffusive transport. */
namespace model_advection_diffusion_2d {
struct Model {
  Model() = delete;

  inline static const std::string name{"advection_diffusion_2d"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------------------\n"
           << name << "\n"
           << "--------------------------------------------------------------\n"
              "\n";
  }
};

using Geometry = Geometry<2, meta::ParallelOptions::Serial>;
using Info = ptof::Info_Absorbed;
using State = State<Geometry::dim, Info, double, double, std::size_t>;
using CTRW = ctrw::CTRW<State>;

struct Solvers {
  Solvers() = delete;

  using Steppers = Steppers_Advection_Euler_Diffusion_Euler;

  struct Parameters {
  private:
    std::string comment_sequence =
        "#"; /**< Sequence of characters marking comment for file parsing. */

  public:
    std::size_t nr_particles;
    double local_time_step_adv;
    double local_time_step_diff;
    double local_time_step_react;
    double global_time_step_adv;
    double global_time_step_diff;
    double global_time_step_react;

    template <typename Geometry, typename TransportParameters,
              typename ReactionParameters>
    Parameters(Directories const &directories, std::string const &name,
               Geometry const &geometry,
               TransportParameters const &params_transport,
               ReactionParameters const &params_reaction) {
      auto input = useful::open_read(directories.dir_parameters +
                                     "/parameters_solvers_" + name + ".dat");
      useful::read_first_from_line(input, nr_particles, comment_sequence);
      useful::read_first_from_line(input, local_time_step_adv,
                                   comment_sequence);
      useful::read_first_from_line(input, local_time_step_diff,
                                   comment_sequence);
      useful::read_first_from_line(input, local_time_step_react,
                                   comment_sequence);
      useful::read_first_from_line(input, global_time_step_adv,
                                   comment_sequence);
      useful::read_first_from_line(input, global_time_step_diff,
                                   comment_sequence);
      useful::read_first_from_line(input, global_time_step_react,
                                   comment_sequence);
      input.close();
    }

    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Solver parameters\n"
             "--------------------------------------------------------------\n"
             "- Number of Lagrangian particles in each injection step\n"
             "- Local time step in terms of cell-based advection time\n"
             "- Local time step in terms of cell-based diffusion time\n"
             "- Local time step in terms of cell-based surface reaction time\n"
             "- Global time step in terms of characteristic advection time\n"
             "- Global time step in terms of characteristic diffusion time\n"
             "- Global time step in terms of surface reaction time\n"
             "  (Note: Minimum between processes and maximum between local\n"
             "         and global is used.\n"
             "         Initial values (e.g., of flow) are used for global\n"
             "         quantities.)\n"
             "--------------------------------------------------------------\n";
    }
  };

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Solvers\n"
              "--------------------------------------------------\n"
              "Advection: Euler\n"
              "Diffusion: Stochastic Euler\n"
              "--------------------------------------------------\n";
  }
};

struct Transport {
  Transport() = delete;

  struct Parameters {
  private:
    std::string comment_sequence =
        "#"; /**< Sequence of characters marking comment for file parsing. */

  public:
    std::string peclet_option;
    double lengthscale;
    double diff_coeff;
    double diffusion_time;
    double advection_time;
    double peclet;
    double mean_velocity;
    double velocity_rescaling_factor;

    template <typename Geometry>
    Parameters(Directories const &directories, std::string const &name,
               Geometry const &geometry) {
      auto input = useful::open_read(directories.dir_parameters +
                                     "/parameters_transport_" + name + ".dat");
      useful::read_first_from_line(input, peclet_option, comment_sequence);
      useful::read_first_from_line(input, lengthscale, comment_sequence);
      if (peclet_option == "rescale_velocity_to_peclet") {
        useful::read_first_from_line(input, peclet, comment_sequence);
        useful::read_first_from_line(input, diff_coeff, comment_sequence);
        diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
        advection_time = 2. * diffusion_time / peclet;
        mean_velocity = lengthscale / advection_time;
      } else if (peclet_option == "rescale_velocity_to_mean") {
        useful::read_first_from_line(input, mean_velocity, comment_sequence);
        useful::read_first_from_line(input, diff_coeff, comment_sequence);
        diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
        peclet = lengthscale * mean_velocity / diff_coeff;
        advection_time = lengthscale / mean_velocity;
      } else if (peclet_option == "rescale_velocity_to_advection_time") {
        useful::read_first_from_line(input, advection_time, comment_sequence);
        useful::read_first_from_line(input, diff_coeff, comment_sequence);
        diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
        peclet = 2. * diffusion_time / advection_time;
        mean_velocity = lengthscale / advection_time;
      } else if (peclet_option == "compute_from_diff_coeff") {
        useful::read_first_from_line(input, diff_coeff, comment_sequence);
        diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      } else if (peclet_option == "compute_from_diff_time") {
        useful::read_first_from_line(input, diffusion_time, comment_sequence);
        diff_coeff = lengthscale * lengthscale / (2. * diffusion_time);
      } else if (peclet_option == "set_diff_coeff") {
        useful::read_first_from_line(input, peclet, comment_sequence);
      } else
        throw std::runtime_error{"Peclet number setting option " +
                                 peclet_option + " not supported"};
      input.close();
    }

    template <typename VelocityField, typename Mesh>
    void rescale(VelocityField &velocity_field, Mesh const &mesh) {
      double current_mean =
          ptof::magnitude_of_average(velocity_field.field(), mesh);
      if (peclet_option == "rescale_velocity_to_peclet" ||
          peclet_option == "rescale_velocity_to_mean" ||
          peclet_option == "rescale_velocity_to_advection_time") {
        velocity_rescaling_factor = mean_velocity / current_mean;
        velocity_field.rescale(velocity_rescaling_factor);
      } else if (peclet_option == "compute_from_diff_coeff" ||
                 peclet_option == "compute_from_diff_time") {
        velocity_rescaling_factor = 1.;
        mean_velocity = current_mean;
        advection_time = lengthscale / mean_velocity;
        peclet = 2. * diffusion_time / advection_time;
      } else if (peclet_option == "set_diff_coeff") {
        velocity_rescaling_factor = 1.;
        mean_velocity = current_mean;
        diff_coeff = lengthscale * mean_velocity / peclet;
        diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
        advection_time = lengthscale / mean_velocity;
      } else
        throw std::runtime_error{"Peclet number setting option " +
                                 peclet_option + " not supported"};
    }

    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Transport parameters\n"
             "--------------------------------------------------------------\n"
             "- How to set the Peclet number:\n"
             "  (Note: Rescaling of the velocity field is not handled when\n"
             "         transport parameters are set.\n"
             "         Rescale method must always be called, both to\n"
             "         rescale and to compute\n"
             "         velocity-dependent quantities when not rescaling)\n"
             "\tcompute_from_diff_coeff: Compute from given diffusion\n"
             "\t                         coefficient (do not rescale velocity\n"
             "\t                         field)\n"
             "\tset_diff_coeff: Compute diffusion coefficient to set given\n"
             "\t                Peclet number (do not rescale velocity field)\n"
             "\tcompute_from_diff_time: Compute from given diffusion time (do\n"
             "\t                        not rescale velocity field)\n"
             "\trescale_velocity_to_peclet: Rescale velocity field according\n"
             "\t                            to imposed peclet number\n"
             "\trescale_velocity_to_mean: Rescale according to imposed mean\n"
             "\t                          flow velocity\n"
             "\trescale_velocity_to_advection_time: Rescale according to\n"
             "\t                                    imposed advection time\n"
             "- Reference lengthscale\n"
             "- Peclet number (do not pass if Peclet option is\n"
             "  compute_from_diff_coeff or compute_from_diff_time)\n"
             "- Mean flow velocity (pass only if Peclet option is\n"
             "  rescale_velocity_to_mean)\n"
             "- Advection time (pass only if Peclet option is\n"
             "  rescale_velocity_to_advection_time)\n"
             "- Diffusion time (pass only if Peclet option is\n"
             "  compute_from_diff_time)\n"
             "- Diffusion coefficient (do not pass if Peclet option is\n"
             "  set_diff_coeff or compute_from_diff_time)\n"
             "--------------------------------------------------------------\n";
    }
  };

  template <typename Solvers, typename VelocityField, typename Geometry,
            typename Boundary, typename ReactionParameters>
  static auto
  makeTransitions(VelocityField const &velocity_field, Geometry const &geometry,
                  Boundary &boundary, Parameters const &params_transport,
                  ReactionParameters const &params_reaction,
                  typename Solvers::Parameters const &params_solvers) {
    return makeTransportTransitions<typename Solvers::Steppers>(
        velocity_field, geometry, boundary, params_transport, params_reaction,
        params_solvers);
  }

  template <typename Geometry, typename TransportParameters>
  static auto makeVelocityInterpolator(Geometry const &geometry,
                                       TransportParameters &params_transport) {
    return makeLinearVelocityInterpolator(
        geometry, ptof::get_velocity_data(geometry.mesh()), params_transport);
  }

  template <typename Geometry, typename Uninterpolated,
            typename TransportParameters>
  static auto makeVelocityInterpolator(Geometry const &geometry,
                                       Uninterpolated &&uninterpolated,
                                       TransportParameters &params_transport) {
    return makeLinearVelocityInterpolator(
        geometry, ptof::get_velocity_data(geometry.mesh()),
        std::forward<Uninterpolated>(uninterpolated), params_transport);
  }

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Transport\n"
           "--------------------------------------------------------------\n"
           "Processes: Advection-diffusion\n"
           "Interpolation: Linear\n"
           "--------------------------------------------------------------\n";
  }
};

struct Reaction {
  Reaction() = delete;

  using BulkReaction = BulkReaction_DoNothing;
  using SurfaceReaction = SurfaceReaction_DoNothing;

  struct Parameters {
    double damkohler{0.};
    double rate_constant{0.};
    double reaction_time{std::numeric_limits<double>::infinity()};

    Parameters() {}

    template <typename Geometry, typename TransportParameters>
    Parameters(Directories const &dir, std::string const &name,
               Geometry const &geometry,
               TransportParameters const &params_transport) {}

    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Reaction parameters\n"
             "--------------------------------------------------------------\n"
             "None\n"
             "--------------------------------------------------------------\n";
    }
  };

  template <typename Geometry, typename TransportParameters,
            typename SolverParameters>
  static auto makeBulkReaction(Geometry const &geometry,
                               Parameters const &parameters,
                               TransportParameters const &params_transport,
                               SolverParameters const &params_solvers) {
    return BulkReaction{};
  }

  template <typename Geometry, typename TransportParameters,
            typename SolverParameters>
  static auto makeSurfaceReaction(Geometry const &geometry,
                                  Parameters const &parameters,
                                  TransportParameters const &params_transport,
                                  SolverParameters const &params_solvers) {
    return SurfaceReaction{};
  }

  template <typename OStream> static void info(OStream &output) {
    BulkReaction::info(output);
    output << "\n";
    SurfaceReaction::info(output);
  }
};

struct InitialCondition {
  InitialCondition() = delete;

  using Parameters = InitialConditionParameters_Cases;

  template <typename CTRW, typename Geometry, typename VelocityField,
            typename SolverParameters, typename Mask = useful::Empty>
  static auto makeInitialCondition(Geometry const &geometry,
                                   VelocityField const &velocity_field,
                                   Parameters const &parameters,
                                   SolverParameters const &params_solvers,
                                   Mask &&mask = {}, double threshold = 0.) {
    if constexpr (meta::has_periodicity_v<typename CTRW::State>)
      return InitialCondition_Cases{
          geometry,
          velocity_field,
          ParticleMaker_Periodic{
              meta::Selector_t<typename CTRW::Particle>{}, geometry.locator,
              geometry.boundary_periodic, parameters.time_min,
              parameters.initial_mass / params_solvers.nr_particles},
          params_solvers.nr_particles,
          parameters,
          std::forward<Mask>(mask),
          threshold};
    else
      return InitialCondition_Cases{
          geometry,
          velocity_field,
          ParticleMaker{meta::Selector_t<typename CTRW::Particle>{},
                        geometry.locator, parameters.time_min,
                        parameters.initial_mass / params_solvers.nr_particles},
          params_solvers.nr_particles,
          parameters,
          std::forward<Mask>(mask),
          threshold};
  }
};

struct Output {
  Output() = delete;

  using Parameters = OutputParameters;

  template <typename Subject, typename VelocityField, typename Geometry,
            typename Mask = useful::Empty>
  static auto
  makeOutput(Subject const &subject, VelocityField const &velocity_field,
             Geometry const &geometry, Directories const &directories,
             Parameters parameters, std::string const &identifier,
             std::vector<std::reference_wrapper<const Mask>> masks = {},
             std::vector<double> thresholds = {}) {
    return Output_Cases{subject,    velocity_field, geometry, directories,
                        parameters, identifier,     masks,    thresholds};
  }

  template <typename Subject, typename VelocityField, typename Geometry,
            typename Mask>
  static auto
  makeOutput(Subject const &subject, VelocityField const &velocity_field,
             Geometry const &geometry, Directories const &directories,
             Parameters parameters, std::string const &identifier,
             std::initializer_list<std::reference_wrapper<const Mask>> masks,
             std::initializer_list<double> thresholds = {}) {
    return Output_Cases{subject,    velocity_field, geometry, directories,
                        parameters, identifier,     masks,    thresholds};
  }

  template <typename Subject, typename VelocityField, typename Geometry,
            typename Mask>
  static auto
  makeOutput(Subject const &subject, VelocityField const &velocity_field,
             Geometry const &geometry, Directories const &directories,
             Parameters parameters, std::string const &identifier,
             std::vector<std::reference_wrapper<const Mask>> masks,
             std::initializer_list<double> thresholds = {}) {
    return Output_Cases{subject,    velocity_field, geometry, directories,
                        parameters, identifier,     masks,    thresholds};
  }

  template <typename Subject, typename VelocityField, typename Geometry,
            typename Mask>
  static auto
  makeOutput(Subject const &subject, VelocityField const &velocity_field,
             Geometry const &geometry, Directories const &directories,
             Parameters parameters, std::string const &identifier,
             std::initializer_list<std::reference_wrapper<const Mask>> masks,
             std::vector<double> thresholds = {}) {
    return Output_Cases{subject,    velocity_field, geometry, directories,
                        parameters, identifier,     masks,    thresholds};
  }
};

template <typename Transport, typename Solvers, typename Reaction,
          typename VelocityField, typename Geometry, typename Boundary,
          std::enable_if_t<std::is_same_v<typename Geometry::ParallelOption,
                                          meta::ParallelOptions::Serial>,
                           std::nullptr_t> = nullptr>
auto makeTransitions(VelocityField const &velocity_field,
                     Geometry const &geometry, Boundary &boundary,
                     typename Transport::Parameters const &params_transport,
                     typename Reaction::Parameters const &params_reaction,
                     typename Solvers::Parameters const &params_solvers) {
  return ctrw::Transitions_CTRW_Transport_Reaction{
      Transport::template makeTransitions<Solvers>(
          velocity_field, geometry, boundary, params_transport, params_reaction,
          params_solvers),
      Reaction::makeBulkReaction(geometry, params_reaction, params_transport,
                                 params_solvers)};
}

template <typename Transport, typename Solvers, typename VelocityField,
          typename Geometry, typename Boundary, typename ReactionParameters,
          typename BulkReaction,
          std::enable_if_t<std::is_same_v<typename Geometry::ParallelOption,
                                          meta::ParallelOptions::Serial>,
                           std::nullptr_t> = nullptr>
auto makeTransitions(VelocityField const &velocity_field,
                     Geometry const &geometry, Boundary &boundary,
                     typename Transport::Parameters const &params_transport,
                     ReactionParameters const &params_reaction,
                     typename Solvers::Parameters const &params_solvers,
                     BulkReaction const &bulk_reaction) {
  return ctrw::Transitions_CTRW_Transport_Reaction{
      Transport::template makeTransitions<Solvers>(
          velocity_field, geometry, boundary, params_transport, params_reaction,
          params_solvers),
      bulk_reaction};
}
} // namespace model_advection_diffusion_2d

/** \namespace ptof::model_advection_2d
 Definitions for 2D advective transport. */
namespace model_advection_2d {
struct Model {
  Model() = delete;

  inline static const std::string name{"advection_2d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d::CTRW;
using model_advection_diffusion_2d::Geometry;
using model_advection_diffusion_2d::Info;
using model_advection_diffusion_2d::InitialCondition;
using model_advection_diffusion_2d::makeTransitions;
using model_advection_diffusion_2d::Output;
using model_advection_diffusion_2d::Reaction;
using model_advection_diffusion_2d::State;

struct Solvers {
  Solvers() = delete;

  using Steppers = Steppers_Advection_Euler_Diffusion_Euler;

  struct Parameters {
  private:
    std::string comment_sequence =
        "#"; /**< Sequence of characters marking comment for file parsing. */

  public:
    std::size_t nr_particles;
    double local_time_step_adv;
    double local_time_step_diff = std::numeric_limits<double>::infinity();
    double local_time_step_react = std::numeric_limits<double>::infinity();
    double global_time_step_adv;
    double global_time_step_diff = std::numeric_limits<double>::infinity();
    double global_time_step_react = std::numeric_limits<double>::infinity();

    template <typename Geometry, typename TransportParameters,
              typename ReactionParameters>
    Parameters(Directories const &directories, std::string const &name,
               Geometry const &geometry,
               TransportParameters const &params_transport,
               ReactionParameters const &params_reaction) {
      auto input = useful::open_read(directories.dir_parameters +
                                     "/parameters_solvers_" + name + ".dat");
      useful::read_first_from_line(input, nr_particles, comment_sequence);
      useful::read_first_from_line(input, local_time_step_adv,
                                   comment_sequence);
      useful::read_first_from_line(input, global_time_step_adv,
                                   comment_sequence);
      input.close();
    }

    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Solver parameters\n"
             "--------------------------------------------------------------\n"
             "- Number of Lagrangian particles in each injection step\n"
             "- Local timestep accuracy in terms of cell-based advection time\n"
             "- Global timestep accuracy in terms of characteristic advection\n"
             "  time\n"
             "  (Note: Minimum between processes and maximum between local\n"
             "         and global is used.\n"
             "         Initial values (e.g., of flow) are used for global\n"
             "         quantities.)\n"
             "--------------------------------------------------------------\n";
    }
  };

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Solvers\n"
           "--------------------------------------------------------------\n"
           "Advection: Euler\n"
           "--------------------------------------------------------------\n";
  }
};

struct Transport {
  Transport() = delete;

  struct Parameters {
  private:
    std::string comment_sequence =
        "#"; /**< Sequence of characters marking comment for file parsing. */

  public:
    std::string rescale_velocity_option;
    double lengthscale;
    const double diff_coeff{0.};
    const double diffusion_time{std::numeric_limits<double>::infinity()};
    double advection_time;
    const double peclet{std::numeric_limits<double>::infinity()};
    double mean_velocity;
    double velocity_rescaling_factor;

    template <typename Geometry>
    Parameters(Directories const &directories, std::string const &name,
               Geometry const &geometry) {
      auto input = useful::open_read(directories.dir_parameters +
                                     "/parameters_transport_" + name + ".dat");
      useful::read_first_from_line(input, rescale_velocity_option,
                                   comment_sequence);
      useful::read_first_from_line(input, lengthscale, comment_sequence);
      if (rescale_velocity_option == "rescale_velocity_to_mean") {
        useful::read_first_from_line(input, mean_velocity, comment_sequence);
        advection_time = lengthscale / mean_velocity;
      } else if (rescale_velocity_option ==
                 "rescale_velocity_to_advection_time") {
        useful::read_first_from_line(input, advection_time, comment_sequence);
        mean_velocity = lengthscale / advection_time;
      } else if (rescale_velocity_option == "no_rescale_velocity") {
      } else
        throw std::runtime_error{"Flow velocity field rescaling option " +
                                 rescale_velocity_option + " not supported"};
      input.close();
    }

    template <typename VelocityField, typename Mesh>
    void rescale(VelocityField &velocity_field, Mesh const &mesh) {
      double current_mean =
          ptof::magnitude_of_average(velocity_field.field(), mesh);
      if (rescale_velocity_option == "rescale_velocity_to_mean" ||
          rescale_velocity_option == "rescale_velocity_to_advection_time") {
        velocity_rescaling_factor = mean_velocity / current_mean;
        velocity_field.rescale(velocity_rescaling_factor);
      } else if (rescale_velocity_option == "no_rescale_velocity") {
        velocity_rescaling_factor = 1.;
        advection_time = lengthscale / current_mean;
        mean_velocity = lengthscale / advection_time;
      } else
        throw std::runtime_error{"Flow velocity field rescaling option " +
                                 rescale_velocity_option + " not supported"};
    }

    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Transport parameters\n"
             "--------------------------------------------------------------\n"
             "- Whether and how to rescale velocity field:\n"
             "  (Note: Rescaling of the velocity field is not handled when\n"
             "         transport parameters are set.\n"
             "         Rescale method must always be called, either to\n"
             "         rescale or to compute velocity-dependent quantities\n"
             "         when not rescaling)\n"
             "\trescale_velocity_to_mean: Rescale according to imposed mean\n"
             "\t                          flow velocity\n"
             "\trescale_velocity_to_advection_time: Rescale according to\n"
             "\t                                    imposed advection time\n"
             "\tno_rescale_velocity: Do not rescale\n"
             "- Reference lengthscale\n"
             "- Mean flow velocity (pass only if rescaling with\n"
             "  rescale_velocity_to_mean)\n"
             "- Advection time (pass only if rescaling with\n"
             "  rescale_velocity_to_advection_time)\n"
             "--------------------------------------------------------------\n";
    }
  };

  template <typename Solvers, typename VelocityField, typename Geometry,
            typename Boundary, typename ReactionParameters>
  static auto
  makeTransitions(VelocityField const &velocity_field, Geometry const &geometry,
                  Boundary &boundary, Parameters const &params_transport,
                  ReactionParameters const &params_reaction,
                  typename Solvers::Parameters const &params_solvers) {
    return makeTransportTransitions_Advection<typename Solvers::Steppers>(
        velocity_field, geometry, boundary, params_transport, params_reaction,
        params_solvers);
  }

  template <typename Geometry, typename TransportParameters>
  static auto makeVelocityInterpolator(Geometry const &geometry,
                                       TransportParameters &params_transport) {
    return makeLinearVelocityInterpolator(
        geometry, ptof::get_velocity_data(geometry.mesh()), params_transport);
  }

  template <typename Geometry, typename Uninterpolated,
            typename TransportParameters>
  static auto makeVelocityInterpolator(Geometry const &geometry,
                                       Uninterpolated &&uninterpolated,
                                       TransportParameters &params_transport) {
    return makeLinearVelocityInterpolator(
        geometry, ptof::get_velocity_data(geometry.mesh()),
        std::forward<Uninterpolated>(uninterpolated), params_transport);
  }

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Transport\n"
           "--------------------------------------------------------------\n"
           "Process: Advection\n"
           "Interpolation: Linear\n"
           "--------------------------------------------------------------\n";
  }
};
}; // namespace model_advection_2d

/** \namespace ptof::model_advection_diffusion_fpt_2d
 Definitions for first-passage times under 2D advective--diffusive transport. */
namespace model_advection_diffusion_fpt_2d {
struct Model {
  Model() = delete;

  inline static const std::string name{"advection_diffusion_fpt_2d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using Geometry =
    Geometry<2, meta::ParallelOptions::Serial, Dynamics::Type::firstpassage>;
using Info = ptof::Info_Absorbed_Reinjections;
using State = State<Geometry::dim, Info, double, double, std::size_t>;
using CTRW = ctrw::CTRW<State>;
using model_advection_diffusion_2d::InitialCondition;
using model_advection_diffusion_2d::makeTransitions;
using model_advection_diffusion_2d::Output;
using model_advection_diffusion_2d::Reaction;
using model_advection_diffusion_2d::Solvers;
using model_advection_diffusion_2d::Transport;
} // namespace model_advection_diffusion_fpt_2d

/** \namespace ptof::model_advection_diffusion_surface_decay_2d
 Definitions for 2D advective--diffusive transport with surface reaction \f$ A_F
 + B_S \to B_S\f$. */
namespace model_advection_diffusion_surface_decay_2d {
struct Model {
  Model() = delete;

  inline static const std::string name{"advection_diffusion_surface_decay_2d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d::CTRW;
using model_advection_diffusion_2d::Geometry;
using model_advection_diffusion_2d::Info;
using model_advection_diffusion_2d::InitialCondition;
using model_advection_diffusion_2d::makeTransitions;
using model_advection_diffusion_2d::Output;
using model_advection_diffusion_2d::Solvers;
using model_advection_diffusion_2d::State;
using model_advection_diffusion_2d::Transport;

struct Reaction {
  Reaction() = delete;

  using BulkReaction = BulkReaction_DoNothing;
  using SurfaceReaction = SurfaceReaction_AFluidPlusASolidtoASolid;

  struct Parameters {
  private:
    std::string comment_sequence =
        "#"; /**< Sequence of characters marking comment for file parsing. */

  public:
    double damkohler;
    std::string initial_distribution;
    double surface_concentration;
    double rate_constant;
    double reaction_time;

    template <typename Geometry, typename TransportParameters>
    Parameters(Directories const &directories, std::string const &name,
               Geometry const &geometry,
               TransportParameters const &params_transport) {
      auto input = useful::open_read(directories.dir_parameters +
                                     "/parameters_reaction_" + name + ".dat");
      useful::read_first_from_line(input, damkohler, comment_sequence);
      useful::read_first_from_line(input, initial_distribution,
                                   comment_sequence);
      useful::read_first_from_line(input, surface_concentration,
                                   comment_sequence);

      rate_constant = params_transport.lengthscale * damkohler /
                      (surface_concentration * params_transport.diffusion_time);
      reaction_time = params_transport.diffusion_time / damkohler;
      input.close();
    }

    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Reaction parameters:\n"
             "--------------------------------------------------------------\n"
             "- Damkohler number\n"
             "- Solid reactant initial distribution type:\n"
             "\tuniform: Homogeneous throughout the domain\n"
             "- Initial solid reactant surface concentration\n"
             "--------------------------------------------------------------\n";
    }
  };

  template <typename Geometry, typename TransportParameters,
            typename SolverParameters>
  static auto makeBulkReaction(Geometry const &geometry,
                               Parameters const &parameters,
                               TransportParameters const &params_transport,
                               SolverParameters const &params_solvers) {
    return BulkReaction{};
  }

  template <typename Geometry, typename TransportParameters,
            typename SolverParameters>
  static auto makeSurfaceReaction(Geometry const &geometry,
                                  Parameters const &parameters,
                                  TransportParameters const &params_transport,
                                  SolverParameters const &params_solvers) {
    if (parameters.initial_distribution == "uniform")
      return SurfaceReaction{
          parameters.rate_constant, params_transport.diff_coeff,
          uniform_solid_reactant_patches(parameters.surface_concentration,
                                         {"wallFluidSolid"}, geometry.mesh())};
    throw std::runtime_error{"Initial reactant distribution " +
                             parameters.initial_distribution +
                             " not supported"};
  }

  template <typename OStream> static void info(OStream &output) {
    BulkReaction::info(output);
    output << "\n";
    SurfaceReaction::info(output);
  }
};
} // namespace model_advection_diffusion_surface_decay_2d

/** \namespace ptof::model_periodic_cartesian_advection_diffusion_2d
 Definitions for 2D advective--diffusive transport with some periodic boundaries
 aligned with the Cartesian axes. */
namespace model_periodic_cartesian_advection_diffusion_2d {
struct Model {
  Model() = delete;

  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_2d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using Geometry = Geometry_Periodic_Cartesian<2, meta::ParallelOptions::Serial>;
using model_advection_diffusion_2d::Info;
using State = State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
using CTRW = ctrw::CTRW<State>;
using model_advection_diffusion_2d::InitialCondition;
using model_advection_diffusion_2d::makeTransitions;
using model_advection_diffusion_2d::Output;
using model_advection_diffusion_2d::Reaction;
using model_advection_diffusion_2d::Solvers;
using model_advection_diffusion_2d::Transport;
} // namespace model_periodic_cartesian_advection_diffusion_2d

/** \namespace ptof::model_periodic_cartesian_advection_2d
 Definitions for 2D advective transport with some periodic boundaries aligned
 with the Cartesian axes. */
namespace model_periodic_cartesian_advection_2d {
struct Model {
  Model() = delete;

  inline static const std::string name{"periodic_cartesian_advection_2d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using model_advection_2d::InitialCondition;
using model_advection_2d::makeTransitions;
using model_advection_2d::Output;
using model_advection_2d::Reaction;
using model_advection_2d::Solvers;
using model_advection_2d::Transport;
using model_periodic_cartesian_advection_diffusion_2d::CTRW;
using model_periodic_cartesian_advection_diffusion_2d::Geometry;
using model_periodic_cartesian_advection_diffusion_2d::Info;
using model_periodic_cartesian_advection_diffusion_2d::State;
}; // namespace model_periodic_cartesian_advection_2d

/** \namespace ptof::model_periodic_cartesian_advection_diffusion_fpt_2d
 Definitions for first-passage times under 2D advective--diffusive transport
 with some periodic boundaries aligned with the Cartesian axes. */
namespace model_periodic_cartesian_advection_diffusion_fpt_2d {
struct Model {
  Model() = delete;

  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_fpt_2d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using Geometry = Geometry_Periodic_Cartesian<2, meta::ParallelOptions::Serial,
                                             Dynamics::Type::firstpassage>;
using model_advection_diffusion_fpt_2d::Info;
using State = State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
using CTRW = ctrw::CTRW<State>;
using model_periodic_cartesian_advection_diffusion_2d::InitialCondition;
using model_periodic_cartesian_advection_diffusion_2d::makeTransitions;
using model_periodic_cartesian_advection_diffusion_2d::Output;
using model_periodic_cartesian_advection_diffusion_2d::Reaction;
using model_periodic_cartesian_advection_diffusion_2d::Solvers;
using model_periodic_cartesian_advection_diffusion_2d::Transport;
} // namespace model_periodic_cartesian_advection_diffusion_fpt_2d

/** \namespace
 ptof::model_periodic_cartesian_advection_diffusion_surface_decay_2d Definitions
 for first-passage times under 2D advective--diffusive transport with surface
 reaction \f$ A_F + B_S \to B_S\f$, with some periodic boundaries aligned with
 the Cartesian axes. */
namespace model_periodic_cartesian_advection_diffusion_surface_decay_2d {
struct Model {
  Model() = delete;

  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_surface_decay_2d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using model_advection_diffusion_surface_decay_2d::Reaction;
using model_periodic_cartesian_advection_diffusion_2d::CTRW;
using model_periodic_cartesian_advection_diffusion_2d::Geometry;
using model_periodic_cartesian_advection_diffusion_2d::Info;
using model_periodic_cartesian_advection_diffusion_2d::InitialCondition;
using model_periodic_cartesian_advection_diffusion_2d::makeTransitions;
using model_periodic_cartesian_advection_diffusion_2d::Output;
using model_periodic_cartesian_advection_diffusion_2d::Solvers;
using model_periodic_cartesian_advection_diffusion_2d::State;
using model_periodic_cartesian_advection_diffusion_2d::Transport;
} // namespace model_periodic_cartesian_advection_diffusion_surface_decay_2d

/** \namespace ptof::model_advection_diffusion_3d
 Definitions for 3D advective--diffusive transport. */
namespace model_advection_diffusion_3d {
struct Model {
  Model() = delete;

  inline static const std::string name{"advection_diffusion_3d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using Geometry = Geometry<3, meta::ParallelOptions::Serial>;
using Info = ptof::Info_Absorbed;
using State = State<Geometry::dim, Info, double, double, std::size_t>;
using CTRW = ctrw::CTRW<State>;
using model_advection_diffusion_2d::InitialCondition;
using model_advection_diffusion_2d::makeTransitions;
using model_advection_diffusion_2d::Output;
using model_advection_diffusion_2d::Reaction;
using model_advection_diffusion_2d::Solvers;
using model_advection_diffusion_2d::Transport;
} // namespace model_advection_diffusion_3d

/** \namespace ptof::model_advection_3d
 Definitions for 3D advective transport. */
namespace model_advection_3d {
struct Model {
  Model() = delete;

  inline static const std::string name{"advection_3d"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------------------\n"
           << name << "\n"
           << "--------------------------------------------------------------\n";
  }
};

using model_advection_2d::Solvers;
using model_advection_2d::Transport;
using model_advection_diffusion_3d::CTRW;
using model_advection_diffusion_3d::Geometry;
using model_advection_diffusion_3d::Info;
using model_advection_diffusion_3d::InitialCondition;
using model_advection_diffusion_3d::makeTransitions;
using model_advection_diffusion_3d::Output;
using model_advection_diffusion_3d::Reaction;
using model_advection_diffusion_3d::State;
}; // namespace model_advection_3d

/** \namespace ptof::model_advection_diffusion_fpt_3d
 Definitions for first-passage times under 3D advective--diffusive transport. */
namespace model_advection_diffusion_fpt_3d {
struct Model {
  Model() = delete;

  inline static const std::string name{"advection_diffusion_fpt_3d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using Geometry =
    Geometry<3, meta::ParallelOptions::Serial, Dynamics::Type::firstpassage>;
using model_advection_diffusion_fpt_2d::Info;
using State = State<Geometry::dim, Info, double, double, std::size_t>;
using CTRW = ctrw::CTRW<State>;
using model_advection_diffusion_3d::InitialCondition;
using model_advection_diffusion_3d::makeTransitions;
using model_advection_diffusion_3d::Output;
using model_advection_diffusion_3d::Reaction;
using model_advection_diffusion_3d::Solvers;
using model_advection_diffusion_3d::Transport;
} // namespace model_advection_diffusion_fpt_3d

/** \namespace ptof::model_advection_diffusion_surface_decay_3d
 Definitions for 3D advective--diffusive transport with surface reaction \f$ A_F
 + B_S \to B_S\f$. */
namespace model_advection_diffusion_surface_decay_3d {
struct Model {
  Model() = delete;

  inline static const std::string name{"advection_diffusion_surface_decay_3d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using model_advection_diffusion_3d::CTRW;
using model_advection_diffusion_3d::Geometry;
using model_advection_diffusion_3d::Info;
using model_advection_diffusion_3d::InitialCondition;
using model_advection_diffusion_3d::makeTransitions;
using model_advection_diffusion_3d::Output;
using model_advection_diffusion_3d::Solvers;
using model_advection_diffusion_3d::State;
using model_advection_diffusion_3d::Transport;
using model_advection_diffusion_surface_decay_2d::Reaction;
} // namespace model_advection_diffusion_surface_decay_3d

/** \namespace ptof::model_periodic_cartesian_advection_diffusion_3d
 Definitions for 3D advective transport with some periodic boundaries aligned
 with the Cartesian axes. */
namespace model_periodic_cartesian_advection_diffusion_3d {
struct Model {
  Model() = delete;

  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_3d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using Geometry = Geometry_Periodic_Cartesian<3, meta::ParallelOptions::Serial>;
using model_advection_diffusion_3d::Info;
using State = State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
using CTRW = ctrw::CTRW<State>;
using model_advection_diffusion_3d::InitialCondition;
using model_advection_diffusion_3d::makeTransitions;
using model_advection_diffusion_3d::Output;
using model_advection_diffusion_3d::Reaction;
using model_advection_diffusion_3d::Solvers;
using model_advection_diffusion_3d::Transport;
} // namespace model_periodic_cartesian_advection_diffusion_3d

/** \namespace ptof::model_periodic_cartesian_advection_3d
 Definitions for 3D advective transport with some periodic boundaries aligned
 with the Cartesian axes. */
namespace model_periodic_cartesian_advection_3d {
struct Model {
  Model() = delete;

  inline static const std::string name{"periodic_cartesian_advection_3d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using model_advection_3d::Reaction;
using model_advection_3d::Solvers;
using model_advection_3d::Transport;
using model_periodic_cartesian_advection_diffusion_3d::CTRW;
using model_periodic_cartesian_advection_diffusion_3d::Geometry;
using model_periodic_cartesian_advection_diffusion_3d::Info;
using model_periodic_cartesian_advection_diffusion_3d::InitialCondition;
using model_periodic_cartesian_advection_diffusion_3d::makeTransitions;
using model_periodic_cartesian_advection_diffusion_3d::Output;
using model_periodic_cartesian_advection_diffusion_3d::State;
}; // namespace model_periodic_cartesian_advection_3d

/** \namespace ptof::model_periodic_cartesian_advection_diffusion_fpt_3d
 Definitions for first-passage times under 3D advective--diffusive transport
 with some periodic boundaries aligned with the Cartesian axes. */
namespace model_periodic_cartesian_advection_diffusion_fpt_3d {
struct Model {
  Model() = delete;

  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_fpt_3d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using Geometry = Geometry_Periodic_Cartesian<3, meta::ParallelOptions::Serial,
                                             Dynamics::Type::firstpassage>;
using model_advection_diffusion_fpt_3d::Info;
using State = State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
using CTRW = ctrw::CTRW<State>;
using model_periodic_cartesian_advection_diffusion_3d::InitialCondition;
using model_periodic_cartesian_advection_diffusion_3d::makeTransitions;
using model_periodic_cartesian_advection_diffusion_3d::Output;
using model_periodic_cartesian_advection_diffusion_3d::Reaction;
using model_periodic_cartesian_advection_diffusion_3d::Solvers;
using model_periodic_cartesian_advection_diffusion_3d::Transport;
} // namespace model_periodic_cartesian_advection_diffusion_fpt_3d

/** \namespace
 ptof::model_periodic_cartesian_advection_diffusion_surface_decay_3d Definitions
 for 3D advective--diffusive transport with surface reaction \f$ A_F + B_S \to
 B_S\f$, with some periodic boundaries aligned with the Cartesian axes. */
namespace model_periodic_cartesian_advection_diffusion_surface_decay_3d {
struct Model {
  Model() = delete;

  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_surface_decay_3d"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using model_advection_diffusion_surface_decay_3d::Reaction;
using model_periodic_cartesian_advection_diffusion_3d::CTRW;
using model_periodic_cartesian_advection_diffusion_3d::Geometry;
using model_periodic_cartesian_advection_diffusion_3d::Info;
using model_periodic_cartesian_advection_diffusion_3d::InitialCondition;
using model_periodic_cartesian_advection_diffusion_3d::makeTransitions;
using model_periodic_cartesian_advection_diffusion_3d::Output;
using model_periodic_cartesian_advection_diffusion_3d::Solvers;
using model_periodic_cartesian_advection_diffusion_3d::State;
using model_periodic_cartesian_advection_diffusion_3d::Transport;
} // namespace model_periodic_cartesian_advection_diffusion_surface_decay_3d

/** \namespace ptof::model_bcc_cartesian_advection_diffusion
 Definitions for advective--diffusive transport in a body centered cubic
 beadpack, based on the primitive unit cell. */
namespace model_bcc_cartesian_advection_diffusion {
struct Model {
  Model() = delete;

  inline static const std::string name{"bcc_cartesian_advection_diffusion"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using Geometry = Geometry_Bcc<meta::ParallelOptions::Serial>;
using model_periodic_cartesian_advection_diffusion_3d::CTRW;
using model_periodic_cartesian_advection_diffusion_3d::Info;
using model_periodic_cartesian_advection_diffusion_3d::InitialCondition;
using model_periodic_cartesian_advection_diffusion_3d::makeTransitions;
using model_periodic_cartesian_advection_diffusion_3d::Output;
using model_periodic_cartesian_advection_diffusion_3d::Reaction;
using model_periodic_cartesian_advection_diffusion_3d::Solvers;
using model_periodic_cartesian_advection_diffusion_3d::State;

struct Transport {
  Transport() = delete;

  struct Parameters {
  private:
    std::string comment_sequence =
        "#"; /**< Sequence of characters marking comment for file parsing. */

  public:
    std::string peclet_option;
    std::string lengthscale_option;
    std::string primitive_cell_location_option;
    double lengthscale;
    double peclet;
    double diff_coeff;
    double diffusion_time;
    double mean_velocity;
    double velocity_rescaling_factor;

    double cell_side;
    std::vector<std::pair<double, double>> primitive_cell_boundaries;
    double advection_time;

    template <typename Geometry>
    Parameters(Directories const &directories, std::string const &name,
               Geometry const &geometry) {
      auto input = useful::open_read(directories.dir_parameters +
                                     "/parameters_transport_" + name + ".dat");
      useful::read_first_from_line(input, peclet_option, comment_sequence);
      useful::read_first_from_line(input, lengthscale_option, comment_sequence);
      double radius = geometry.radius;
      cell_side = 4. / std::sqrt(3.) * radius;
      if (lengthscale_option == "radius")
        lengthscale = radius;
      else if (lengthscale_option == "diameter")
        lengthscale = 2. * radius;
      else if (lengthscale_option == "cell_side")
        lengthscale = cell_side;
      else if (lengthscale_option == "custom")
        useful::read_first_from_line(input, lengthscale, comment_sequence);
      else
        throw std::runtime_error{"Lengthscale definition option " +
                                 lengthscale_option + " not supported"};
      if (peclet_option == "rescale_velocity_to_peclet") {
        useful::read_first_from_line(input, peclet, comment_sequence);
        useful::read_first_from_line(input, diff_coeff, comment_sequence);
        diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
        advection_time = 2. * diffusion_time / peclet;
        mean_velocity = lengthscale / advection_time;
      } else if (peclet_option == "rescale_velocity_to_mean") {
        useful::read_first_from_line(input, mean_velocity, comment_sequence);
        useful::read_first_from_line(input, diff_coeff, comment_sequence);
        diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
        peclet = lengthscale * mean_velocity / diff_coeff;
        advection_time = lengthscale / mean_velocity;
      } else if (peclet_option == "rescale_velocity_to_advection_time") {
        useful::read_first_from_line(input, advection_time, comment_sequence);
        useful::read_first_from_line(input, diff_coeff, comment_sequence);
        diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
        peclet = 2. * diffusion_time / advection_time;
        mean_velocity = lengthscale / advection_time;
      } else if (peclet_option == "compute_from_diff_coeff") {
        useful::read_first_from_line(input, diff_coeff, comment_sequence);
        diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      } else if (peclet_option == "compute_from_diff_time") {
        useful::read_first_from_line(input, diffusion_time, comment_sequence);
        diff_coeff = lengthscale * lengthscale / (2. * diffusion_time);
      } else if (peclet_option == "set_diff_coeff") {
        useful::read_first_from_line(input, peclet, comment_sequence);
      } else
        throw std::runtime_error{"Peclet number setting option " +
                                 peclet_option + " not supported"};
      input.close();
    }

    template <typename VelocityField, typename Mesh>
    void rescale(VelocityField &velocity_field, Mesh const &mesh) {
      double current_mean =
          ptof::magnitude_of_average(velocity_field.field(), mesh);
      if (peclet_option == "rescale_velocity_to_peclet" ||
          peclet_option == "rescale_velocity_to_mean" ||
          peclet_option == "rescale_velocity_to_advection_time") {
        velocity_rescaling_factor = mean_velocity / current_mean;
        velocity_field.rescale(velocity_rescaling_factor);
      } else if (peclet_option == "compute_from_diff_coeff" ||
                 peclet_option == "compute_from_diff_time") {
        velocity_rescaling_factor = 1.;
        mean_velocity = current_mean;
        advection_time = lengthscale / mean_velocity;
        peclet = 2. * diffusion_time / advection_time;
      } else if (peclet_option == "set_diff_coeff") {
        velocity_rescaling_factor = 1.;
        mean_velocity = current_mean;
        diff_coeff = lengthscale * mean_velocity / peclet;
        diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
        advection_time = lengthscale / mean_velocity;
      } else
        throw std::runtime_error{"Peclet number setting option " +
                                 peclet_option + " not supported"};
    }

    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Transport parameters\n"
             "--------------------------------------------------------------\n"
             "- How to set the Peclet number:\n"
             "  (Note: Rescaling of the velocity field is not handled when\n"
             "         transport parameters are set.\n"
             "         Rescale method must always be called, either to\n"
             "         rescale or to compute velocity-dependent quantities\n"
             "         when not rescaling)\n"
             "\tcompute_from_diff_coeff: Compute from given diffusion\n"
             "                           coefficient (do not rescale velocity\n"
             "                           field)\n"
             "\tset_diff_coeff: Compute diffusion coefficient to set given\n"
             "\t                Peclet number (do not rescale velocity field)\n"
             "\tcompute_from_diff_time: Compute from given diffusion time\n"
             "\t                        (do not rescale velocity field)\n"
             "\trescale_velocity_to_peclet: Rescale velocity field\n"
             "\t                            according to imposed Peclet\n"
             "\t                            number\n"
             "\trescale_velocity_to_mean: Rescale according to imposed mean\n"
             "\t                          flow velocity\n"
             "\trescale_velocity_to_advection_time: Rescale according to\n"
             "\t                                    imposed advection time\n"
             "- Reference length scale definition:\n"
             "\tradius: bead radius\n"
             "\tdiameter: bead diameter\n"
             "\tcell_side: primitive cubic cell side\n"
             "\tcustom: custom value\n"
             "- Reference length scale value (pass only if reference\n"
             "  lengthscale is custom)\n"
             "- Peclet number (do not pass if Peclet option is\n"
             "  compute_from_diff_coeff or compute_from_diff_time)\n"
             "- Mean flow velocity (pass only if Peclet option is\n"
             "  rescale_velocity_to_mean)\n"
             "- Advection time (pass only if Peclet option is\n"
             "  rescale_velocity_to_advection_time)\n"
             "- Diffusion time (pass only if Peclet option is\n"
             "  compute_from_diff_time)\n"
             "- Diffusion coefficient (do not pass if Peclet option is\n"
             "  set_diff_coeff or compute_from_diff_time)\n"
             "--------------------------------------------------------------\n";
    }
  };

  template <typename Solvers, typename VelocityField, typename Geometry,
            typename Boundary, typename ReactionParameters>
  static auto
  makeTransitions(VelocityField const &velocity_field, Geometry const &geometry,
                  Boundary &boundary, Parameters const &params_transport,
                  ReactionParameters const &params_reaction,
                  typename Solvers::Parameters const &params_solvers) {
    return makeTransportTransitions<typename Solvers::Steppers>(
        velocity_field, geometry, boundary, params_transport, params_reaction,
        params_solvers);
  }

  template <typename Geometry, typename TransportParameters>
  static auto makeVelocityInterpolator(Geometry const &geometry,
                                       TransportParameters &params_transport) {
    return makeLinearVelocityInterpolator(
        geometry, ptof::get_velocity_data(geometry.mesh()), params_transport);
  }

  template <typename Geometry, typename Uninterpolated,
            typename TransportParameters>
  static auto makeVelocityInterpolator(Geometry const &geometry,
                                       Uninterpolated &&uninterpolated,
                                       TransportParameters &params_transport) {
    return makeLinearVelocityInterpolator(
        geometry, ptof::get_velocity_data(geometry.mesh()),
        std::forward<Uninterpolated>(uninterpolated), params_transport);
  }

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Transport\n"
           "--------------------------------------------------------------\n"
           "Process: Advection-diffusion\n"
           "Interpolation: Linear\n"
           "--------------------------------------------------------------\n";
  }
};
} // namespace model_bcc_cartesian_advection_diffusion

/** \namespace ptof::model_bcc_cartesian_advection_diffusion_fpt
 Definitions for first-passage times under advective--diffusive transport in a
 body centered cubic beadpack, based on the primitive unit cell. */
namespace model_bcc_cartesian_advection_diffusion_fpt {
struct Model {
  Model() = delete;

  inline static const std::string name{"bcc_cartesian_advection_diffusion_fpt"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry =
    Geometry_Bcc<meta::ParallelOptions::Serial, Periodicity::Type::cartesian,
                 Dynamics::Type::firstpassage>;
using model_bcc_cartesian_advection_diffusion::InitialCondition;
using model_bcc_cartesian_advection_diffusion::makeTransitions;
using model_bcc_cartesian_advection_diffusion::Output;
using model_bcc_cartesian_advection_diffusion::Reaction;
using model_bcc_cartesian_advection_diffusion::Solvers;
using model_bcc_cartesian_advection_diffusion::Transport;
using model_periodic_cartesian_advection_diffusion_fpt_3d::CTRW;
using model_periodic_cartesian_advection_diffusion_fpt_3d::Info;
using model_periodic_cartesian_advection_diffusion_fpt_3d::State;
} // namespace model_bcc_cartesian_advection_diffusion_fpt

/** \namespace ptof::model_bcc_cartesian_advection
 Definitions for advective transport in a body centered cubic beadpack, based on
 the primitive unit cell. */
namespace model_bcc_cartesian_advection {
struct Model {
  Model() = delete;

  inline static const std::string name{"bcc_cartesian_advection"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_bcc_cartesian_advection_diffusion::CTRW;
using model_bcc_cartesian_advection_diffusion::Geometry;
using model_bcc_cartesian_advection_diffusion::Info;
using model_bcc_cartesian_advection_diffusion::InitialCondition;
using model_bcc_cartesian_advection_diffusion::makeTransitions;
using model_bcc_cartesian_advection_diffusion::Output;
using model_bcc_cartesian_advection_diffusion::State;
using model_periodic_cartesian_advection_3d::Reaction;
using model_periodic_cartesian_advection_3d::Solvers;

struct Transport {
  Transport() = delete;

  struct Parameters {
  private:
    std::string comment_sequence =
        "#"; /**< Sequence of characters marking comment for file parsing. */

  public:
    std::string rescale_velocity_option;
    std::string lengthscale_option;
    std::string primitive_cell_location_option;
    double lengthscale;
    const double diff_coeff{0.};
    const double diffusion_time{std::numeric_limits<double>::infinity()};
    const double peclet{std::numeric_limits<double>::infinity()};
    double mean_velocity;
    double velocity_rescaling_factor;

    double cell_side;
    std::vector<std::pair<double, double>> primitive_cell_boundaries;
    double advection_time;

    template <typename Geometry>
    Parameters(Directories const &directories, std::string const &name,
               Geometry const &geometry) {
      auto input = useful::open_read(directories.dir_parameters +
                                     "/parameters_transport_" + name + ".dat");
      useful::read_first_from_line(input, rescale_velocity_option,
                                   comment_sequence);
      useful::read_first_from_line(input, lengthscale_option, comment_sequence);
      double radius = geometry.radius;
      cell_side = 4. / std::sqrt(3.) * radius;
      if (lengthscale_option == "radius")
        lengthscale = radius;
      else if (lengthscale_option == "diameter")
        lengthscale = 2. * radius;
      else if (lengthscale_option == "cell_side")
        lengthscale = cell_side;
      else if (lengthscale_option == "custom")
        useful::read_first_from_line(input, lengthscale, comment_sequence);
      else
        throw std::runtime_error{"Lengthscale definition " +
                                 lengthscale_option + " not supported"};
      if (rescale_velocity_option == "rescale_velocity_to_mean") {
        useful::read_first_from_line(input, mean_velocity, comment_sequence);
        advection_time = lengthscale / mean_velocity;
      } else if (rescale_velocity_option ==
                 "rescale_velocity_to_advection_time") {
        useful::read_first_from_line(input, advection_time, comment_sequence);
        mean_velocity = lengthscale / advection_time;
      } else if (rescale_velocity_option == "no_rescale_velocity") {
      } else
        throw std::runtime_error{"Flow velocity field rescaling option " +
                                 rescale_velocity_option + " not supported"};
      input.close();
    }

    template <typename VelocityField, typename Mesh>
    void rescale(VelocityField &velocity_field, Mesh const &mesh) {
      double current_mean =
          ptof::magnitude_of_average(velocity_field.field(), mesh);
      if (rescale_velocity_option == "rescale_velocity_to_mean" ||
          rescale_velocity_option == "rescale_velocity_to_advection_time") {
        velocity_rescaling_factor = mean_velocity / current_mean;
        velocity_field.rescale(velocity_rescaling_factor);
      } else if (rescale_velocity_option == "no_rescale_velocity") {
        velocity_rescaling_factor = 1.;
        advection_time = lengthscale / current_mean;
        mean_velocity = lengthscale / advection_time;
      } else
        throw std::runtime_error{"Flow velocity field rescaling option " +
                                 rescale_velocity_option + " not supported"};
    }

    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Transport parameters\n"
             "--------------------------------------------------------------\n"
             "- Whether and how to rescale velocity field:\n"
             "  (Note: Rescaling of the velocity field is not handled when\n"
             "         transport parameters are set.\n"
             "         Rescale method must always be called, both to\n"
             "         rescale and to compute\n velocity-dependent quantities\n"
             "         when not rescaling)\n"
             "\trescale_velocity_to_mean: Rescale according to imposed mean\n"
             "\t                          flow velocity\n"
             "\trescale_velocity_to_advection_time: Rescale according to\n"
             "                                      imposed advection time\n"
             "\tno_rescale_velocity: Do not rescale\n"
             "- Reference length scale definition:\n"
             "\tradius: bead radius\n"
             "\tdiameter: bead diameter\n"
             "\tcell_side: primitive cubic cell side\n"
             "\tcustom: custom value\n"
             "- Reference length scale value (pass only if reference\n"
             "  lengthscale is custom)\n"
             "- Peclet number (do not pass if Peclet option is\n"
             "  compute_from_diff_coeff or compute_from_diff_time)\n"
             "- Mean flow velocity (pass only if rescaling with\n"
             "  rescale_velocity_to_mean)\n"
             "- Advection time (pass only if rescaling with\n"
             "  rescale_velocity_to_advection_time)\n"
             "--------------------------------------------------------------\n";
    }
  };

  template <typename Solvers, typename VelocityField, typename Geometry,
            typename Boundary, typename ReactionParameters>
  static auto
  makeTransitions(VelocityField const &velocity_field, Geometry const &geometry,
                  Boundary &boundary, Parameters const &params_transport,
                  ReactionParameters const &params_reaction,
                  typename Solvers::Parameters const &params_solvers) {
    return makeTransportTransitions_Advection<typename Solvers::Steppers>(
        velocity_field, geometry, boundary, params_transport, params_reaction,
        params_solvers);
  }

  template <typename Geometry, typename TransportParameters>
  static auto makeVelocityInterpolator(Geometry const &geometry,
                                       TransportParameters &params_transport) {
    return makeLinearVelocityInterpolator(
        geometry, ptof::get_velocity_data(geometry.mesh()), params_transport);
  }

  template <typename Geometry, typename Uninterpolated,
            typename TransportParameters>
  static auto makeVelocityInterpolator(Geometry const &geometry,
                                       Uninterpolated &&uninterpolated,
                                       TransportParameters &params_transport) {
    return makeLinearVelocityInterpolator(
        geometry, ptof::get_velocity_data(geometry.mesh()),
        std::forward<Uninterpolated>(uninterpolated), params_transport);
  }

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Transport\n"
           "--------------------------------------------------------------\n"
           "Process: Advection\n"
           "Interpolation: Linear\n"
           "--------------------------------------------------------------\n";
  }
};
} // namespace model_bcc_cartesian_advection

/** \namespace ptof::model_bcc_cartesian_advection_diffusion_surface_decay
 Definitions for advective--diffusive transport with surface reaction \f$ A_F +
 B_S \to B_S\f$ in a body centered cubic beadpack, based on the primitive unit
 cell. */
namespace model_bcc_cartesian_advection_diffusion_surface_decay {
struct Model {
  Model() = delete;

  inline static const std::string name{
      "bcc_cartesian_advection_diffusion_surface_decay"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using model_bcc_cartesian_advection_diffusion::CTRW;
using model_bcc_cartesian_advection_diffusion::Geometry;
using model_bcc_cartesian_advection_diffusion::Info;
using model_bcc_cartesian_advection_diffusion::InitialCondition;
using model_bcc_cartesian_advection_diffusion::makeTransitions;
using model_bcc_cartesian_advection_diffusion::Output;
using model_bcc_cartesian_advection_diffusion::Solvers;
using model_bcc_cartesian_advection_diffusion::State;
using model_bcc_cartesian_advection_diffusion::Transport;
using model_periodic_cartesian_advection_diffusion_surface_decay_3d::Reaction;
} // namespace model_bcc_cartesian_advection_diffusion_surface_decay

/** \namespace ptof::model_bcc_symmetryplanes_advection_diffusion
 Definitions for advective--diffusive transport in a body centered cubic
 beadpack, based on the minimal periodic unit cell. */
namespace model_bcc_symmetryplanes_advection_diffusion {
struct Model {
  Model() = delete;

  inline static const std::string name{
      "bcc_symmetryplanes_advection_diffusion"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using Geometry = Geometry_Bcc<meta::ParallelOptions::Serial,
                              Periodicity::Type::symmetryplanes>;
using model_bcc_cartesian_advection_diffusion::CTRW;
using model_bcc_cartesian_advection_diffusion::Info;
using model_bcc_cartesian_advection_diffusion::InitialCondition;
using model_bcc_cartesian_advection_diffusion::makeTransitions;
using model_bcc_cartesian_advection_diffusion::Output;
using model_bcc_cartesian_advection_diffusion::Reaction;
using model_bcc_cartesian_advection_diffusion::Solvers;
using model_bcc_cartesian_advection_diffusion::State;
using model_bcc_cartesian_advection_diffusion::Transport;
} // namespace model_bcc_symmetryplanes_advection_diffusion

/** \namespace ptof::model_bcc_symmetryplanes_advection_diffusion_fpt
 Definitions for first-passage times under advective--diffusive transport in a
 body centered cubic beadpack, based on the minimal periodic unit cell. */
namespace model_bcc_symmetryplanes_advection_diffusion_fpt {
struct Model {
  Model() = delete;

  inline static const std::string name{"bcc_cartesian_advection_diffusion_fpt"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using Geometry = Geometry_Bcc<meta::ParallelOptions::Serial,
                              Periodicity::Type::symmetryplanes,
                              Dynamics::Type::firstpassage>;
using model_bcc_cartesian_advection_diffusion_fpt::CTRW;
using model_bcc_cartesian_advection_diffusion_fpt::Info;
using model_bcc_cartesian_advection_diffusion_fpt::State;
using model_bcc_symmetryplanes_advection_diffusion::InitialCondition;
using model_bcc_symmetryplanes_advection_diffusion::makeTransitions;
using model_bcc_symmetryplanes_advection_diffusion::Output;
using model_bcc_symmetryplanes_advection_diffusion::Reaction;
using model_bcc_symmetryplanes_advection_diffusion::Solvers;
using model_bcc_symmetryplanes_advection_diffusion::Transport;
} // namespace model_bcc_symmetryplanes_advection_diffusion_fpt

/** \namespace ptof::model_bcc_symmetryplanes_advection
 Definitions for advective transport in a body centered cubic beadpack, based on
 the minimal periodic unit cell. */
namespace model_bcc_symmetryplanes_advection {
struct Model {
  Model() = delete;

  inline static const std::string name{"bcc_symmetryplanes_advection"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using model_bcc_cartesian_advection::CTRW;
using model_bcc_cartesian_advection::Info;
using model_bcc_cartesian_advection::InitialCondition;
using model_bcc_cartesian_advection::makeTransitions;
using model_bcc_cartesian_advection::Output;
using model_bcc_cartesian_advection::Reaction;
using model_bcc_cartesian_advection::Solvers;
using model_bcc_cartesian_advection::State;
using model_bcc_cartesian_advection::Transport;
using model_bcc_symmetryplanes_advection_diffusion::Geometry;
} // namespace model_bcc_symmetryplanes_advection

/** \namespace ptof::model_bcc_symmetryplanes_advection_diffusion_surface_decay
 Definitions for advective--diffusive transport with surface reaction \f$ A_F +
 B_S \to B_S\f$ in a body centered cubic beadpack, based on the minimal periodic
 unit cell. */
namespace model_bcc_symmetryplanes_advection_diffusion_surface_decay {
struct Model {
  Model() = delete;

  inline static const std::string name{
      "bcc_symmetryplanes_advection_diffusion_surface_decay"};

  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Model\n"
           "--------------------------------------------------------------\n"
        << name << "\n"
        << "--------------------------------------------------------------\n";
  }
};

using model_bcc_cartesian_advection_diffusion_surface_decay::Reaction;
using model_bcc_symmetryplanes_advection_diffusion::CTRW;
using model_bcc_symmetryplanes_advection_diffusion::Geometry;
using model_bcc_symmetryplanes_advection_diffusion::Info;
using model_bcc_symmetryplanes_advection_diffusion::InitialCondition;
using model_bcc_symmetryplanes_advection_diffusion::makeTransitions;
using model_bcc_symmetryplanes_advection_diffusion::Output;
using model_bcc_symmetryplanes_advection_diffusion::Solvers;
using model_bcc_symmetryplanes_advection_diffusion::State;
using model_bcc_symmetryplanes_advection_diffusion::Transport;
} // namespace model_bcc_symmetryplanes_advection_diffusion_surface_decay
} // namespace ptof

#endif /* PTOF_MODELS_H */
