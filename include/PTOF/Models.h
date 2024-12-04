/**
   \file PTOF/Models.h
   \author Tomás Aquino
   \date 22/02/2022
   \brief Type definitions to implement models.
*/

#ifndef PTOF_MODELS_H
#define PTOF_MODELS_H

#include "CTRW/CTRW.h"
#include "CTRW/Meta.h"
#include "General/IO.h"
#include "General/Meta.h"
#include "General/Operations.h"
#include "General/Useful.h"
#include "PTOF/Advection.h"
#include "PTOF/Geometry.h"
#include "PTOF/Info.h"
#include "PTOF/InitialCondition_Cases.h"
#include "PTOF/Locator.h"
#include "PTOF/Output_Cases.h"
#include "PTOF/ParticleMaker.h"
#include "PTOF/Reaction.h"
#include "PTOF/State.h"
#include "PTOF/Steppers.h"
#include "PTOF/Transitions.h"
#include "PTOF/Useful.h"
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <meshSearch.H>
#include <stdexcept>
#include <string>
#include <type_traits>

/** \namespace ptof Objects and methods for ParticleTrackingOF. */
namespace ptof {
/**
   \namespace ptof::model_advection_diffusion_2d Definitions for 2D
   advective--diffusive transport.
*/
namespace model_advection_diffusion_2d {
template <typename ParallelOption> struct Definitions {
  Definitions() = delete;

  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_diffusion_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport in 2D\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Geometry = Geometry_Generic<2, ParallelOption>;
  using Info = Info_Absorbed;
  using State = State_Generic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;

  struct Solvers {
    Solvers() = delete;

    using Steppers = Steppers_Advection_Euler_Diffusion_Euler;

    struct Parameters {
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
      Parameters(Directories const &directories,
                 std::string const &parameter_set_name,
                 Geometry const &geometry,
                 TransportParameters const &params_transport,
                 ReactionParameters const &params_reaction) {
        std::string filename = directories.dir_parameters +
                               "/parameters_solvers_" + parameter_set_name +
                               ".dat";
        auto input = io::open_read(filename);
        std::string in_file = std::string{"In file "} + filename + " : ";

        auto split_line = io::split_line(input);
        std::size_t param_index = 0;
        io::read(split_line, param_index,
                 in_file + "Could not parse number of particles", nr_particles);

        split_line = io::split_line(input);
        param_index = 0;
        io::read(split_line, param_index,
                 in_file +
                     "Could not parse local and global advective, diffusive, "
                     "and reactive time step accuracy",
                 local_time_step_adv, local_time_step_diff,
                 local_time_step_react, global_time_step_adv,
                 global_time_step_diff, global_time_step_react);
      }

      /**
         \brief Output generic information about object.
         \param output Output stream.
      */
      template <typename OStream> static void info(OStream &output) {
        output
            << "--------------------------------------------------------------"
               "\n"
               "Solver parameters\n"
               "--------------------------------------------------------------"
               "\n"
               "- Number of Lagrangian particles in each injection step\n"
               "- Local and global time step accuracy:\n"
               "  (Note:\n"
               "    - Minimum between processes and maximum between local and\n"
               "      global is used.\n"
               "    - Initial values (e.g., of flow) are used for global\n"
               "      quantities)\n"
               "  - Pass on same line:\n"
               "    - Time step accuracy with respect to local advection time\n"
               "    - Time step accuracy with respect to local diffusion time\n"
               "    - Time step accuracy with respect to local reaction time\n"
               "    - Time step accuracy with respect to global advection\n"
               "      time\n"
               "    - Time step accuracy with respect to global diffusion\n"
               "      time\n"
               "    - Time step accuracy with respect to global reaction time\n"
               "--------------------------------------------------------------"
               "\n";
      }
    };

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Solvers\n"
             "--------------------------------------------------------------\n"
             "Advection: Euler\n"
             "Diffusion: Stochastic Euler\n"
             "--------------------------------------------------------------\n";
    }
  };

  struct Transport {
    Transport() = delete;

    struct Parameters {
    public:
      std::string peclet_option;
      double lengthscale;
      double diff_coeff;
      double diffusion_time;
      double advection_time;
      double peclet;
      double mean_velocity;
      double velocity_rescaling_factor;

      template <typename Geometry, typename VelocityField>
      Parameters(Directories const &directories,
                 std::string const &parameter_set_name,
                 Geometry const &geometry, VelocityField &velocity_field) {
        std::string filename = directories.dir_parameters +
                               "/parameters_transport_" + parameter_set_name +
                               ".dat";
        auto input = io::open_read(filename);
        std::string in_file = std::string{"In file "} + filename + " : ";

        auto split_line = io::split_line(input);
        std::size_t param_index = 0;
        io::read(split_line, param_index,
                 in_file + "Could not parse reference lengthscale",
                 lengthscale);

        split_line = io::split_line(input);
        param_index = 0;
        io::read(split_line, param_index,
                 in_file + "Could not parse Peclet setting option",
                 peclet_option);
        std::string for_peclet_option =
            std::string{"Peclet option "} + peclet_option + " : ";
        if (peclet_option == "rescale_velocity_to_peclet") {
          io::read(
              split_line, param_index,
              in_file + for_peclet_option +
                  "Could not parse Peclet number and diffusion coefficient",
              peclet, diff_coeff);
          diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
          advection_time = 2. * diffusion_time / peclet;
          mean_velocity = lengthscale / advection_time;
        } else if (peclet_option == "rescale_velocity_to_mean") {
          io::read(
              split_line, param_index,
              in_file + for_peclet_option +
                  "Could not parse mean velocity and diffusion coefficient",
              mean_velocity, diff_coeff);
          diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
          peclet = lengthscale * mean_velocity / diff_coeff;
          advection_time = lengthscale / mean_velocity;
        } else if (peclet_option == "rescale_velocity_to_advection_time") {
          io::read(
              split_line, param_index,
              in_file + for_peclet_option +
                  "Could not parse advection time and diffusion coefficient",
              advection_time, diff_coeff);
          diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
          peclet = 2. * diffusion_time / advection_time;
          mean_velocity = lengthscale / advection_time;
        } else if (peclet_option == "compute_from_diff_coeff") {
          io::read(split_line, param_index,
                   in_file + for_peclet_option +
                       "Could not parse diffusion coefficient",
                   diff_coeff);
          diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
        } else if (peclet_option == "compute_from_diff_time") {
          io::read(split_line, param_index,
                   in_file + for_peclet_option +
                       "Could not parse diffusion time",
                   diffusion_time);
          diff_coeff = lengthscale * lengthscale / (2. * diffusion_time);
        } else if (peclet_option == "set_diff_coeff") {
          io::read(split_line, param_index,
                   in_file + for_peclet_option +
                       "Could not parse Peclet number",
                   peclet);
        } else {
          throw std::runtime_error{in_file + for_peclet_option +
                                   "Not supported"};
        }

        rescale(velocity_field, geometry.mesh());
      }

      template <typename VelocityField, typename Mesh>
      void rescale(VelocityField &velocity_field, Mesh const &mesh) {
        double current_mean =
            magnitude_of_average(velocity_field.field(), mesh);
        std::cout << "Mean velocity = " << current_mean << std::endl;
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

      /**
         \brief Output generic information about object.
         \param output Output stream.
      */
      template <typename OStream> static void info(OStream &output) {
        output
            << "--------------------------------------------------------------"
               "\n"
               "Transport parameters\n"
               "--------------------------------------------------------------"
               "\n"
               "- Reference lengthscale\n"
               "- How to set the Peclet number:\n"
               "  - compute_from_diff_coeff\n"
               "    - Compute from given diffusion coefficient (do not\n"
               "      rescale velocity field)\n"
               "    - Pass on same line:\n"
               "      - Diffusion coefficient\n"
               "  - set_diff_coeff\n"
               "    - Set diffusion coefficient to impose given Peclet number\n"
               "      (do not rescale velocity field)\n"
               "    - Pass on same line:\n"
               "      - Peclet number\n"
               "  - compute_from_diff_time\n"
               "    - Compute from given diffusion time (do not rescale\n"
               "      velocity field)\n"
               "    - Pass on same line:\n"
               "      - Diffusion coefficient\n"
               "  - rescale_velocity_to_peclet\n"
               "    - Rescale velocity field according to given peclet number\n"
               "    - Pass on same line:\n"
               "      - Peclet number\n"
               "      - Diffusion coefficient\n"
               "  - rescale_velocity_to_mean\n"
               "    - Rescale according to given mean flow velocity\n"
               "    - Pass on same line:\n"
               "      - Absolute value of mean velocity vector\n"
               "      - Diffusion coefficient\n"
               "  - rescale_velocity_to_advection_time\n"
               "    - Rescale according to given advection time\n"
               "    - Pass on same line:\n"
               "      - Advection time, based on absolute value of mean\n"
               "        velocity vector\n"
               "      - Diffusion coefficient\n"
               "--------------------------------------------------------------"
               "\n";
      }
    };

    template <typename Solvers, typename VelocityField, typename Geometry,
              typename Boundary, typename ReactionParameters>
    static auto
    makeTransitions(VelocityField const &velocity_field,
                    Geometry const &geometry, Boundary &boundary,
                    Parameters const &params_transport,
                    ReactionParameters const &params_reaction,
                    typename Solvers::Parameters const &params_solvers) {
      return makeTransportTransitions<typename Solvers::Steppers>(
          velocity_field, geometry, boundary, params_transport, params_reaction,
          params_solvers);
    }

    template <typename Geometry>
    static auto makeVelocityInterpolator(Geometry const &geometry) {
      return makeLinearVelocityInterpolator(geometry,
                                            get_velocity_data(geometry.mesh()));
    }

    template <typename Geometry, typename Uninterpolated>
    static auto makeVelocityInterpolator(Geometry const &geometry,
                                         Uninterpolated &&uninterpolated) {
      return makeLinearVelocityInterpolator(
          geometry, get_velocity_data(geometry.mesh()),
          std::forward<Uninterpolated>(uninterpolated));
    }

    template <typename Geometry>
    static auto
    makeVelocityInterpolator_WithUninterpolated(Geometry const &geometry) {
      return makeLinearVelocityInterpolator(
          geometry, get_velocity_data(geometry.mesh()),
          Foam::volVectorField{
              Foam::IOobject{"", geometry.mesh().time().timeName(),
                             geometry.mesh(), Foam::IOobject::NO_READ,
                             Foam::IOobject::NO_WRITE},
              geometry.mesh(),
              Foam::dimensionedVector{"", Foam::dimVelocity, Foam::zero{}}});
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
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
      Parameters(Directories const &dir, std::string const &parameter_set_name,
                 Geometry const &geometry,
                 TransportParameters const &params_transport) {}

      /**
         \brief Output generic information about object.
         \param output Output stream.
      */
      template <typename OStream> static void info(OStream &output) {
        output << "------------------------------------------------------------"
                  "--\n"
                  "Reaction parameters\n"
                  "------------------------------------------------------------"
                  "--\n"
                  "None\n"
                  "------------------------------------------------------------"
                  "--\n";
      }
    };

    template <typename Geometry, typename TransportParameters,
              typename SolverParameters>
    static auto makeBulkReaction(Geometry const &geometry,
                                 Parameters const &params_reaction,
                                 TransportParameters const &params_transport,
                                 SolverParameters const &params_solvers) {
      return BulkReaction{};
    }

    template <typename Geometry, typename TransportParameters,
              typename SolverParameters>
    static auto makeSurfaceReaction(Geometry const &geometry,
                                    Parameters const &params_reaction,
                                    TransportParameters const &params_transport,
                                    SolverParameters const &params_solvers) {
      return SurfaceReaction{};
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      BulkReaction::info(output);
      output << "\n";
      SurfaceReaction::info(output);
    }
  };

  struct InitialCondition {
    InitialCondition() = delete;

    using Parameters = InitialConditionParameters_Cases;

    template <typename Particle, typename Geometry, typename VelocityField,
              typename SolverParameters, typename Mask = useful::Empty>
    static auto makeInitialCondition(Geometry const &geometry,
                                     VelocityField const &velocity_field,
                                     Parameters const &params_initial_condition,
                                     SolverParameters const &params_solvers,
                                     Mask &&mask = {}, double threshold = 0.) {
      return InitialCondition_Cases{meta::Selector_t<Particle>{},
                                    geometry,
                                    velocity_field,
                                    params_solvers.nr_particles,
                                    params_initial_condition,
                                    std::forward<Mask>(mask),
                                    threshold};
    }
  };

  struct Output {
    Output() = delete;

    using Parameters = OutputParameters_Cases;

    template <typename Subject, typename VelocityField, typename Geometry,
              typename Mask = useful::Empty>
    static auto
    makeOutput(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters const &params_output, std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks = {},
               std::vector<double> thresholds = {}) {
      return Output_Cases{subject,       velocity_field, geometry, directories,
                          params_output, identifier,     masks,    thresholds};
    }

    template <typename Subject, typename VelocityField, typename Geometry,
              typename Mask>
    static auto
    makeOutput(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters const &params_output, std::string const &identifier,
               std::initializer_list<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {}) {
      return Output_Cases{subject,       velocity_field, geometry, directories,
                          params_output, identifier,     masks,    thresholds};
    }

    template <typename Subject, typename VelocityField, typename Geometry,
              typename Mask>
    static auto
    makeOutput(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters const &params_output, std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {}) {
      return Output_Cases{subject,       velocity_field, geometry, directories,
                          params_output, identifier,     masks,    thresholds};
    }

    template <typename Subject, typename VelocityField, typename Geometry,
              typename Mask>
    static auto
    makeOutput(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters const &params_output, std::string const &identifier,
               std::initializer_list<std::reference_wrapper<const Mask>> masks,
               std::vector<double> thresholds = {}) {
      return Output_Cases{subject,       velocity_field, geometry, directories,
                          params_output, identifier,     masks,    thresholds};
    }
  };
};
} // namespace model_advection_diffusion_2d

/**
   \namespace ptof::model_advection_2d Definitions for 2D advective transport.
*/
namespace model_advection_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output << "------------------------------------------------------------"
                "--\n"
                "Model\n"
                "------------------------------------------------------------"
                "--\n"
             << "Name: " << name << "\n"
             << "Description: Advective transport in 2D\n"
             << "Parallelization: "
             << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                     ? "Serial"
                     : "Parallel")
             << "\n"
             << "------------------------------------------------------------"
                "--\n";
    }
  };

  using CTRW =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::CTRW;
  using Geometry = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Geometry;
  using Info =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Reaction;
  using State =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::State;

  struct Solvers {
    Solvers() = delete;

    using Steppers = Steppers_Advection_Euler_Diffusion_Euler;

    struct Parameters {
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
      Parameters(Directories const &directories,
                 std::string const &parameter_set_name,
                 Geometry const &geometry,
                 TransportParameters const &params_transport,
                 ReactionParameters const &params_reaction) {
        std::string filename = directories.dir_parameters +
                               "/parameters_solvers_" + parameter_set_name +
                               ".dat";
        auto input = io::open_read(filename);
        std::string in_file = std::string{"In file "} + filename + " : ";

        auto split_line = io::split_line(input);
        std::size_t param_index = 0;
        io::read(split_line, param_index,
                 in_file + "Could not parse number of particles", nr_particles);

        split_line = io::split_line(input);
        param_index = 0;
        io::read(
            split_line, param_index,
            in_file +
                "Could not parse local and global advective time step accuracy",
            local_time_step_adv, global_time_step_adv);
      }

      /**
         \brief Output generic information about object.
         \param output Output stream.
      */
      template <typename OStream> static void info(OStream &output) {
        output
            << "--------------------------------------------------------------"
               "\n"
               "Solver parameters\n"
               "--------------------------------------------------------------"
               "\n"
               "- Number of Lagrangian particles in each injection time step\n"
               "- Local and global time step accuracy:\n"
               "  (Note:\n"
               "    - Maximum between local and global is used.\n"
               "    - Initial values (e.g., of flow) are used for global\n"
               "      quantities)\n"
               "  - Pass on same line:\n"
               "    - Time step accuracy with respect to local advection time\n"
               "    - Time step accuracy with respect to global advection\n"
               "      time\n"
               "--------------------------------------------------------------"
               "\n";
      }
    };

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output << "------------------------------------------------------------"
                "--\n"
                "Solvers\n"
                "------------------------------------------------------------"
                "--\n"
                "Advection: Euler\n"
                "------------------------------------------------------------"
                "--\n";
    }
  };

  struct Transport {
    Transport() = delete;

    struct Parameters {
    public:
      std::string rescale_velocity_option;
      double lengthscale;
      const double diff_coeff{0.};
      const double diffusion_time{std::numeric_limits<double>::infinity()};
      double advection_time;
      const double peclet{std::numeric_limits<double>::infinity()};
      double mean_velocity;
      double velocity_rescaling_factor;

      template <typename Geometry, typename VelocityField>
      Parameters(Directories const &directories,
                 std::string const &parameter_set_name,
                 Geometry const &geometry, VelocityField &velocity_field) {
        std::string filename = directories.dir_parameters +
                               "/parameters_transport_" + parameter_set_name +
                               ".dat";
        auto input = io::open_read(filename);
        std::string in_file = std::string{"In file "} + filename + " : ";

        auto split_line = io::split_line(input);
        std::size_t param_index = 0;
        io::read(split_line, param_index,
                 in_file + "Could not parse reference length scale",
                 lengthscale);

        split_line = io::split_line(input);
        param_index = 0;
        io::read(split_line, param_index,
                 in_file + "Could not flow velocity field rescaling option",
                 rescale_velocity_option);
        std::string for_rescale_velocity_option =
            std::string{"Velocity rescaling option "} +
            rescale_velocity_option + " : ";
        if (rescale_velocity_option == "rescale_velocity_to_mean") {
          io::read(split_line, param_index,
                   in_file + for_rescale_velocity_option +
                       "Could not parse mean velocity",
                   mean_velocity);
          advection_time = lengthscale / mean_velocity;
        } else if (rescale_velocity_option ==
                   "rescale_velocity_to_advection_time") {
          io::read(split_line, param_index,
                   in_file + for_rescale_velocity_option +
                       "Could not parse advection time",
                   advection_time);
          mean_velocity = lengthscale / advection_time;
        } else if (rescale_velocity_option == "no_rescale_velocity") {
        } else {
          throw std::runtime_error{in_file + for_rescale_velocity_option +
                                   "Not supported"};
        }

        rescale(velocity_field, geometry.mesh());
      }

      template <typename VelocityField, typename Mesh>
      void rescale(VelocityField &velocity_field, Mesh const &mesh) {
        double current_mean =
            magnitude_of_average(velocity_field.field(), mesh);
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

      /**
         \brief Output generic information about object.
         \param output Output stream.
      */
      template <typename OStream> static void info(OStream &output) {
        output
            << "--------------------------------------------------------------"
               "\n"
               "Transport parameters\n"
               "--------------------------------------------------------------"
               "\n"
               "- Reference lengthscale\n"
               "- Whether and how to rescale velocity field:\n"
               "  - rescale_velocity_to_mean\n"
               "    - Rescale according to given mean\n"
               "    - Pass on same line:\n"
               "      - Absolute value of mean velocity vector\n"
               "  - rescale_velocity_to_advection_time\n"
               "    - Rescale according to given advection time\n"
               "    - Pass on same line:\n"
               "      - Advection time, based on absolute value of mean "
               "        velocity vector\n"
               "  - no_rescale_velocity\n"
               "    - Do not rescale\n"
               "--------------------------------------------------------------"
               "\n";
      }
    };

    template <typename Solvers, typename VelocityField, typename Geometry,
              typename Boundary, typename ReactionParameters>
    static auto
    makeTransitions(VelocityField const &velocity_field,
                    Geometry const &geometry, Boundary &boundary,
                    Parameters const &params_transport,
                    ReactionParameters const &params_reaction,
                    typename Solvers::Parameters const &params_solvers) {
      return makeTransportTransitions_Advection<typename Solvers::Steppers>(
          velocity_field, geometry, boundary, params_transport, params_reaction,
          params_solvers);
    }

    template <typename Geometry>
    static auto makeVelocityInterpolator(Geometry const &geometry) {
      return makeLinearVelocityInterpolator(geometry,
                                            get_velocity_data(geometry.mesh()));
    }

    template <typename Geometry, typename Uninterpolated>
    static auto makeVelocityInterpolator(Geometry const &geometry,
                                         Uninterpolated &&uninterpolated) {
      return makeLinearVelocityInterpolator(
          geometry, get_velocity_data(geometry.mesh()),
          std::forward<Uninterpolated>(uninterpolated));
    }

    template <typename Geometry>
    static auto
    makeVelocityInterpolator_WithUninterpolated(Geometry const &geometry) {
      return makeLinearVelocityInterpolator(
          geometry, get_velocity_data(geometry.mesh()),
          Foam::volVectorField{
              Foam::IOobject{"", geometry.mesh().time().timeName(),
                             geometry.mesh(), Foam::IOobject::NO_READ,
                             Foam::IOobject::NO_WRITE},
              geometry.mesh(),
              Foam::dimensionedVector{"", Foam::dimVelocity, Foam::zero{}}});
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output << "------------------------------------------------------------"
                "--\n"
                "Transport\n"
                "------------------------------------------------------------"
                "--\n"
                "Process: Advection\n"
                "Interpolation: Linear\n"
                "------------------------------------------------------------"
                "--\n";
    }
  };
};
} // namespace model_advection_2d

/**
   \namespace ptof::model_advection_diffusion_fpt_2d Definitions for
   first-passage times under 2D advective--diffusive transport.
*/
namespace model_advection_diffusion_fpt_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_diffusion_fpt_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: First passage times under advective-diffusive\n"
             "             transport in 2D\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Geometry =
      Geometry_Generic<2, ParallelOption, Dynamics::Type::firstpassage>;
  using Info = Info_Absorbed_Reinjections;
  using State = State_Generic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Solvers;
  using Transport = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Transport;
};
} // namespace model_advection_diffusion_fpt_2d

/**
   \namespace ptof::model_advection_diffusion_surface_decay_2d Definitions for
   2D advective--diffusive transport with surface reaction \f$ A_F
   + B_S \to B_S\f$.
*/
namespace model_advection_diffusion_surface_decay_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "advection_diffusion_surface_decay_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport with surface\n"
             "             reaction A_F + B_S -> B_S\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using CTRW =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::CTRW;
  using Geometry = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Geometry;
  using Info =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Output;
  using Solvers = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Solvers;
  using State =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::State;
  using Transport = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Transport;

  struct Reaction {
    Reaction() = delete;

    using BulkReaction = BulkReaction_DoNothing;
    using SurfaceReaction = SurfaceReaction_AFluidPlusASolidtoASolid;

    struct Parameters {
    public:
      double damkohler;
      std::string initial_distribution;
      std::vector<std::string> patch_names;
      std::vector<double> surface_concentrations;
      double rate_constant;
      double reaction_time;

      template <typename Geometry, typename TransportParameters>
      Parameters(Directories const &directories,
                 std::string const &parameter_set_name,
                 Geometry const &geometry,
                 TransportParameters const &params_transport) {
        std::string filename = directories.dir_parameters +
                               "/parameters_reaction_" + parameter_set_name +
                               ".dat";
        auto input = io::open_read(filename);
        std::string in_file = std::string{"In file "} + filename + " : ";

        auto split_line = io::split_line(input);
        std::size_t param_index = 0;
        io::read(
            split_line, param_index,
            in_file +
                "Could not parse initial surface reactant distribution type",
            initial_distribution);
        std::string for_initial_distribution =
            std::string{"Initial surface reactant distribution type "} +
            initial_distribution + " : ";
        if (initial_distribution == "uniform") {
          auto rate_option =
              io::read<std::string>(split_line, param_index,
                                    in_file + for_initial_distribution +
                                        "Could not parse rate format option");
          if (rate_option == "damkohler") {
            io::read(split_line, param_index,
                     in_file + for_initial_distribution +
                         "Could not parse damkohler number",
                     damkohler);
          } else if (rate_option == "rate_constant") {
            io::read(split_line, param_index,
                     in_file + for_initial_distribution +
                         "Could not parse rate constant",
                     rate_constant);

          } else {
            throw std::runtime_error{in_file + for_initial_distribution +
                                     "Not supported"};
          }

          std::size_t nr_patches = (split_line.size() - param_index) / 2;
          if (nr_patches == 0)
            throw std::runtime_error{
                in_file + for_initial_distribution +
                "Could not parse patch names(s) and surface concentration(s) : "
                "At least one patch and surface concentration required"};
          patch_names.reserve(nr_patches);
          surface_concentrations.reserve(nr_patches);
          for (std::size_t ii = 0; ii < nr_patches; ++ii) {
            patch_names.push_back(
                io::read<std::string>(split_line, param_index,
                                      in_file + for_initial_distribution +
                                          "Could not parse patch name"));
            surface_concentrations.push_back(
                io::read<double>(split_line, param_index,
                                 in_file + for_initial_distribution +
                                     "Could not parse initial surface "
                                     "concentration for patch " +
                                     patch_names[ii]));
          }
          auto areas = patch_areas(patch_names, geometry.mesh());

          if (rate_option == "damkohler") {
            rate_constant = params_transport.lengthscale * damkohler /
                            (op::sum(op::times(surface_concentrations, areas)) /
                             op::sum(areas) * params_transport.diffusion_time);
          } else if (rate_option == "rate_constant") {
            damkohler = op::sum(op::times(surface_concentrations, areas)) /
                        op::sum(areas) * params_transport.diffusion_time *
                        rate_constant / params_transport.lengthscale;
          }
        } else {
          throw std::runtime_error{
              std::string{"Initial surface reactant distribution type "} +
              initial_distribution + " not supported"};
        }
      }

      /**
         \brief Output generic information about object.
         \param output Output stream.
      */
      template <typename OStream> static void info(OStream &output) {
        output
            << "--------------------------------------------------------------"
               "\n"
               "Reaction parameters:\n"
               "--------------------------------------------------------------"
               "\n"
               "- Solid reactant initial distribution type:\n"
               "  - uniform\n"
               "    - Homogeneous in specified boundary patches\n"
               "    - Pass on same line:\n"
               "      - Rate specification format:\n"
               "        - damkohler\n"
               "          - Pass on same line:\n"
               "            - Damkohler number\n"
               "        - rate_constant\n"
               "          - Pass on same line:\n"
               "            - Surface rate constant\n"
               "      - Pairs of boundary patch names and surface\n"
               "        concentrations\n"
               "--------------------------------------------------------------"
               "\n";
      }
    };

    template <typename Geometry, typename TransportParameters,
              typename SolverParameters>
    static auto makeBulkReaction(Geometry const &geometry,
                                 Parameters const &params_reactions,
                                 TransportParameters const &params_transport,
                                 SolverParameters const &params_solvers) {
      return BulkReaction{};
    }

    template <typename Geometry, typename TransportParameters,
              typename SolverParameters>
    static auto makeSurfaceReaction(Geometry const &geometry,
                                    Parameters const &params_reaction,
                                    TransportParameters const &params_transport,
                                    SolverParameters const &params_solvers) {
      if (params_reaction.initial_distribution == "uniform")
        return SurfaceReaction{
            params_reaction.rate_constant, params_transport.diff_coeff,
            uniform_solid_reactant(params_reaction.patch_names,
                                   params_reaction.surface_concentrations,
                                   geometry.mesh())};
      throw std::runtime_error{"Initial reactant distribution " +
                               params_reaction.initial_distribution +
                               " not supported"};
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      BulkReaction::info(output);
      output << "\n";
      SurfaceReaction::info(output);
    }
  };
};
} // namespace model_advection_diffusion_surface_decay_2d

/**
   \namespace ptof::model_periodic_cartesian_advection_diffusion_2d Definitions
   for 2D advective--diffusive transport with some periodic boundaries aligned
   with the Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport in 2D, with some\n"
             "             periodic boundaries aligned with the Cartesian\n"
             "             axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption>;
  using Info =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::Info;
  using State =
      State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Solvers;
  using Transport = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Transport;
};
} // namespace model_periodic_cartesian_advection_diffusion_2d

/**
   \namespace ptof::model_periodic_cartesian_advection_2d Definitions for 2D
   advective transport with some periodic boundaries aligned with the Cartesian
   axes.
*/
namespace model_periodic_cartesian_advection_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"periodic_cartesian_advection_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advective transport in 2D, with some periodic\n"
             "             boundaries aligned with the Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using InitialCondition = typename model_advection_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output =
      typename model_advection_2d::Definitions<ParallelOption>::Output;
  using Reaction =
      typename model_advection_2d::Definitions<ParallelOption>::Reaction;
  using Solvers =
      typename model_advection_2d::Definitions<ParallelOption>::Solvers;
  using Transport =
      typename model_advection_2d::Definitions<ParallelOption>::Transport;
  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Geometry;
  using mInfo =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Info;
  using State =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::State;
};
} // namespace model_periodic_cartesian_advection_2d

/**
   \namespace ptof::model_periodic_cartesian_advection_diffusion_fpt_2d
   Definitions for first-passage times under 2D advective--diffusive transport
   with some periodic boundaries aligned with the Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_fpt_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_fpt_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: First-passage times under advective-diffusive\n"
             "             transport in 2D, with some periodic boundaries\n"
             "             aligned with the Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption,
                                               Dynamics::Type::firstpassage>;
  using Info = typename model_advection_diffusion_fpt_2d::Definitions<
      ParallelOption>::Info;
  using State =
      State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Output;
  using Reaction =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Reaction;
  using Solvers =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Solvers;
  using Transport =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_periodic_cartesian_advection_diffusion_fpt_2d

/**
   \namespace
   ptof::model_periodic_cartesian_advection_diffusion_surface_decay_2d
   Definitions for first-passage times under 2D advective--diffusive transport
   with surface reaction \f$ A_F + B_S \to B_S\f$, with some periodic boundaries
   aligned with the Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_surface_decay_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_decay_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: First-passage times under advection-diffusion in\n"
             "             2D, with surface reaction A_F + B_S -> B_S, with\n"
             "             some periodic boundaries aligned with the\n"
             "             Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Reaction =
      typename model_advection_diffusion_surface_decay_2d::Definitions<
          ParallelOption>::Reaction;
  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Geometry;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Output;
  using Solvers =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Solvers;
  using State =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::State;
  using Transport =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_periodic_cartesian_advection_diffusion_surface_decay_2d

/**
   \namespace ptof::model_advection_diffusion_3d Definitions for 3D
   advective--diffusive transport.
*/
namespace model_advection_diffusion_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_diffusion_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in 3D\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Geometry = Geometry_Generic<3, ParallelOption>;
  using Info = Info_Absorbed;
  using State = State_Generic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Solvers;
  using Transport = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Transport;
};
} // namespace model_advection_diffusion_3d

/**
   \namespace ptof::model_advection_3d Definitions for 3D advective transport.
*/
namespace model_advection_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection in 3D\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Solvers =
      typename model_advection_2d::Definitions<ParallelOption>::Solvers;
  using Transport =
      typename model_advection_2d::Definitions<ParallelOption>::Transport;
  using CTRW =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::CTRW;
  using Geometry = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Geometry;
  using Info =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Reaction;
  using State =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::State;
};
} // namespace model_advection_3d

/**
   \namespace ptof::model_advection_diffusion_fpt_3d Definitions for
   first-passage times under 3D advective--diffusive transport.
*/
namespace model_advection_diffusion_fpt_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_diffusion_fpt_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: First-passage times under advection-diffusion in\n"
             "             3D\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Geometry =
      Geometry_Generic<3, ParallelOption, Dynamics::Type::firstpassage>;
  using Info = typename model_advection_diffusion_fpt_2d::Definitions<
      ParallelOption>::Info;
  using State = State_Generic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Solvers;
  using Transport = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Transport;
};
} // namespace model_advection_diffusion_fpt_3d

/**
   \namespace ptof::model_advection_diffusion_surface_decay_3d Definitions for
   3D advective--diffusive transport with surface reaction \f$ A_F + B_S \to
   B_S\f$.
*/
namespace model_advection_diffusion_surface_decay_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "advection_diffusion_surface_decay_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Description: Advection-diffusion in 3D, with surface reaction\n"
             "             A_F + B_S -> B_S\n"
          << "Name: " << name << "\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using CTRW =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::CTRW;
  using Geometry = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Geometry;
  using Info =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Output;
  using Solvers = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Solvers;
  using State =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::State;
  using Transport = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Transport;
  using Reaction =
      typename model_advection_diffusion_surface_decay_2d::Definitions<
          ParallelOption>::Reaction;
};
} // namespace model_advection_diffusion_surface_decay_3d

/**
   \namespace ptof::model_periodic_cartesian_advection_diffusion_3d Definitions
   for 3D advective transport with some periodic boundaries aligned with the
   Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in 3D, with some periodic\n"
             "             boundaries aligned with the Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption>;
  using Info =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::Info;
  using State =
      State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Solvers;
  using Transport = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Transport;
};
} // namespace model_periodic_cartesian_advection_diffusion_3d

/**
   \namespace ptof::model_periodic_cartesian_advection_3d
   Definitions for 3D advective transport with some periodic boundaries aligned
   with the Cartesian axes.
*/
namespace model_periodic_cartesian_advection_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"periodic_cartesian_advection_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection in 3D, with some periodic boundaries\n"
             "             aligned with the Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Reaction =
      typename model_advection_3d::Definitions<ParallelOption>::Reaction;
  using Solvers =
      typename model_advection_3d::Definitions<ParallelOption>::Solvers;
  using Transport =
      typename model_advection_3d::Definitions<ParallelOption>::Transport;
  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Geometry;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Output;
  using State =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::State;
};
} // namespace model_periodic_cartesian_advection_3d

/**
   \namespace ptof::model_periodic_cartesian_advection_diffusion_fpt_3d
   Definitions for first-passage times under 3D advective--diffusive transport
   with some periodic boundaries aligned with the Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_fpt_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_fpt_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: First-passage times under advection-diffusion in\n"
             "             3D, with some periodic boundaries aligned with the\n"
             "             Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption,
                                               Dynamics::Type::firstpassage>;
  using Info = typename model_advection_diffusion_fpt_3d::Definitions<
      ParallelOption>::Info;
  using State =
      State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Output;
  using Reaction =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Reaction;
  using Solvers =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Solvers;
  using Transport =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_periodic_cartesian_advection_diffusion_fpt_3d

/**
   \namespace
   ptof::model_periodic_cartesian_advection_diffusion_surface_decay_3d
   Definitions for 3D advective--diffusive transport with surface reaction \f$
   A_F + B_S \to B_S\f$, with some periodic boundaries aligned with the
   Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_surface_decay_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_decay_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in 3D, with surface reaction\n"
             "             A_F + B_S -> B_S, with some periodic boundaries\n"
             "             aligned with the Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Reaction =
      typename model_advection_diffusion_surface_decay_3d::Definitions<
          ParallelOption>::Reaction;
  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Geometry;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Output;
  using Solvers =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Solvers;
  using State =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::State;
  using Transport =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_periodic_cartesian_advection_diffusion_surface_decay_3d

/**
   \namespace ptof::model_bcc_cartesian_advection_diffusion Definitions for
   advective--diffusive transport in a body centered cubic beadpack, based on
   the primitive unit cell.
*/
namespace model_bcc_cartesian_advection_diffusion {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"bcc_cartesian_advection_diffusion"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in a body centered cubic\n"
             "             beadpack, based on the primitive unit cell\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Geometry = Geometry_Bcc<ParallelOption>;
  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::CTRW;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Output;
  using Reaction =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Reaction;
  using Solvers =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Solvers;
  using State =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::State;

  struct Transport {
    Transport() = delete;

    struct Parameters {
    public:
      std::string lengthscale_option;
      std::string peclet_option;
      double lengthscale;
      double peclet;
      double diff_coeff;
      double diffusion_time;
      double mean_velocity;
      double velocity_rescaling_factor;

      double cell_side;
      std::vector<std::pair<double, double>> primitive_cell_boundaries;
      double advection_time;

      template <typename Geometry, typename VelocityField>
      Parameters(Directories const &directories,
                 std::string const &parameter_set_name,
                 Geometry const &geometry, VelocityField &velocity_field) {
        std::string filename = directories.dir_parameters +
                               "/parameters_transport_" + parameter_set_name +
                               ".dat";
        auto input = io::open_read(filename);
        std::string in_file = std::string{"In file "} + filename + " : ";

        auto split_line = io::split_line(input);
        std::size_t param_index = 0;
        io::read(split_line, param_index,
                 in_file +
                     "Could not parse reference lengthscale definition option",
                 lengthscale_option);
        std::string for_lengthscale_option =
            std::string{"Lengthscale definition option "} + lengthscale_option +
            " : ";

        double radius = geometry.radius;
        cell_side = 4. / std::sqrt(3.) * radius;
        if (lengthscale_option == "radius") {
          lengthscale = radius;
        } else if (lengthscale_option == "diameter")
          lengthscale = 2. * radius;
        else if (lengthscale_option == "cell_side") {
          lengthscale = cell_side;
        } else if (lengthscale_option == "custom") {
          io::read(split_line, param_index,
                   in_file + for_lengthscale_option +
                       "Could not parse reference lengthscale",
                   lengthscale);
          lengthscale = std::stod(split_line[param_index++]);
        } else {
          throw std::runtime_error{in_file + for_lengthscale_option +
                                   "Not supported"};
        }

        split_line = io::split_line(input);
        param_index = 0;
        io::read(split_line, param_index,
                 in_file + "Could not parse Peclet setting option",
                 peclet_option);
        std::string for_peclet_option =
            std::string{"Peclet option "} + peclet_option + " : ";
        if (peclet_option == "rescale_velocity_to_peclet") {
          io::read(
              split_line, param_index,
              in_file + for_peclet_option +
                  "Could not parse Peclet number and diffusion coefficient",
              peclet, diff_coeff);
          diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
          advection_time = 2. * diffusion_time / peclet;
          mean_velocity = lengthscale / advection_time;
        } else if (peclet_option == "rescale_velocity_to_mean") {
          io::read(
              split_line, param_index,
              in_file + for_peclet_option +
                  "Could not parse mean velocity and diffusion coefficient",
              mean_velocity, mean_velocity);
          diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
          peclet = lengthscale * mean_velocity / diff_coeff;
          advection_time = lengthscale / mean_velocity;
        } else if (peclet_option == "rescale_velocity_to_advection_time") {
          io::read(
              split_line, param_index,
              in_file + for_peclet_option +
                  "Could not parse advection time and diffusion coefficient",
              advection_time, diff_coeff);
          diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
          peclet = 2. * diffusion_time / advection_time;
          mean_velocity = lengthscale / advection_time;
        } else if (peclet_option == "compute_from_diff_coeff") {
          io::read(split_line, param_index,
                   in_file + for_peclet_option +
                       "Could not parse diffusion coefficient",
                   diff_coeff);
          diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
        } else if (peclet_option == "compute_from_diff_time") {
          io::read(split_line, param_index,
                   in_file + for_peclet_option +
                       "Could not parse diffusion time",
                   diffusion_time);
          diff_coeff = lengthscale * lengthscale / (2. * diffusion_time);
        } else if (peclet_option == "set_diff_coeff") {
          io::read(split_line, param_index,
                   in_file + for_peclet_option +
                       "Could not parse Peclet number",
                   peclet);
        } else {
          throw std::runtime_error{in_file + for_peclet_option +
                                   "Not supported"};
        }

        rescale(velocity_field, geometry.mesh());
      }

      template <typename VelocityField, typename Mesh>
      void rescale(VelocityField &velocity_field, Mesh const &mesh) {
        double current_mean =
            magnitude_of_average(velocity_field.field(), mesh);
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
        } else {
          throw std::runtime_error{"Peclet number setting option " +
                                   peclet_option + " not supported"};
        }
      }

      /**
         \brief Output generic information about object.
         \param output Output stream.
      */
      template <typename OStream> static void info(OStream &output) {
        output
            << "--------------------------------------------------------------"
               "\n"
               "Transport parameters\n"
               "--------------------------------------------------------------"
               "\n"
               "- How to set the reference lengthscale:\n"
               "  - radius\n"
               "    - Equal to bead radius\n"
               "  - diameter\n"
               "    - Equal to bead diameter\n"
               "  - cell_side\n"
               "    - Equal to primitive cell side\n"
               "  - custom\n"
               "    - Equal to specified value\n"
               "    - Pass on same line:\n"
               "      - Reference lengthscale value\n"
               "- How to set the Peclet number:\n"
               "  - compute_from_diff_coeff\n"
               "    - Compute from given diffusion coefficient (do not\n"
               "      rescale velocity field)\n"
               "    - Pass on same line:\n"
               "      - Diffusion coefficient\n"
               "  - set_diff_coeff\n"
               "    - Set diffusion coefficient to impose given Peclet number\n"
               "      (do not rescale velocity field)\n"
               "    - Pass on same line:\n"
               "      - Peclet number\n"
               "  - compute_from_diff_time\n"
               "    - Compute from given diffusion time (do not rescale\n"
               "      velocity field)\n"
               "    - Pass on same line:\n"
               "      - Diffusion coefficient\n"
               "  - rescale_velocity_to_peclet\n"
               "    - Rescale velocity field according to given peclet number\n"
               "    - Pass on same line:\n"
               "      - Peclet number\n"
               "      - Diffusion coefficient\n"
               "  - rescale_velocity_to_mean\n"
               "    - Rescale according to given mean flow velocity\n"
               "    - Pass on same line:\n"
               "      - Absolute value of mean velocity vector\n"
               "      - Diffusion coefficient\n"
               "  - rescale_velocity_to_advection_time\n"
               "    - Rescale according to given advection time\n"
               "    - Pass on same line:\n"
               "      - Advection time, based on absolute value of mean\n"
               "        velocity vector\n"
               "      - Diffusion coefficient\n"
               "--------------------------------------------------------------"
               "\n";
      }
    };

    template <typename Solvers, typename VelocityField, typename Geometry,
              typename Boundary, typename ReactionParameters>
    static auto
    makeTransitions(VelocityField const &velocity_field,
                    Geometry const &geometry, Boundary &boundary,
                    Parameters const &params_transport,
                    ReactionParameters const &params_reaction,
                    typename Solvers::Parameters const &params_solvers) {
      return makeTransportTransitions<typename Solvers::Steppers>(
          velocity_field, geometry, boundary, params_transport, params_reaction,
          params_solvers);
    }

    template <typename Geometry>
    static auto makeVelocityInterpolator(Geometry const &geometry) {
      return makeLinearVelocityInterpolator(geometry,
                                            get_velocity_data(geometry.mesh()));
    }

    template <typename Geometry, typename Uninterpolated>
    static auto makeVelocityInterpolator(Geometry const &geometry,
                                         Uninterpolated &&uninterpolated) {
      return makeLinearVelocityInterpolator(
          geometry, get_velocity_data(geometry.mesh()),
          std::forward<Uninterpolated>(uninterpolated));
    }

    template <typename Geometry>
    static auto
    makeVelocityInterpolator_WithUninterpolated(Geometry const &geometry) {
      return makeLinearVelocityInterpolator(
          geometry, get_velocity_data(geometry.mesh()),
          Foam::volVectorField{
              Foam::IOobject{"", geometry.mesh().time().timeName(),
                             geometry.mesh(), Foam::IOobject::NO_READ,
                             Foam::IOobject::NO_WRITE},
              geometry.mesh(),
              Foam::dimensionedVector{"", Foam::dimVelocity, Foam::zero{}}});
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
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
};
} // namespace model_bcc_cartesian_advection_diffusion

/**
   \namespace ptof::model_bcc_cartesian_advection_diffusion_fpt Definitions for
   first-passage times under advective--diffusive transport in a body centered
   cubic beadpack, based on the primitive unit cell.
*/
namespace model_bcc_cartesian_advection_diffusion_fpt {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_fpt"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Description: First-passage times under advection-diffusion in\n"
             "             body centered cubic beadpack, based on the\n"
             "             primitive unit cell\n"
          << "Name: " << name << "\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "-------------------------------------------------------------"
             "-\n";
    }
  };

  using Geometry = Geometry_Bcc<ParallelOption, Periodicity::Type::cartesian,
                                Dynamics::Type::firstpassage>;
  using InitialCondition =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Output;
  using Reaction =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Reaction;
  using Solvers = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Solvers;
  using Transport =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Transport;
  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_fpt_3d::Definitions<
          ParallelOption>::CTRW;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_fpt_3d::Definitions<
          ParallelOption>::Info;
  using State =
      typename model_periodic_cartesian_advection_diffusion_fpt_3d::Definitions<
          ParallelOption>::State;
};
} // namespace model_bcc_cartesian_advection_diffusion_fpt

/**
   \namespace ptof::model_bcc_cartesian_advection Definitions for advective
   transport in a body centered cubic beadpack, based on the primitive unit
   cell.
*/
namespace model_bcc_cartesian_advection {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"bcc_cartesian_advection"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Description: Advection in a body centered cubic beadpack, based\n"
             "             on the primitive unit cell\n"
          << "Name: " << name << "\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using CTRW = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::CTRW;
  using Geometry =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Geometry;
  using Info = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Info;
  using InitialCondition =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Output;
  using State = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::State;
  using Reaction = typename model_periodic_cartesian_advection_3d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_periodic_cartesian_advection_3d::Definitions<
      ParallelOption>::Solvers;

  struct Transport {
    Transport() = delete;

    struct Parameters {
    public:
      std::string lengthscale_option;
      std::string rescale_velocity_option;
      double lengthscale;
      const double diff_coeff{0.};
      const double diffusion_time{std::numeric_limits<double>::infinity()};
      const double peclet{std::numeric_limits<double>::infinity()};
      double mean_velocity;
      double velocity_rescaling_factor;

      double cell_side;
      std::vector<std::pair<double, double>> primitive_cell_boundaries;
      double advection_time;

      template <typename Geometry, typename VelocityField>
      Parameters(Directories const &directories,
                 std::string const &parameter_set_name,
                 Geometry const &geometry, VelocityField &velocity_field) {
        std::string filename = directories.dir_parameters +
                               "/parameters_transport_" + parameter_set_name +
                               ".dat";
        auto input = io::open_read(filename);
        std::string in_file = std::string{"In file "} + filename + " : ";

        auto split_line = io::split_line(input);
        std::size_t param_index = 0;
        io::read(split_line, param_index,
                 in_file + "Could not parse lengthscale definition option",
                 lengthscale_option);
        std::string for_lengthscale_option =
            std::string{"Lengthscale definition option "} + lengthscale_option +
            " : ";

        double radius = geometry.radius;
        cell_side = 4. / std::sqrt(3.) * radius;
        if (lengthscale_option == "radius") {
          lengthscale = radius;
        } else if (lengthscale_option == "diameter")
          lengthscale = 2. * radius;
        else if (lengthscale_option == "cell_side") {
          lengthscale = cell_side;
        } else if (lengthscale_option == "custom") {
          io::read(split_line, param_index,
                   in_file + for_lengthscale_option +
                       "Could not parse reference lengthscale",
                   lengthscale);
          lengthscale = std::stod(split_line[param_index++]);
        } else {
          throw std::runtime_error{"Lengthscale definition option " +
                                   lengthscale_option + " not supported"};
        }

        split_line = io::split_line(input);
        param_index = 0;
        io::read(split_line, param_index,
                 in_file + "Could not flow velocity field rescaling option",
                 rescale_velocity_option);
        std::string for_rescale_velocity_option =
            std::string{"Velocity rescaling option "} +
            rescale_velocity_option + " : ";
        if (rescale_velocity_option == "rescale_velocity_to_mean") {
          io::read(split_line, param_index,
                   in_file + for_rescale_velocity_option +
                       "Could not parse mean velocity",
                   mean_velocity);
          advection_time = lengthscale / mean_velocity;
        } else if (rescale_velocity_option ==
                   "rescale_velocity_to_advection_time") {
          io::read(split_line, param_index,
                   in_file + for_rescale_velocity_option +
                       "Could not parse advection time",
                   advection_time);
          mean_velocity = lengthscale / advection_time;
        } else if (rescale_velocity_option == "no_rescale_velocity") {
        } else {
          throw std::runtime_error{in_file + for_rescale_velocity_option +
                                   "Not supported"};
        }

        rescale(velocity_field, geometry.mesh());
      }

      template <typename VelocityField, typename Mesh>
      void rescale(VelocityField &velocity_field, Mesh const &mesh) {
        double current_mean =
            magnitude_of_average(velocity_field.field(), mesh);
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

      /**
         \brief Output generic information about object.
         \param output Output stream.
      */
      template <typename OStream> static void info(OStream &output) {
        output
            << "--------------------------------------------------------------"
               "\n"
               "Transport parameters\n"
               "--------------------------------------------------------------"
               "\n"
               "- How to set the reference lengthscale:\n"
               "  - radius\n"
               "    - Equal to bead radius\n"
               "  - diameter\n"
               "    - Equal to bead diameter\n"
               "  - cell_side\n"
               "    - Equal to primitive cell side\n"
               "  - custom\n"
               "    - Equal to specified value\n"
               "    - Pass on same line:\n"
               "      - Reference lengthscale value\n"
               "- Whether and how to rescale velocity field:\n"
               "  - rescale_velocity_to_mean\n"
               "    - Rescale according to given mean\n"
               "    - Pass on same line:\n"
               "      - Absolute value of mean velocity vector\n"
               "  - rescale_velocity_to_advection_time\n"
               "    - Rescale according to given advection time\n"
               "    - Pass on same line:\n"
               "      - Advection time, based on absolute value of mean "
               "        velocity vector\n"
               "  - no_rescale_velocity\n"
               "    - Do not rescale\n"
               "--------------------------------------------------------------"
               "\n";
      }
    };

    template <typename Solvers, typename VelocityField, typename Geometry,
              typename Boundary, typename ReactionParameters>
    static auto
    makeTransitions(VelocityField const &velocity_field,
                    Geometry const &geometry, Boundary &boundary,
                    Parameters const &params_transport,
                    ReactionParameters const &params_reaction,
                    typename Solvers::Parameters const &params_solvers) {
      return makeTransportTransitions_Advection<typename Solvers::Steppers>(
          velocity_field, geometry, boundary, params_transport, params_reaction,
          params_solvers);
    }

    template <typename Geometry>
    static auto makeVelocityInterpolator(Geometry const &geometry) {
      return makeLinearVelocityInterpolator(geometry,
                                            get_velocity_data(geometry.mesh()));
    }

    template <typename Geometry, typename Uninterpolated>
    static auto makeVelocityInterpolator(Geometry const &geometry,
                                         Uninterpolated &&uninterpolated) {
      return makeLinearVelocityInterpolator(
          geometry, get_velocity_data(geometry.mesh()),
          std::forward<Uninterpolated>(uninterpolated));
    }

    template <typename Geometry>
    static auto
    makeVelocityInterpolator_WithUninterpolated(Geometry const &geometry) {
      return makeLinearVelocityInterpolator(
          geometry, get_velocity_data(geometry.mesh()),
          Foam::volVectorField{
              Foam::IOobject{"", geometry.mesh().time().timeName(),
                             geometry.mesh(), Foam::IOobject::NO_READ,
                             Foam::IOobject::NO_WRITE},
              geometry.mesh(),
              Foam::dimensionedVector{"", Foam::dimVelocity, Foam::zero{}}});
    }

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
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
};
} // namespace model_bcc_cartesian_advection

/**
   \namespace ptof::model_bcc_cartesian_advection_diffusion_surface_decay
   Definitions for advective--diffusive transport with surface reaction \f$ A_F
   + B_S \to B_S\f$ in a body centered cubic beadpack, based on the primitive
   unit cell.
*/
namespace model_bcc_cartesian_advection_diffusion_surface_decay {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_surface_decay"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in a body centered cubic\n"
             "             beadpack, based on the primitive unit cell, with\n"
             "             surface reaction A_F + B_S -> B_S\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using CTRW = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::CTRW;
  using Geometry =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Geometry;
  using Info = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Info;
  using InitialCondition =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Output;
  using Solvers = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Solvers;
  using State = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::State;
  using Transport =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Transport;
  using Reaction =
      typename model_periodic_cartesian_advection_diffusion_surface_decay_3d::
          Definitions<ParallelOption>::Reaction;
};
} // namespace model_bcc_cartesian_advection_diffusion_surface_decay

/**
   \namespace ptof::model_bcc_symmetryplanes_advection_diffusion Definitions for
   advective--diffusive transport in a body centered cubic beadpack, based on
   the minimal unit cell.
*/
namespace model_bcc_symmetryplanes_advection_diffusion {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in a body centered cubic\n"
             "             beadpack, based on the minimal unit cell\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Geometry =
      Geometry_Bcc<ParallelOption, Periodicity::Type::symmetryplanes>;
  using CTRW = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::CTRW;
  using Info = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Info;
  using InitialCondition =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Output;
  using Reaction =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Reaction;
  using Solvers = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Solvers;
  using State = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::State;
  using Transport =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_bcc_symmetryplanes_advection_diffusion

/**
   \namespace ptof::model_bcc_symmetryplanes_advection_diffusion_fpt Definitions
   for first-passage times under advective--diffusive transport in a body
   centered cubic beadpack, based on the minimal unit cell.
*/
namespace model_bcc_symmetryplanes_advection_diffusion_fpt {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_fpt"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Description: First-passage times under advection-diffusion in\n"
             "             body centered cubic beadpack, based on the minimal\n"
             "             unit cell\n"
          << "Name: " << name << "\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Geometry =
      Geometry_Bcc<ParallelOption, Periodicity::Type::symmetryplanes,
                   Dynamics::Type::firstpassage>;
  using CTRW =
      typename model_bcc_cartesian_advection_diffusion_fpt::Definitions<
          ParallelOption>::CTRW;
  using Info =
      typename model_bcc_cartesian_advection_diffusion_fpt::Definitions<
          ParallelOption>::Info;
  using State =
      typename model_bcc_cartesian_advection_diffusion_fpt::Definitions<
          ParallelOption>::State;
  using InitialCondition =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Output;
  using Reaction =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Reaction;
  using Solvers =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Solvers;
  using Transport =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_bcc_symmetryplanes_advection_diffusion_fpt

/**
   \namespace ptof::model_bcc_symmetryplanes_advection Definitions for advective
   transport in a body centered cubic beadpack, based on the minimal unit cell.
*/
namespace model_bcc_symmetryplanes_advection {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"bcc_symmetryplanes_advection"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using CTRW =
      typename model_bcc_cartesian_advection::Definitions<ParallelOption>::CTRW;
  using Info =
      typename model_bcc_cartesian_advection::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::Solvers;
  using State = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::State;
  using Transport = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::Transport;
  using Geometry =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Geometry;
};
} // namespace model_bcc_symmetryplanes_advection

/**
   \namespace ptof::model_bcc_symmetryplanes_advection_diffusion_surface_decay
   Definitions for advective--diffusive transport with surface reaction \f$ A_F
   + B_S \to B_S\f$ in a body centered cubic beadpack, based on the minimal unit
   cell.
*/
namespace model_bcc_symmetryplanes_advection_diffusion_surface_decay {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion_surface_decay"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    template <typename OStream> static void info(OStream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in a body centered cubic\n"
             "             beadpack, based on the minimal unit cell, with\n"
             "             surface reaction A_F + B_S -> B_S\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
    }
  };

  using Reaction =
      typename model_bcc_cartesian_advection_diffusion_surface_decay::
          Definitions<ParallelOption>::Reaction;
  using CTRW =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Geometry;
  using Info =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Output;
  using Solvers =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Solvers;
  using State =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::State;
  using Transport =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_bcc_symmetryplanes_advection_diffusion_surface_decay
} // namespace ptof

#endif /* PTOF_MODELS_H */
