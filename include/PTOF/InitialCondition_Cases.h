/**
   \file PTOF/InitialCondition_Cases.h
   \author Tomás Aquino
   \date 24/02/2022
   \brief Create particles according to different initial conditions.
*/

#ifndef PTOF_INITIALCONDITION_CASES_H
#define PTOF_INITIALCONDITION_CASES_H

#include "CTRW/Meta.h"
#include "General/IO.h"
#include "General/Meta.h"
#include "General/Operations.h"
#include "PTOF/Boundary.h"
#include "PTOF/InitialCondition.h"
#include "PTOF/InitialConditionList.h"
#include "PTOF/ParticleMaker.h"
#include "PTOF/TimeUnits.h"
#include "PTOF/Useful.h"
#include "Stochastic/Random.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fieldTypes.H>
#include <limits>
#include <map>
#include <memory>
#include <point.H>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ptof {
/**
   \class InitialConditionParameters_Cases PTOF/InitialCondition_Cases.h
   "PTOF/InitialCondition_Cases.h"
   \brief Initial condition condition parameters to handle all initial condition
   types in \c InitialCondition_Cases.
*/
struct InitialConditionParameters_Cases {
public:
  std::string name;
  InitialConditionList::Type type;
  std::string continuity;

  struct SpecificParameters {
    virtual ~SpecificParameters() {}
  };

  struct SpecificParameters_Patches : public SpecificParameters {
    std::vector<std::string> patch_names;

    SpecificParameters_Patches(std::vector<std::string> patch_names)
        : patch_names{patch_names} {}
  };

  struct SpecificParameters_Position : public SpecificParameters {
    std::vector<double> position;

    SpecificParameters_Position(std::vector<double> position)
        : position{position} {}
  };

  struct SpecificParameters_Patches_Distance
      : public SpecificParameters_Patches {
    double distance;
    std::size_t nr_tries;

    SpecificParameters_Patches_Distance(std::vector<std::string> patch_names,
                                        double distance, std::size_t nr_tries)
        : SpecificParameters_Patches{patch_names}, distance{distance},
          nr_tries{nr_tries} {
      if (distance < 0.)
        throw std::runtime_error{"Negative distance to patches"};
    }
  };

  struct SpecificParameters_Filename : public SpecificParameters {
    std::string filename;

    SpecificParameters_Filename(std::string filename) : filename{filename} {}
  };

  struct SpecificParameters_Boundaries : public SpecificParameters {
    std::vector<std::pair<double, double>> region_boundaries;

    SpecificParameters_Boundaries(
        std::vector<std::pair<double, double>> region_boundaries)
        : region_boundaries{region_boundaries} {
      for (auto const &bounds : region_boundaries)
        if (bounds.second < bounds.first)
          throw std::runtime_error{
              "Cartesian region upper bounds smaller than lower bounds"};
    }
  };

  struct SpecificParameters_Cylinder : public SpecificParameters {
    std::vector<double> center_1;
    std::vector<double> center_2;
    double radius;
    std::size_t nr_tries;
    double inner_radius{0.};

    SpecificParameters_Cylinder(std::vector<double> center_1,
                                std::vector<double> center_2, double radius,
                                std::size_t nr_tries, double inner_radius = 0.)
        : center_1{center_1}, center_2{center_2}, radius{radius},
          nr_tries{nr_tries}, inner_radius{inner_radius} {
      if (inner_radius >= radius)
        throw std::runtime_error{"Inner radius not smaller than radius"};
    }
  };

  struct SpecificParameters_Parallelipiped : public SpecificParameters {
    std::vector<double> corner;
    std::vector<std::vector<double>> sides;
    std::size_t nr_tries;

    SpecificParameters_Parallelipiped(std::vector<double> corner,
                                      std::vector<std::vector<double>> sides,
                                      std::size_t nr_tries)
        : corner{corner}, sides{sides}, nr_tries{nr_tries} {}
  };

  struct SpecificParameters_Sphere : public SpecificParameters {
    std::vector<double> center;
    double radius;
    std::size_t nr_tries;
    double inner_radius{0.};

    SpecificParameters_Sphere(std::vector<double> center, double radius,
                              std::size_t nr_tries, double inner_radius = 0.)
        : center{center}, radius{radius}, nr_tries{nr_tries},
          inner_radius{inner_radius} {
      if (inner_radius >= radius)
        throw std::runtime_error{"Inner radius not smaller than radius"};
    }
  };

  struct InjectionParameters {
    virtual ~InjectionParameters() {}
  };

  struct InjectionParameters_PulseNoMass : public InjectionParameters {
    double time;

    InjectionParameters_PulseNoMass(double time) : time{time} {}
  };

  struct InjectionParameters_Pulse : public InjectionParameters_PulseNoMass {
    double mass;

    InjectionParameters_Pulse(double time, double mass)
        : InjectionParameters_PulseNoMass{time}, mass{mass} {}

    InjectionParameters_Pulse(
        InjectionParameters_PulseNoMass const &specific_parameters_pulse,
        double mass)
        : InjectionParameters_Pulse{specific_parameters_pulse.time, mass} {}
  };

  struct InjectionParameters_ContinuousNoMass
      : public InjectionParameters_PulseNoMass {
    double time_end;
    double time_step;

    InjectionParameters_ContinuousNoMass(double time_start, double time_end,
                                         double time_step)
        : InjectionParameters_PulseNoMass{time_start}, time_end{time_end},
          time_step{time_step} {
      if (time_end <= time_start)
        throw std::runtime_error{
            "Injection end time not larger than injection start time"};
      if (time_step <= 0)
        throw std::runtime_error{"Injection time step not positive"};
    }

    InjectionParameters_ContinuousNoMass(
        InjectionParameters_PulseNoMass const &specific_parameters_pulse,
        double time_end, double time_step)
        : InjectionParameters_ContinuousNoMass{specific_parameters_pulse.time,
                                               time_end, time_step} {}
  };

  struct InjectionParameters_Continuous
      : public InjectionParameters_ContinuousNoMass {
    double mass_per_time_step;

    InjectionParameters_Continuous(double time_start, double time_end,
                                   double time_step, double mass_per_time_step)
        : InjectionParameters_ContinuousNoMass{time_start, time_end, time_step},
          mass_per_time_step{mass_per_time_step} {}

    InjectionParameters_Continuous(
        InjectionParameters_ContinuousNoMass const &specific_parameters_pulse,
        double time_end, double time_step, double mass_per_time_step)
        : InjectionParameters_ContinuousNoMass{specific_parameters_pulse.time,
                                               time_end, time_step},
          mass_per_time_step{mass_per_time_step} {}

    InjectionParameters_Continuous(
        InjectionParameters_PulseNoMass const &specific_parameters_pulse,
        double time_end, double time_step, double mass_per_time_step)
        : InjectionParameters_Continuous{specific_parameters_pulse.time,
                                         time_end, time_step,
                                         mass_per_time_step} {}

    InjectionParameters_Continuous(
        InjectionParameters_Pulse const &specific_parameters_pulse,
        double time_end, double time_step)
        : InjectionParameters_Continuous{specific_parameters_pulse.time,
                                         time_end, time_step,
                                         specific_parameters_pulse.mass} {}
  };

  std::unique_ptr<SpecificParameters> specific_parameters;
  std::unique_ptr<InjectionParameters> injection_parameters;

  template <typename Geometry, typename TransportParameters,
            typename ReactionParameters>
  InitialConditionParameters_Cases(Directories const &directories,
                                   std::string const &parameter_set_name,
                                   Geometry const &geometry,
                                   TransportParameters const &params_transport,
                                   ReactionParameters const &params_reaction) {
    std::string filename = directories.dir_parameters +
                           "/parameters_initial_condition_" +
                           parameter_set_name + ".dat";
    auto input = io::open_read(filename);
    specific_parameters =
        set_specific_parameters(input, filename, directories, geometry,
                                params_transport, params_reaction);
    injection_parameters = set_injection_parameters(
        input, filename, params_transport, params_reaction);
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  template <typename OStream> static std::ostream &info(OStream &output) {
    output
        << io::line() << "Initial condition parameters\n"
        << io::line()
        << "- Initial condition type:\n"
           "  - point\n"
           "    - All particles at specified position\n"
           "    - Pass on same line:\n"
           "      - position\n"
           "  - uniform_all_cells\n"
           "    - Homogeneous over cell centers throughout the domain\n"
           "  - fluxweighted_all_cells\n"
           "    - Flux-weighted over cell centers throughout the domain\n"
           "  - uniform_patch_faces\n"
           "    - Homogeneous over specified boundary patch face centers\n"
           "    - Pass on same line:\n"
           "      - Patch name(s)\n"
           "  - fluxweighted_patch_faces\n"
           "    - Flux-weighted over specified boundary patch face centers\n"
           "    - Pass on same line:\n"
           "      - Patch name(s)\n"
           "  - uniform_near_patch\n"
           "    - Homogeneous at a fixed distance from specified boundary\n"
           "      patch faces\n"
           "    - Pass on same line:\n"
           "      - Number of tries to place each particle in valid position\n"
           "      - Distance to patch\n"
           "      - Patch name(s)\n"
           "  - fluxweighted_near_patch\n"
           "    - Flux-weighted at a fixed distance from specified boundary\n"
           "      patch faces\n"
           "    - Pass on same line:\n"
           "      - Number of tries to place each particle in valid position\n"
           "      - Distance to patch\n"
           "      - Patch name(s)\n"
           "  - uniform_cartesian_cells\n"
           "    - Homogeneous over cell centers within a Cartesian region\n"
           "    - Pass on same line:\n"
           "      - Number of tries to place each particle in valid position\n"
           "  - fluxweighted_cartesian_cells\n"
           "    - Pass on same line:\n"
           "      - Lower and upper region boundaries along each Cartesian\n"
           "        dimension\n"
           "  - prescribed_positions_masses\n"
           "    - Prescribed particle positions\n"
           "    - Pass on same line:\n"
           "      - Name of file with positions (one particle per line)\n"
           "      - Path to file (optional [Parameters directory])\n"
           "  - prescribed_positions_masses\n"
           "    - Prescribed particle positions and masses\n"
           "    - Pass on same line:\n"
           "      - Name of file with positions and masses (one particle per\n"
           "        line)\n"
           "      - Path to file (optional [Parameters directory])\n"
           "  - prescribed_positions_masses_tags\n"
           "    - Prescribed particle positions and masses\n"
           "    - Pass on same line:\n"
           "      - Name of file with particle positions, masses, and tags\n"
           "        (one particle per line)\n"
           "        per line\n"
           "      - Path to file (optional [Parameters directory])\n"
           "  - prescribed_positions_masses_tags_times\n"
           "    - Prescribed particle positions and masses\n"
           "    - Pass on same line:\n"
           "      - File name for file with particle positions, masses, tags,\n"
           "        and times (one particle per line)\n"
           "      - Path to file (optional [Parameters directory])\n"
           "  - uniform_cylinder\n"
           "    - Homogeneous inside cylinder\n"
           "    - Pass on same line:\n"
           "      - Number of tries to place each particle in valid position\n"
           "      - Base centers\n"
           "      - Radius\n"
           "      - Excluded inner radius (optional [0.])\n"
           "  - uniform_parallelipiped\n"
           "    - Homogeneous inside parallelipiped\n"
           "    - Pass on same line:\n"
           "      - Number of tries to place each particle in valid position\n"
           "      - Corner (side vector origin)\n"
           "      - Side vectors (up to one per spatial dimension, degenerate\n"
           "        if fewer)\n"
           "  - uniform_sphere\n"
           "    - Homogeneous inside sphere\n"
           "    - Pass on same line:\n"
           "      - Number of tries to place each particle in valid position\n"
           "      - Centers\n"
           "      - Radius\n"
           "      - Excluded inner radius (optional [0.])\n"
           "- Injection continuity type:\n"
           "  (Note:\n"
           "    - Line ignored if initial condition type prescribes times\n"
           "    - Total injected mass ignored if initial condition type\n"
           "      prescribes masses)\n"
           "  - pulse\n"
           "    - Instantaneous pulse injection\n"
           "    - Pass on same line:\n"
           "      - Time units for injection time:\n"
           "        - diffusion\n"
           "          - Diffusion time units\n"
           "        - advection\n"
           "          - Advection time units\n"
           "        - reaction\n"
           "          - Reaction time units\n"
           "        - arbitrary\n"
           "          - Arbitary units (no rescaling)\n"
           "      - Injection time\n"
           "      - Total injected mass\n"
           "  - continuous\n"
           "    - Continuous injection\n"
           "    - Pass on same line:\n"
           "      - Time units for injection-related times:\n"
           "        - diffusion\n"
           "          - Diffusion time units\n"
           "        - advection\n"
           "          - Advection time units\n"
           "        - reaction\n"
           "          - Reaction time units\n"
           "        - arbitrary\n"
           "          - Arbitary units (no rescaling)\n"
           "      - Injection start time\n"
           "      - Injection end time\n"
           "      - Injection discretization option:\n"
           "        - time_step\n"
           "          - Prescribed injection time step\n"
           "          - Pass on same line:\n"
           "            - Injection time step\n"
           "        - accuracy\n"
           "          - Choose the minimum time step based on prescribed\n"
           "            accuracies\n"
           "          - Pass on same line:\n"
           "            - Injection time step accuracy with respect to\n"
           "              advection time\n"
           "            - Injection time step accuracy with respect to\n"
           "              diffusion time\n"
           "            - Injection time step accuracy with respect to\n"
           "              reaction time\n"
           "      - Total injected mass per injection time step\n"
        << io::line();
    return output;
  }

private:
  /** \brief Read initial condition type and specific parameters. */
  template <typename Geometry, typename TransportParameters,
            typename ReactionParameters, typename IStream>
  std::unique_ptr<SpecificParameters>
  set_specific_parameters(IStream &input, std::string const &filename,
                          Directories const &directories,
                          Geometry const &geometry,
                          TransportParameters const &params_transport,
                          ReactionParameters const &params_reaction) {
    std::string in_file = std::string{"In file "} + filename + " : ";
    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse initial condition type", name);
    if (!InitialConditionList{}.contains(name))
      throw std::runtime_error{"Initial condition type " + name +
                               " not supported"};
    type = InitialConditionList::type(name);
    std::string for_initial_condition_type =
        "Initial condition type " + name + " : ";
    switch (type) {
    case (InitialConditionList::Type::point): {
      return std::make_unique<SpecificParameters_Position>(io::read<double>(
          split_line, param_index, Geometry::dim,
          in_file + for_initial_condition_type + "Could not parse position"));
    }
    case (InitialConditionList::Type::uniform_patch_faces):
    case (InitialConditionList::Type::fluxweighted_patch_faces): {
      std::size_t nr_patches = split_line.size() - param_index;
      if (nr_patches == 0)
        throw std::runtime_error{
            in_file + for_initial_condition_type +
            "Could not parse patch names(s) : At least one patch required"};
      return std::make_unique<SpecificParameters_Patches>(io::read<std::string>(
          split_line, param_index, split_line.size() - param_index,
          in_file + for_initial_condition_type +
              "Could not parse patch names"));
    }
    case (InitialConditionList::Type::uniform_near_patch):
    case (InitialConditionList::Type::fluxweighted_near_patch): {
      auto nr_tries =
          io::read<std::size_t>(split_line, param_index,
                                in_file + for_initial_condition_type +
                                    "Could not parse number of tries");
      auto distance = io::read<double>(split_line, param_index,
                                       in_file + for_initial_condition_type +
                                           "Could not parse distance to patch");
      std::size_t nr_patches = split_line.size() - param_index;
      if (nr_patches == 0)
        throw std::runtime_error{
            in_file + for_initial_condition_type +
            "Could not parse patch names(s) : At least one patch required"};
      auto patch_names =
          io::read<std::string>(split_line, param_index, nr_patches,
                                in_file + for_initial_condition_type +
                                    "Could not parse patch name(s)");
      return std::make_unique<SpecificParameters_Patches_Distance>(
          patch_names, distance, nr_tries);
    }
    case (InitialConditionList::Type::uniform_cartesian_cells):
    case (InitialConditionList::Type::fluxweighted_cartesian_cells): {
      return std::make_unique<SpecificParameters_Boundaries>(
          io::read<std::pair<double, double>>(
              split_line, param_index, Geometry::dim,
              in_file + for_initial_condition_type +
                  "Could not parse region boundaries"));
    }
    case (InitialConditionList::Type::prescribed_positions):
    case (InitialConditionList::Type::prescribed_positions_masses):
    case (InitialConditionList::Type::prescribed_positions_masses_tags):
    case (InitialConditionList::Type::prescribed_positions_masses_tags_times): {
      auto filename_data = io::read<std::string>(
          split_line, param_index,
          in_file + for_initial_condition_type + "Could not parse file name");
      auto dir = param_index < split_line.size()
                     ? io::expand_env(io::expand_home_dir(io::read<std::string>(
                           split_line, param_index,
                           in_file + for_initial_condition_type +
                               "Could not file directory")))
                     : directories.dir_parameters;
      return std::make_unique<SpecificParameters_Filename>(dir + "/" +
                                                           filename_data);
    }
    case (InitialConditionList::Type::uniform_cylinder): {
      auto nr_tries =
          io::read<std::size_t>(split_line, param_index,
                                in_file + for_initial_condition_type +
                                    "Could not parse number of tries");
      std::vector<std::vector<double>> centers;
      centers.reserve(Geometry::dim);
      for (std::size_t ii = 0; ii < Geometry::dim; ++ii)
        centers.push_back(io::read<double>(
            split_line, param_index, Geometry::dim,
            in_file + for_initial_condition_type + "Could not parse " +
                io::ordinal(ii) + " center position"));
      auto radius = io::read<double>(split_line, param_index,
                                     in_file + for_initial_condition_type +
                                         "Could not parse radius");
      auto inner_radius =
          param_index < split_line.size()
              ? io::read<double>(split_line, param_index,
                                 in_file + for_initial_condition_type +
                                     "Could not parse inner radius")
              : 0.;
      return std::make_unique<SpecificParameters_Cylinder>(
          centers[0], centers[1], radius, nr_tries, inner_radius);
    }
    case (InitialConditionList::Type::uniform_parallelipiped): {
      auto nr_tries =
          io::read<std::size_t>(split_line, param_index,
                                in_file + for_initial_condition_type +
                                    "Could not parse number of tries");
      auto corner = io::read<double>(split_line, param_index, Geometry::dim,
                                     in_file + for_initial_condition_type +
                                         "Could not parse corner");
      std::vector<std::vector<double>> sides;
      sides.reserve(Geometry::dim);
      for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
        sides.push_back(io::read<double>(split_line, param_index, Geometry::dim,
                                         in_file + for_initial_condition_type +
                                             "Could not parse " +
                                             io::ordinal(dd) + " side"));
      return std::make_unique<SpecificParameters_Parallelipiped>(corner, sides,
                                                                 nr_tries);
    }
    case (InitialConditionList::Type::uniform_sphere): {
      auto nr_tries =
          io::read<std::size_t>(split_line, param_index,
                                in_file + for_initial_condition_type +
                                    "Could not parse number of tries");
      auto center = io::read<double>(split_line, param_index, Geometry::dim,
                                     in_file + for_initial_condition_type +
                                         "Could not parse center");
      auto radius = io::read<double>(split_line, param_index,
                                     in_file + for_initial_condition_type +
                                         " : "
                                         "Could not parse radius");
      auto inner_radius =
          param_index < split_line.size()
              ? io::read<double>(split_line, param_index,
                                 in_file + for_initial_condition_type +
                                     "Could not parse inner radius")
              : 0.;
      return std::make_unique<SpecificParameters_Sphere>(
          center, radius, nr_tries, inner_radius);
    }
    default:
      return std::make_unique<SpecificParameters>();
    }
  }

  /** \brief Read initial condition info from input stream. */
  template <typename IStream, typename TransportParameters,
            typename ReactionParameters>
  std::unique_ptr<InjectionParameters>
  set_injection_parameters(IStream &input, std::string const &filename,
                           TransportParameters const &params_transport,
                           ReactionParameters const &params_reaction) {
    if (type ==
        InitialConditionList::Type::prescribed_positions_masses_tags_times) {
      continuity = "prescribed";
      return std::make_unique<InjectionParameters>();
    }

    std::string in_file = std::string{"In file "} + filename + " : ";
    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             std::string{"In file "} + filename +
                 " : "
                 "Could not parse injection continuity type",
             continuity);
    std::string for_injection_continuity_type =
        std::string{"Injection continuity type "} + continuity + " : ";

    auto time_units =
        io::read<std::string>(split_line, param_index,
                              in_file + for_injection_continuity_type +
                                  "Could not parse injection time units");
    if (!TimeUnits{}.contains(time_units))
      throw std::runtime_error{in_file + for_injection_continuity_type +
                               "Time units " + time_units + " not supported"};
    double time_unit_factor =
        ptof::time_unit_factor(time_units, params_transport, params_reaction);

    auto injection_time_start =
        io::read<double>(split_line, param_index,
                         in_file + for_injection_continuity_type +
                             "Could not parse injection start time");
    injection_time_start *= time_unit_factor;

    if (continuity == "pulse") {
      if (type != InitialConditionList::Type::prescribed_positions_masses &&
          type !=
              InitialConditionList::Type::prescribed_positions_masses_tags) {
        return std::make_unique<InjectionParameters_Pulse>(
            injection_time_start,
            io::read<double>(split_line, param_index,
                             in_file + for_injection_continuity_type +
                                 "Could not parse total mass"));
      }
      return std::make_unique<InjectionParameters_PulseNoMass>(
          injection_time_start);
    } else if (continuity == "continuous") {
      auto injection_time_end =
          io::read<double>(split_line, param_index,
                           in_file + for_injection_continuity_type +
                               "Could not parse injection end time");
      injection_time_end *= time_unit_factor;
      auto discretization_option = io::read<std::string>(
          split_line, param_index,
          in_file + for_injection_continuity_type +
              "Could not parse time discretization option");
      double injection_time_step;
      if (discretization_option == "time_step") {
        io::read(split_line, param_index,
                 in_file + for_injection_continuity_type +
                     "Could not parse discretization time step",
                 injection_time_step);
        injection_time_step *= time_unit_factor;
      } else if (discretization_option == "accuracy") {
        double time_step_adv, time_step_diff, time_step_react;
        io::read(split_line, param_index,
                 in_file + for_injection_continuity_type +
                     "Could not parse time step in units of advective, "
                     "diffusive, and reactive time",
                 time_step_adv, time_step_diff, time_step_react);
        injection_time_step =
            std::min({time_step_adv * params_transport.advection_time,
                      time_step_diff * params_transport.diffusion_time,
                      time_step_react * params_reaction.reaction_time});
      } else {
        throw std::runtime_error{std::string{"Injection continuity type "} +
                                 continuity + " not supported"};
      }

      if (type != InitialConditionList::Type::prescribed_positions_masses &&
          type !=
              InitialConditionList::Type::prescribed_positions_masses_tags) {
        return std::make_unique<InjectionParameters_Continuous>(
            injection_time_start, injection_time_end, injection_time_step,
            io::read<double>(
                split_line, param_index,
                in_file + for_injection_continuity_type +
                    "Could not parse total mass per discretization time step"));
      }
      return std::make_unique<InjectionParameters_ContinuousNoMass>(
          injection_time_start, injection_time_end, injection_time_step);
    } else {
      throw std::runtime_error{std::string{"Injection continuity type "} +
                               continuity + " not supported"};
    }
  }
};

/**
   \class InitialCondition_Cases PTOF/InitialCondition_Cases.h
   "PTOF/InitialCondition_Cases.h"
   \brief InitialCondition object to handle implemented initial condition types.
*/
template <typename Particle, typename Geometry> class InitialCondition_Cases {
public:
  using Parameters =
      InitialConditionParameters_Cases; /**< Initial condition parameters.*/

  std::size_t const
      nr_particles; /**< Number of particles per injection step. */
  InitialConditionList::Type type; /**< Type of initial condition. */

  /**
     \brief Constructor.
     \param geometry Domain geometry info and utilities.
     \param velocity_field Velocity field as a function of state.
     \param nr_particles Number of particles to make.
     \param parameters Output parameters.
     \param mask Scalar field.
     \param threshold Threshold for the mask, such that cells where the mask is
     above or equal to the threshold are considered.
  */
  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(Geometry const &geometry,
                         VelocityField const &velocity_field,
                         std::size_t nr_particles, Parameters const &parameters,
                         Mask const &mask = {}, double threshold = 0.)
      : nr_particles{nr_particles}, type{parameters.type},
        _particle_maker{make_particle_maker(geometry, parameters)},
        _initial_condition{
            set_type(geometry, velocity_field, parameters, mask, threshold)} {}

  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(Geometry &&, VelocityField const &, std::size_t,
                         Parameters const &, Mask const & = {},
                         double = 0.) = delete;

  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(meta::Selector_t<Particle>, Geometry const &geometry,
                         VelocityField const &velocity_field,
                         std::size_t nr_particles, Parameters const &parameters,
                         Mask const &mask = {}, double threshold = 0.)
      : nr_particles{nr_particles}, type{parameters.type},
        _particle_maker{make_particle_maker(geometry, parameters)},
        _initial_condition{
            set_type(geometry, velocity_field, parameters, mask, threshold)} {}

  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  InitialCondition_Cases(meta::Selector_t<Particle>, Geometry &&,
                         VelocityField const &, std::size_t, Parameters const &,
                         Mask const & = {}, double = 0.) = delete;

  /**
     \brief Make particles according to prescribed initial condition and
     particle maker, with \c nr_particles particles as specified in \c
     parameters.
     \return Container with particles.
  */
  auto operator()() { return make_particles(nr_particles); }

  /**
     \brief Make particles.
     \param nr_particles Number of particles to make.
     \return Container with particles.
  */
  auto operator()(std::size_t nr_particles) {
    return make_particles(nr_particles);
  }

  /**
     \brief Make a single particle.
     \return One particle.
  */
  auto make_particle() { return _initial_condition->make_particle(); }

  /**
     \brief Make a single position.
     \return Position.
  */
  auto make_position() { return _initial_condition->make_position(); }

  /**
     \brief Make a single position along with location hint.
     \return Position and hint.
  */
  auto make_position_and_hint() {
    return _initial_condition->make_position_and_hint();
  }

  /**
     \brief Make particles.
     \param nr_particles Number of particles to make.
     \return Container with particles.
  */
  auto make_particles(std::size_t nr_particles) {
    return _initial_condition->make_particles(nr_particles);
  };

  /**
     \brief Make positions.
     \param nr_positions Number of positions to make.
     \return Container with positions.
  */
  auto make_positions(std::size_t nr_positions) {
    return _initial_condition->make_positions(nr_particles);
  };

private:
  template <typename TT, bool periodic = false> struct BP {
    using type = meta::Empty;
  };

  template <typename TT> struct BP<TT, true> {
    using type = typename TT::BoundaryPeriodic;
  };

  static constexpr bool periodic =
      meta::has_periodicity_v<typename Particle::State>;

  using BoundaryPeriodic =
      typename std::conditional_t<periodic, BP<Geometry, true>,
                                  BP<Geometry, false>>::type;

  using ParticleMaker = std::conditional_t<
      periodic,
      ParticleMaker_Periodic<Particle, typename Geometry::Locator const &,
                             BoundaryPeriodic const &>,
      ParticleMaker_Generic<Particle, typename Geometry::Locator const &>>;
  using IC = InitialCondition<ParticleMaker &, Geometry const &>;

  ParticleMaker _particle_maker;
  std::unique_ptr<IC>
      _initial_condition; /**< Pointer to initial condition object. */

  ParticleMaker make_particle_maker(Geometry const &geometry,
                                    Parameters const &parameters) {
    double initial_particle_mass = 0.;
    double initial_time = 0.;
    if (parameters.continuity == "pulse" &&
        type != InitialConditionList::Type::prescribed_positions_masses &&
        type != InitialConditionList::Type::prescribed_positions_masses_tags) {
      using InjectionParameters = Parameters::InjectionParameters_Pulse;
      InjectionParameters *injection_parameters =
          cast_injection<InjectionParameters>(parameters);
      initial_particle_mass = injection_parameters->mass / nr_particles;
      initial_time = injection_parameters->time;
    }
    if constexpr (periodic)
      return ParticleMaker{meta::Selector_t<Particle>{}, geometry.locator,
                           geometry.boundary_periodic, initial_time,
                           initial_particle_mass};
    else
      return ParticleMaker{meta::Selector_t<Particle>{}, geometry.locator,
                           initial_time, initial_particle_mass};
  }

  template <typename VelocityField, typename Mask = meta::Empty>
  std::unique_ptr<IC> set_type(Geometry const &geometry,
                               VelocityField const &velocity_field,
                               Parameters const &parameters,
                               Mask const &mask = {}, double threshold = 0.) {
    switch (type) {
    case InitialConditionList::Type::point: {
      using InitialCondition =
          InitialCondition_Point<ParticleMaker &, Geometry const &>;
      using SpecificParameters = Parameters::SpecificParameters_Position;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      return set_injection(InitialCondition{_particle_maker, geometry,
                                            ParticleMaker::State::make_position(
                                                specific_parameters->position)},
                           parameters);
    }
    case InitialConditionList::Type::uniform_all_cells: {
      using InitialCondition =
          InitialCondition_UniformCellCenters<ParticleMaker &, Geometry const &,
                                              std::vector<Foam::label>>;
      return set_injection(
          InitialCondition{_particle_maker, geometry,
                           cell_ids_masked(all_cell_ids(geometry.mesh()),
                                           geometry, mask, threshold)},
          parameters);
    }
    case InitialConditionList::Type::fluxweighted_all_cells: {
      using InitialCondition = InitialCondition_FluxweightedCellCenters<
          ParticleMaker &, Geometry const &, std::vector<Foam::label>>;
      if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
        return set_injection(
            InitialCondition{_particle_maker, geometry,
                             cell_ids_masked(all_cell_ids(geometry.mesh()),
                                             geometry, mask, threshold),
                             velocity_field},
            parameters);
      else
        throw std::runtime_error{std::string("Initial condition type ") +
                                 InitialConditionList::name(type) + ": " +
                                 "Velocity field not provided"};
    }
    case InitialConditionList::Type::uniform_patch_faces: {
      using InitialCondition =
          InitialCondition_UniformFaceCenters<ParticleMaker &, Geometry const &,
                                              std::vector<Foam::label>>;
      using SpecificParameters = Parameters::SpecificParameters_Patches;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      return set_injection(
          InitialCondition{
              _particle_maker, geometry,
              patch_face_ids_masked(specific_parameters->patch_names, geometry,
                                    mask, threshold)},
          parameters);
    }
    case InitialConditionList::Type::fluxweighted_patch_faces: {
      using InitialCondition = InitialCondition_FluxweightedFaceCenters<
          ParticleMaker &, Geometry const &, std::vector<Foam::label>>;
      using SpecificParameters = Parameters::SpecificParameters_Patches;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
        return set_injection(
            InitialCondition{
                _particle_maker, geometry,
                patch_face_ids_masked(specific_parameters->patch_names,
                                      geometry, mask, threshold),
                velocity_field},
            parameters);
      else
        throw std::runtime_error{std::string("Initial condition type ") +
                                 InitialConditionList::name(type) + ": " +
                                 "Velocity field not provided"};
    }
    case InitialConditionList::Type::uniform_near_patch: {
      using InitialCondition =
          InitialCondition_UniformNearFaces<ParticleMaker &, Geometry const &,
                                            std::vector<Foam::label>>;
      using SpecificParameters =
          Parameters::SpecificParameters_Patches_Distance;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      return set_injection(
          InitialCondition{
              _particle_maker, geometry,
              patch_face_ids_masked(specific_parameters->patch_names, geometry,
                                    mask, threshold),
              specific_parameters->distance, specific_parameters->nr_tries},
          parameters);
    }
    case InitialConditionList::Type::fluxweighted_near_patch: {
      using InitialCondition = InitialCondition_FluxweightedNearFaces<
          ParticleMaker &, Geometry const &, std::vector<Foam::label>>;
      using SpecificParameters =
          Parameters::SpecificParameters_Patches_Distance;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
        return set_injection(
            InitialCondition{
                _particle_maker, geometry,
                patch_face_ids_masked(specific_parameters->patch_names,
                                      geometry, mask, threshold),
                velocity_field, specific_parameters->distance,
                specific_parameters->nr_tries},
            parameters);
      else
        throw std::runtime_error{std::string("Initial condition type ") +
                                 InitialConditionList::name(type) + ": " +
                                 "Velocity field not provided"};
    }
    case InitialConditionList::Type::uniform_cartesian_cells: {
      using InitialCondition =
          InitialCondition_UniformCellCenters<ParticleMaker &, Geometry const &,
                                              std::vector<Foam::label>>;
      using SpecificParameters = Parameters::SpecificParameters_Boundaries;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      return set_injection(
          InitialCondition{
              _particle_maker, geometry,
              cell_ids_masked(
                  cell_ids_region_cartesian(
                      specific_parameters->region_boundaries, geometry.locator),
                  geometry, mask, threshold)},
          parameters);
    }
    case InitialConditionList::Type::fluxweighted_cartesian_cells: {
      using InitialCondition = InitialCondition_FluxweightedCellCenters<
          ParticleMaker &, Geometry const &, std::vector<Foam::label>>;
      using SpecificParameters = Parameters::SpecificParameters_Boundaries;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
        return set_injection(
            InitialCondition{
                _particle_maker, geometry,
                cell_ids_masked(cell_ids_region_cartesian(
                                    specific_parameters->region_boundaries,
                                    geometry.locator),
                                geometry, mask, threshold),
                velocity_field},
            parameters);
      else
        throw std::runtime_error{std::string("Initial condition type ") +
                                 InitialConditionList::name(type) + ": " +
                                 "Velocity field not provided"};
    }
    case InitialConditionList::Type::prescribed_positions: {
      using InitialCondition =
          InitialCondition_PrescribedPositions<ParticleMaker &,
                                               Geometry const &>;
      using SpecificParameters = Parameters::SpecificParameters_Filename;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      return set_injection(InitialCondition{_particle_maker, geometry,
                                            specific_parameters->filename},
                           parameters);
    }
    case InitialConditionList::Type::prescribed_positions_masses: {
      using InitialCondition =
          InitialCondition_PrescribedPositionsMasses<ParticleMaker &,
                                                     Geometry const &>;
      using SpecificParameters = Parameters::SpecificParameters_Filename;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      return set_injection(InitialCondition{_particle_maker, geometry,
                                            specific_parameters->filename},
                           parameters);
    }
    case InitialConditionList::Type::prescribed_positions_masses_tags: {
      using InitialCondition =
          InitialCondition_PrescribedPositionsMassesTags<ParticleMaker &,
                                                         Geometry const &>;
      using SpecificParameters = Parameters::SpecificParameters_Filename;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      return set_injection(InitialCondition{_particle_maker, geometry,
                                            specific_parameters->filename},
                           parameters);
    }
    case InitialConditionList::Type::prescribed_positions_masses_tags_times: {
      using InitialCondition =
          InitialCondition_PrescribedPositionsMassesTagsTimes<ParticleMaker &,
                                                              Geometry const &>;
      using SpecificParameters = Parameters::SpecificParameters_Filename;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      return set_injection(InitialCondition{_particle_maker, geometry,
                                            specific_parameters->filename},
                           parameters);
    }
    case InitialConditionList::Type::uniform_cylinder: {
      using Distribution = stochastic::uniform_cylinder_distribution;
      using InitialCondition =
          InitialCondition_DistributedPosition<ParticleMaker &,
                                               Geometry const &, Distribution>;
      using SpecificParameters = Parameters::SpecificParameters_Cylinder;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      return set_injection(
          InitialCondition{_particle_maker, geometry,
                           Distribution{specific_parameters->center_1,
                                        specific_parameters->center_2,
                                        specific_parameters->radius,
                                        specific_parameters->inner_radius},
                           specific_parameters->nr_tries},
          parameters);
    }
    case InitialConditionList::Type::uniform_parallelipiped: {
      using Distribution = stochastic::uniform_parallelipiped_distribution;
      using InitialCondition =
          InitialCondition_DistributedPosition<ParticleMaker &,
                                               Geometry const &, Distribution>;
      using SpecificParameters = Parameters::SpecificParameters_Parallelipiped;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      return set_injection(
          InitialCondition{_particle_maker, geometry,
                           Distribution{specific_parameters->corner,
                                        specific_parameters->sides},
                           specific_parameters->nr_tries},
          parameters);
    }
    case InitialConditionList::Type::uniform_sphere: {
      using Distribution = stochastic::uniform_sphere_distribution;
      using InitialCondition =
          InitialCondition_DistributedPosition<ParticleMaker &,
                                               Geometry const &, Distribution>;
      using SpecificParameters = Parameters::SpecificParameters_Sphere;
      SpecificParameters *specific_parameters =
          cast<SpecificParameters>(parameters);
      return set_injection(
          InitialCondition{_particle_maker, geometry,
                           Distribution{specific_parameters->center,
                                        specific_parameters->radius,
                                        specific_parameters->inner_radius},
                           specific_parameters->nr_tries},
          parameters);
    }

    default:
      throw std::runtime_error{std::string{"Initial condition type "} +
                               InitialConditionList::name(parameters.type) +
                               " not supported"};
    }
  }

  template <typename Container, typename Mask = meta::Empty>
  auto cell_ids_masked(Container cell_ids, Geometry const &geometry,
                       Mask const &mask = {}, double threshold = 0.) {
    if constexpr (!std::is_same_v<Mask, meta::Empty>)
      apply_mask_cells_inplace(cell_ids, mask, threshold);
    return cell_ids;
  }

  template <typename Container, typename Mask = meta::Empty>
  auto patch_face_ids_masked(Container face_ids, Geometry const &geometry,
                             Mask const &mask = {}, double threshold = 0.) {
    if constexpr (!std::is_same_v<Mask, meta::Empty>)
      apply_mask_patch_faces_inplace(face_ids, geometry.mesh(), mask,
                                     threshold);
    return face_ids;
  }

  template <typename Mask = meta::Empty>
  auto patch_face_ids_masked(std::vector<std::string> const &patch_names,
                             Geometry const &geometry, Mask const &mask = {},
                             double threshold = 0.) {
    return patch_face_ids_masked(patch_face_ids(patch_names, geometry.mesh()),
                                 geometry);
  }

  template <typename InitialCondition>
  std::unique_ptr<IC> set_injection(InitialCondition ic,
                                    Parameters const &parameters) {
    if (parameters.continuity == "continuous") {
      using InjectionParameters =
          Parameters::InjectionParameters_ContinuousNoMass;
      InjectionParameters *injection_parameters =
          cast_injection<InjectionParameters>(parameters);
      return std::make_unique<InitialCondition_Continuous<InitialCondition>>(
          std::move(ic), injection_parameters->time,
          injection_parameters->time_end, injection_parameters->time_step);
    }
    return std::make_unique<InitialCondition>(std::move(ic));
  }

  template <typename SpecificParameters>
  SpecificParameters *cast(Parameters const &parameters) {
    return dynamic_cast<SpecificParameters *>(
        parameters.specific_parameters.get());
  }

  template <typename InjectionParameters>
  InjectionParameters *cast_injection(Parameters const &parameters) {
    return dynamic_cast<InjectionParameters *>(
        parameters.injection_parameters.get());
  }

public:
  /**
     \brief Output information about current object.
     \param output Output stream.
  */
  template <typename OStream>
  std::ostream &info_runtime(OStream &output) const {
    output << io::line() << "Initial condition\n"
           << io::line() << "Type: " + InitialConditionList::name(type) + "\n"
           << io::line();
    return output;
  }
};
template <typename Particle, typename Geometry, typename VelocityField,
          typename Mask>
InitialCondition_Cases(meta::Selector_t<Particle>, Geometry const &,
                       VelocityField const &, std::size_t nr_particles,
                       InitialConditionParameters_Cases const &, Mask const &,
                       double) -> InitialCondition_Cases<Particle, Geometry>;
template <typename Particle, typename Geometry, typename VelocityField,
          typename Mask>
InitialCondition_Cases(meta::Selector_t<Particle>, Geometry const &,
                       VelocityField const &,
                       InitialConditionParameters_Cases const &, Mask const &)
    -> InitialCondition_Cases<Particle, Geometry>;
template <typename Particle, typename Geometry, typename VelocityField>
InitialCondition_Cases(meta::Selector_t<Particle>, Geometry const &,
                       VelocityField const &, std::size_t nr_particles,
                       InitialConditionParameters_Cases const &)
    -> InitialCondition_Cases<Particle, Geometry>;
} // namespace ptof

#endif /* PTOF_INITIALCONDITION_CASES_H */
