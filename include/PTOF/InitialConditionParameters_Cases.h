/**
   \file PTOF/InitialConditionParameters_Cases.h
   \author Tomás Aquino
   \date 08/05/2025
   \brief Parameters for initial condition.
*/

#ifndef PTOF_INITIALCONDITIONPARAMETERS_CASES_H
#define PTOF_INITIALCONDITIONPARAMETERS_CASES_H

#include "General/IO.h"
#include "PTOF/Directories.h"
#include "PTOF/InitialConditionList.h"
#include "PTOF/TimeUnits.h"
#include <ostream>
#include <stdexcept>
#include <string>
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

    virtual std::ostream &info_runtime(std::ostream &output) const {
      return output;
    }
  };

  struct SpecificParameters_Patches : public SpecificParameters {
    std::vector<std::string> patch_names;

    SpecificParameters_Patches(std::vector<std::string> patch_names)
        : patch_names{patch_names} {}

    std::ostream &info_runtime(std::ostream &output) const override {
      SpecificParameters::info_runtime(output);
      output << "Patches:\n";
      for (auto const &patch : patch_names)
        output << "  - " << patch << "\n";
      return output;
    }
  };

  struct SpecificParameters_Position : public SpecificParameters {
    std::vector<double> position;

    SpecificParameters_Position(std::vector<double> position)
        : position{position} {}

    std::ostream &info_runtime(std::ostream &output) const override {
      SpecificParameters::info_runtime(output);
      output << "Position: ";
      io::print(output, position, false, " ") << "\n";
      return output;
    }
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

    std::ostream &info_runtime(std::ostream &output) const override {
      SpecificParameters_Patches::info_runtime(output);
      output << "Distance: " << distance << "\n";
      output << "Number of tries: " << nr_tries << "\n";
      return output;
    }
  };

  struct SpecificParameters_Filename : public SpecificParameters {
    std::string filename;

    SpecificParameters_Filename(std::string filename) : filename{filename} {}

    std::ostream &info_runtime(std::ostream &output) const override {
      SpecificParameters::info_runtime(output);
      output << "Filename: " << filename << "\n";
      return output;
    }
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

    std::ostream &info_runtime(std::ostream &output) const override {
      SpecificParameters::info_runtime(output);
      output << "Region boundaries:";
      for (auto const &bounds : region_boundaries)
        output << " (" << bounds.first << " " << bounds.second << ")";
      output << "\n";
      return output;
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

    std::ostream &info_runtime(std::ostream &output) const override {
      SpecificParameters::info_runtime(output);
      output << "Center 1: ";
      io::print(output, center_1, false, " ") << "\n";
      output << "Center 2: ";
      io::print(output, center_2, false, " ") << "\n";
      output << "Radius: " << radius << "\n"
             << "Number of tries: " << nr_tries << "\n";
      if (inner_radius != 0.)
        output << "Inner radius: " << inner_radius << "\n";
      return output;
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

    std::ostream &info_runtime(std::ostream &output) const override {
      SpecificParameters::info_runtime(output);
      output << "Corner: ";
      io::print(output, corner, false, " ") << "\n";
      output << "Sides:";
      for (auto const &side : sides) {
        output << " (";
        io::print(output, side, false, " ");
        output << ")";
      }
      output << "\n";
      output << "Number of tries: " << nr_tries;
      return output;
    }
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

    std::ostream &info_runtime(std::ostream &output) const override {
      SpecificParameters::info_runtime(output);
      output << "Center: ";
      io::print(output, center, false, " ") << "\n";
      output << "Radius: " << radius << "\n"
             << "Number of tries: " << nr_tries << "\n";
      if (inner_radius != 0.)
        output << "Inner radius: " << inner_radius << "\n";
      return output;
    }
  };

  struct InjectionParameters {
    virtual ~InjectionParameters() {}

    virtual std::ostream &info_runtime(std::ostream &output) const {
      return output;
    }
  };

  struct InjectionParameters_PulseNoMass : public InjectionParameters {
    double time;

    InjectionParameters_PulseNoMass(double time) : time{time} {}

    std::ostream &info_runtime(std::ostream &output) const override {
      InjectionParameters::info_runtime(output);
      output << "Injection time: " << time << "\n";
      return output;
    }
  };

  struct InjectionParameters_Pulse : public InjectionParameters_PulseNoMass {
    double mass;

    InjectionParameters_Pulse(double time, double mass)
        : InjectionParameters_PulseNoMass{time}, mass{mass} {}

    InjectionParameters_Pulse(
        InjectionParameters_PulseNoMass const &specific_parameters_pulse,
        double mass)
        : InjectionParameters_Pulse{specific_parameters_pulse.time, mass} {}

    std::ostream &info_runtime(std::ostream &output) const override {
      InjectionParameters_PulseNoMass::info_runtime(output);
      output << "Injection mass: " << mass << "\n";
      return output;
    }
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

    std::ostream &info_runtime(std::ostream &output) const override {
      InjectionParameters_PulseNoMass::info_runtime(output);
      output << "Injection end time: " << time_end << "\n"
             << "Injection time step: " << time_step << "\n";
      return output;
    }
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

    std::ostream &info_runtime(std::ostream &output) const override {
      InjectionParameters_ContinuousNoMass::info_runtime(output);
      output << "Mass per time step: " << time_end << "\n";
      return output;
    }
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
  static std::ostream &info(std::ostream &output) {
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
} // namespace ptof

#endif /* PTOF_INITIALCONDITIONPARAMETERS_CASES_H */
