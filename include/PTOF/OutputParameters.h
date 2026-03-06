/**
   \file PTOF/OutputParameters.h
   \author Tomas Aquino
   \date 08/05/2025
   \brief Parameters for output.
*/

#ifndef PTOF_OUTPUTPARAMETERS_H
#define PTOF_OUTPUTPARAMETERS_H

#include "General/IO.h"
#include "PTOF/Directories.h"
#include "PTOF/EndCriterionList.h"
#include "PTOF/MeasurementList.h"
#include "PTOF/MeasurementSpacingList.h"
#include "PTOF/TimeUnitsList.h"
#include <cstddef>
#include <fstream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace ptof {
/**
   \struct OutputParameters_Cases PTOF/OutputParameters.h
   "PTOF/OutputParameters.h"
   \brief Parameters for output.
*/
struct OutputParameters_Cases {
public:
  std::string time_units;
  double time_unit_factor;
  std::string end_criterion_logical_combination;
  std::vector<std::string> end_criteria;
  std::unordered_map<EndCriterionList::Type, double> end_values;
  std::string measurement_spacing;
  double time_min;
  double time_max;
  double time_increment;

  struct Measurement {
    Measurement(std::string name, int precision = 8)
        : name{name}, precision{precision} {}

    virtual ~Measurement() = default;

    std::string name;
    int precision;

    std::ostream &info_runtime(std::ostream &output) const {
      output << "  - " << name << "\n"
             << "    - Precision: " << precision << "\n";
      return output;
    }
  };

  struct Measurement_Field : public Measurement {
    Measurement_Field(std::string name, std::string field_name,
                      int precision = 8)
        : Measurement{name, precision}, field_name{field_name} {}

    std::string field_name;

    std::ostream &info_runtime(std::ostream &output) const {
      Measurement::info_runtime(output);
      output << "    - Field: " << field_name << "\n";
      return output;
    }
  };

  struct Measurement_Dim_Position : public Measurement {
    Measurement_Dim_Position(std::string name, std::size_t dim, double position,
                             int precision = 8)
        : Measurement{name, precision}, dim{dim}, position{position} {}

    std::size_t dim;
    double position;

    std::ostream &info_runtime(std::ostream &output) const {
      Measurement::info_runtime(output);
      output << "    - Dimension: " << dim << "\n"
             << "    - Position: " << position << "\n";
      return output;
    }
  };

  struct Measurement_Order : public Measurement {
    Measurement_Order(std::string name, double order, int precision = 8)
        : Measurement{name, precision}, order{order} {}

    double order;

    std::ostream &info_runtime(std::ostream &output) const {
      Measurement::info_runtime(output);
      output << "    - Order: " << order << "\n";
      return output;
    }
  };

  struct Measurement_Orders : public Measurement {
    Measurement_Orders(std::string name, std::vector<double> orders,
                       int precision = 8)
        : Measurement{name, precision}, orders{orders} {}

    std::vector<double> orders;

    std::ostream &info_runtime(std::ostream &output) const {
      Measurement::info_runtime(output);
      output << "    - Orders: ";
      io::print(output, orders, false, " ") << "\n";
      return output;
    }
  };

  std::vector<std::unique_ptr<Measurement>> measurements;

  /**
     \brief Constructor.
     \param directories Current case directory information.
     \param parameter_set_name Name of parameter set.
     \param geometry Domain geometry info and utilities.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
  */
  template <typename Geometry, typename TransportParameters,
            typename ReactionParameters, typename SolverParameters>
  OutputParameters_Cases(Directories const &directories,
                         std::string const &parameter_set_name,
                         Geometry const &geometry,
                         TransportParameters const &params_transport,
                         ReactionParameters const &params_reaction,
                         SolverParameters const &params_solvers) {
    std::string filename =
        directories.dir_parameters + "/output_" + parameter_set_name + ".param";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    std::size_t param_index = 0;
    io::read(split_line, param_index, in_file + "Could not parse time units",
             time_units);
    std::string for_time_units =
        std::string{"For time units "} + time_units + " : ";
    if (!TimeUnitsList::contains(time_units)) {
      throw std::runtime_error{in_file + for_time_units + "Not supported"};
    }
    if (TimeUnitsList::type(time_units) == TimeUnitsList::Type::advection &&
        !params_solvers.advection) {
      throw std::runtime_error{in_file + for_time_units + "No advection"};
    }
    if (TimeUnitsList::type(time_units) == TimeUnitsList::Type::diffusion &&
        !params_solvers.diffusion) {
      throw std::runtime_error{in_file + for_time_units + "No diffusion"};
    }
    if (TimeUnitsList::type(time_units) == TimeUnitsList::Type::reaction &&
        !params_solvers.reaction) {
      throw std::runtime_error{in_file + for_time_units + "No reaction"};
    }
    time_unit_factor =
        ptof::time_unit_factor(time_units, params_transport, params_reaction);
    read_end_criteria(input, filename);
    read_measurement_spacing(input, filename, params_transport,
                             params_reaction);
    read_measurement_types(input, filename, geometry);
  }

  /** \brief Output generic information about object. */
  static std::ostream &info(std::ostream &output) {
    output << io::line() << "Output parameters\n"
           << io::line() << R"(- Time units for measurement times:
  - advection
    - Advection time units (valid only if there is advection).
  - diffusion
    - Diffusion time units (valid only if there is diffusion).
  - reaction
    - Reaction time units (valid only if there is reaction).
  - arbitrary
    - Arbitary units (no rescaling).
- End criterion to finish dynamics:
  (Note:
    - To combine criteria pass 'and' or 'or', followed by the number of criteria
      to combine on the same line, then specify one criterion per subsequent
      line.)
  - time
    - Specified time.
    - Pass on same line:
      - End time.
  - time_max
    - Maximum output time.
  - mass_below
    - Total mass below or equal to value.
    - Pass on same line:
      - Threshold mass.
  - mass_above
    - Total mass above or equal to value.
    - Pass on same line:
      - Threshold mass.
  - all_absorbed
    - All particles absorbed.
  - one_absorbed
    - One particle absorbed.
  - fraction_not_absorbed
    - Fraction of particles not absorbed below or equal to value.
    - Pass on same line:
      - Threshold fraction.
- Measurement spacing:
  - linear
    - Linear spacing, specified maximum time and number of measurements.
    - Pass on same line:
      - Minimum measurement time.
      - Maximum measurement time.
      - Number of measurements.
  - linear_step
    - Linear spacing, specified time step (no maximum time).
    - Pass on same line:
      - Minimum measurement time.
      - Measurement time step.
  - log
    - Log spacing, specified maximum time and number of measurements
    - Pass on same line:
      - Minimum measurement time.
      - Maximum measurement time.
      - Number of measurements.
  - log_step
    - Log spacing, specified time step factor (no maximum time).
    - Pass on same line:
      - Minimum measurement time.
      - Measurement time step factor.
- Measurement types:
  (Note:
    - Pass any number, one per line.
    - Precision [8] can be passed at end of line.
    - Masks are passed on construction of output handler; the same masks are
      used for all measurements involving masks.)
  - position
    - Particle tags, masses, cells, and positions, as well as times and number
      of particles output per time.
  - position_in_regions
    - Particle tags, masses, cells, and positions within regions specified by
      masks, as well as times and number of particles output per time,
  - position_mean
    - Time and mean position.
  - position_abs_mean
    - Time and mean of absolute value of position.
  - position_second_moment
    - Time and position second moment.
  - position_nth_moment
    - Time and position nth moment.
  - position_moment
    - Time and position moment with orders specified along each dimension.
  - position_variance
    - Time and position variance.
  - mass
    - Time and total mass.
  - mass_absorbed
    - Time and total absorbed mass.
  - mass_adsorbed
    - Time and total adsorbed mass.
  - mass_in_regions
    - Time and total mass within regions specified by masks.
  - velocity
    - Particle tags, masses, and velocities, as well as times and number of
      particles output per time to a separate file.
  - velocity_mean
    - Time and mean of velocity field over particles.
  - scalar_field
    - Particle tags, masses, and scalar field values, as well as times and
      number of particles output per time to a separate file.
    - Pass on same line:
      - Name of field to be read from OF file.
  - scalar_field_mean
    - Time and mean of scalar field over particles.
    - Pass on same line:
      - Name of field to be read from OF file.
  - vector_field
    - Particle tags, masses, and vector field values, as well as times and
      number of particles output per time to a separate file.
    - Pass on same line:
      - Name of field to be read from OF file.
  - vector_field_mean
    - Time and mean of vector field over particles.
    - Pass on same line:
      - Name of field to be read from OF file.
  - tensor_field
    - Particle tags, masses, and tensor field values, as well as times and
      number of particles output per time to a separate file.
    - Pass on same line:
      - Name of field to be read from OF file.
  - tensor_field_mean
    - Time and mean of tensor field over particles.
    - Pass on same line:
      - Name of field to be read from OF file.
  - position_periodic
    - Particle tags, masses, cells, and true positions accounting for
      periodicity, as well as times and number of particles output per time to a
      separate file.
  - position_in_regions_periodic
    - Particle tags, masses, cells, and true positions accounting for
      periodicity in regions specified by masks.
  - position_mean_periodic
    - Time and true mean position accounting for periodicity.
  - position_abs_mean_periodic
    - Time and true mean of absolute value of position accounting for
      periodicity.
  - position_second_moment_periodic
    - Time and true position second moment accounting for periodicity.
  - position_nth_moment_periodic
    - Time and true position nth moment accounting periodicity.
  - position_moment_periodic
    - Time and true position moment with orders specified along each dimension
      accounting for periodicity.
  - position_variance_periodic
    - Time and true position variance accounting for periodicity.
  - first_crossing_time
    - Particle crossing times, tags, and masses at end of dynamics.
    - Pass on same line:
      - Dimension.
      - Crossing position along dimension.
  - position_adsorbed
    - Time, tag, position, and mass of adsorbed particles.
  - position_adsorbed_periodic
    - Time, tag, position accounting for periodicity, and mass of adsorbed
      particles.
  - mass_adsorbed_face
    - Boundary faces and corresponding adsorbed masses, as well as times and
      number of faces output per time to a separate file.
  - mass_adsorbed_face_periodic
    - Boundary faces and corresponding adsorbed masses accouting for
      periodicity, as well as times and number of faces output per time to a
      separate file.
  - mass_reacted_face
    - Boundary faces and corresponding net reacted masses, as well as times and
      number of faces output per time to a separate file.
  - mass_reacted_face_periodic
    - Boundary faces and corresponding net reacted masses accouting for
      periodicity, as well as times and number of faces output per time to a
      separate file.
  - absorption_time
    - Particle absorption times, tags, and masses at end of dynamics.
  - absorption_time_patch
    (Note: Patch info is not guaranteed accurate unless particle states store
           boundary face.)
    - Particle absorption times, tags, particle masses, and absorption patch at
      end of dynamics.
  - absorption_time_position
    - Particle absorption times, tags, particle masses, and positions at end of
      dynamics.
  - absorption_time_patch_position
    (Note: Patch info is not guaranteed accurate unless particle
           states store boundary face.)
    - Particle absorption times, tags, masses, absorption patch, and positions
      at end of dynamics.
  - absorption_time_position_periodic
    - Particle absorption times, tags, masses, and positions accounting for
      periodicity at end of dynamics.
  - absorption_time_patch_position_periodic
    (Note: Patch info is not guaranteed accurate unless particle states store
           boundary face.)
    - Particle absorption times, tags, masses, absorption patch, and positions
      accounting for periodicity at end of dynamics.
)" << io::line();
    return output;
  }

  /** \brief Read all end criteria from input stream. */
  void read_end_criteria(std::ifstream &input, std::string const &filename) {
    std::string in_file = std::string{"In file "} + filename + " : ";
    auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    std::size_t param_index = 0;
    auto criterion = io::read<std::string>(
        split_line, param_index, in_file + "Could not parse end criterion");

    if (criterion == "and" || criterion == "or") {
      end_criterion_logical_combination = criterion;
      std::string for_end_criterion =
          std::string{"End criterion logical combination "} + criterion + " : ";
      auto nr_criteria = io::read<std::size_t>(
          split_line, param_index,
          in_file + for_end_criterion +
              "Could not parse number of end criteria to combine");
      for (std::size_t ii = 0; ii < nr_criteria; ++ii) {
        auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
        std::size_t param_index = 0;
        io::read(split_line, param_index, 1,
                 in_file + "Could not parse end criterion", end_criteria);
        parse_end_criterion_parameters(split_line, param_index,
                                       end_criteria.back(), in_file);
      }
    } else {
      end_criterion_logical_combination = "single";
      end_criteria.push_back(criterion);
      parse_end_criterion_parameters(split_line, param_index,
                                     end_criteria.back(), in_file);
    }

    if (end_criteria.empty()) {
      throw std::runtime_error{"No end criterion"};
    }

    if (std::unordered_set<std::string>(end_criteria.begin(),
                                        end_criteria.end())
            .size() != end_criteria.size()) {
      throw std::runtime_error{
          "The same end criterion type was specified more than once"};
    }
  }

  /** \brief Parse one end criterion's parameters from split line string. */
  void parse_end_criterion_parameters(
      std::vector<std::string> const &split_line, std::size_t &param_index,
      std::string const &end_criterion, std::string const &in_file = "") {
    std::string for_end_criterion =
        std::string{"End criterion "} + end_criterion + " : ";
    if (!EndCriterionList::contains(end_criterion)) {
      throw std::runtime_error{for_end_criterion + end_criterion +
                               " Not supported"};
    }
    auto criterion_type = EndCriterionList::type(end_criterion);
    switch (criterion_type) {
    case EndCriterionList::Type::time: {
      end_values[criterion_type] =
          io::read<double>(split_line, param_index,
                           in_file + for_end_criterion +
                               "Could not parse end time") *
          time_unit_factor;
      break;
    }
    case EndCriterionList::Type::time_max:
      break;
    case EndCriterionList::Type::mass_below:
    case EndCriterionList::Type::mass_above: {
      end_values[criterion_type] = io::read<double>(
          split_line, param_index,
          in_file + for_end_criterion + "Could not parse mass threshold");
      break;
    }
    case EndCriterionList::Type::all_absorbed:
      break;
    case EndCriterionList::Type::one_absorbed:
      break;
    case EndCriterionList::Type::fraction_not_absorbed: {
      end_values[criterion_type] =
          io::read<double>(split_line, param_index,
                           in_file + for_end_criterion +
                               "Could not parse absorbed fraction threshold");
      break;
    }
    default: {
      throw std::runtime_error{in_file + for_end_criterion + "Not supported"};
    }
    }
  }

  /** \brief Read measurement spacing type from input stream. */
  template <typename TransportParameters, typename ReactionParameters>
  void read_measurement_spacing(std::ifstream &input,
                                std::string const &filename,
                                TransportParameters const &params_transport,
                                ReactionParameters const &params_reaction) {
    std::string in_file = std::string{"In file "} + filename + " : ";
    auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse measurement spacing",
             measurement_spacing);
    std::string for_measurement_spacing =
        std::string{"Measurement spacing "} + measurement_spacing + " : ";
    if (!MeasurementSpacingList::contains(measurement_spacing)) {
      throw std::runtime_error{in_file + for_measurement_spacing +
                               "Not supported"};
    }
    io::read(split_line, param_index,
             in_file + for_measurement_spacing +
                 "Could not parse minimum measurement time",
             time_min);
    time_min *= time_unit_factor;
    switch (MeasurementSpacingList::type(measurement_spacing)) {
    case MeasurementSpacingList::Type::linear_step: {
      io::read(split_line, param_index, in_file + "Could not parse time step",
               time_increment);
      time_increment *= time_unit_factor;
      time_max = std::numeric_limits<double>::infinity();
      if (!(time_increment > 0.))
        throw std::runtime_error{in_file + for_measurement_spacing +
                                 "Non-positive time increment"};
      break;
    }
    case MeasurementSpacingList::Type::log_step: {
      io::read(split_line, param_index,
               in_file + for_measurement_spacing +
                   "Could not parse time step factor",
               time_increment);
      time_max = std::numeric_limits<double>::infinity();
      if (!(time_increment > 1.))
        throw std::runtime_error{in_file + for_measurement_spacing +
                                 "Time increment factor not greater than one"};
      break;
    }
    case MeasurementSpacingList::Type::linear: {
      std::size_t nr_measurements;
      io::read(split_line, param_index,
               in_file + for_measurement_spacing +
                   "Could not parse maximum time and number of measurements",
               time_max, nr_measurements);
      time_max *= time_unit_factor;
      if (nr_measurements < 2)
        throw std::runtime_error{in_file + for_measurement_spacing +
                                 "Number of measurements not greater than one"};
      time_increment = (time_max - time_min) / (nr_measurements - 1);
      if (!(time_increment > 0.))
        throw std::runtime_error{in_file + for_measurement_spacing +
                                 "Non-positive time increment"};
      break;
    }
    case MeasurementSpacingList::Type::log: {
      if (!(time_min > 0.))
        throw std::runtime_error{in_file + for_measurement_spacing +
                                 "Non-positive minimum measurement time"};
      std::size_t nr_measurements;
      io::read(split_line, param_index,
               in_file + for_measurement_spacing +
                   "Could not parse maximum time and number of measurements",
               time_max, nr_measurements);
      time_max *= time_unit_factor;
      if (nr_measurements < 2)
        throw std::runtime_error{in_file + for_measurement_spacing +
                                 "Number of measurements not greater than one"};
      time_increment =
          std::pow(time_max / time_min, 1. / (nr_measurements - 1));
      if (!(time_increment > 1.))
        throw std::runtime_error{in_file + for_measurement_spacing +
                                 "Time increment factor not greater than one"};
      break;
    }
    default: {
      throw std::runtime_error{in_file + for_measurement_spacing +
                               "Not supported"};
    }
    }
  }

  /** \brief Read measurement type from input stream. */
  template <typename Geometry>
  void read_measurement_types(std::ifstream &input, std::string const &filename,
                              Geometry const &geometry) {
    std::string in_file = std::string{"In file "} + filename + " : ";
    for (std::vector<std::string> split_line;
         io::split_line(input, split_line, "#", "\t,|\r()[]{} ");
         split_line.clear()) {
      std::size_t param_index = 0;
      auto name =
          io::read<std::string>(split_line, param_index,
                                in_file + "Could not parse measurement type");
      std::string for_measurement_type =
          std::string{"Measurement type "} + name + " : ";
      if (!MeasurementList::contains(name)) {
        throw std::runtime_error{in_file + for_measurement_type +
                                 "Not supported"};
      }
      switch (MeasurementList::type(name)) {
      case (MeasurementList::Type::scalar_field):
      case (MeasurementList::Type::vector_field):
      case (MeasurementList::Type::tensor_field): {
        auto field_name = io::read<std::string>(
            split_line, param_index,
            in_file + for_measurement_type + "Could not parse field name");
        measurements.emplace_back(
            std::make_unique<Measurement_Field>(name, field_name));
        break;
      }
      case (MeasurementList::Type::first_crossing_time): {
        auto dim = io::read<double>(split_line, param_index,
                                    in_file + for_measurement_type +
                                        "Could not parse crossing dimension");
        auto position =
            io::read<double>(split_line, param_index,
                             in_file + for_measurement_type +
                                 "Could not parse crossing position");
        if (dim >= Geometry::dim) {
          throw std::runtime_error{
              in_file + for_measurement_type +
              "Crossing dimension higher than simulation dimension"};
        }
        measurements.emplace_back(
            std::make_unique<Measurement_Dim_Position>(name, dim, position));
        break;
      }
      case (MeasurementList::Type::position_nth_moment):
      case (MeasurementList::Type::position_nth_moment_periodic): {
        auto order = io::read<double>(split_line, param_index,
                                      in_file + for_measurement_type +
                                          "Could not parse moment order");
        measurements.emplace_back(
            std::make_unique<Measurement_Order>(name, order));
        break;
      }
      case (MeasurementList::Type::position_moment):
      case (MeasurementList::Type::position_moment_periodic): {
        auto orders = io::read<double>(split_line, param_index, Geometry::dim,
                                       in_file + for_measurement_type +
                                           "Could not parse moment coordinate "
                                           "orders along each dimension");
        measurements.emplace_back(
            std::make_unique<Measurement_Orders>(name, orders));
        break;
      }
      default:
        measurements.emplace_back(std::make_unique<Measurement>(name));
      }

      if (param_index < split_line.size()) {
        io::read(split_line, param_index,
                 in_file + for_measurement_type + "Could not parse precision",
                 measurements.back()->precision);
      }
    }
  }
};
} // namespace ptof

#endif /* PTOF_OUTPUTPARAMETERS_H */
