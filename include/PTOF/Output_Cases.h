/**
   \file PTOF/Output_Cases.h
   \author Tomás Aquino
   \date 07/03/2022
   \brief Handle different types of measurements and output.
*/

#ifndef PTOF_OUTPUT_CASES_H
#define PTOF_OUTPUT_CASES_H

#include "CTRW/Meta.h"
#include "General/IO.h"
#include "General/Meta.h"
#include "PTOF/BoundaryConditionList.h"
#include "PTOF/Criterion.h"
#include "PTOF/Directories.h"
#include "PTOF/Field.h"
#include "PTOF/InitialCondition_Cases.h"
#include "PTOF/MeasurementList.h"
#include "PTOF/MeasurementSpacing.h"
#include "PTOF/Measurer.h"
#include "PTOF/MeasurerTime.h"
#include "PTOF/NextMeasurement.h"
#include "PTOF/TimeUnits.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <fvcGrad.H>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeinfo>
#include <vector>

namespace ptof {
/**
   \struct OutputParameters_Cases_Cases PTOF/Output_Cases.h
   "PTOF/Output_Cases.h"
   \brief Parameters for output.
*/
struct OutputParameters_Cases {
public:
  std::string time_units;
  double time_unit_factor;
  std::string end_criterion;
  double end_value;
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

    inline void info_runtime(std::ostream &output) const {
      output << name << "\n";
      output << "  - Precision: " << precision << "\n";
    }
  };

  struct Measurement_Field : public Measurement {
    Measurement_Field(std::string name, std::string field_name,
                      int precision = 8)
        : Measurement{name, precision}, field_name{field_name} {}

    std::string field_name;

    inline void info_runtime(std::ostream &output) const {
      Measurement::info_runtime(output);
      output << "  - Field: " << field_name << "\n";
    }
  };

  struct Measurement_Dim_Position : public Measurement {
    Measurement_Dim_Position(std::string name, std::size_t dim, double position,
                             int precision = 8)
        : Measurement{name, precision}, dim{dim}, position{position} {}

    std::size_t dim;
    double position;

    inline void info_runtime(std::ostream &output) const {
      Measurement::info_runtime(output);
      output << "  - Dimension: " << dim << "\n"
             << "  - Position: " << position << "\n";
    }
  };

  struct Measurement_Order : public Measurement {
    Measurement_Order(std::string name, double order, int precision = 8)
        : Measurement{name, precision}, order{order} {}

    double order;

    inline void info_runtime(std::ostream &output) const {
      Measurement::info_runtime(output);
      output << "  - Order: " << order << "\n";
    }
  };

  struct Measurement_Orders : public Measurement {
    Measurement_Orders(std::string name, std::vector<double> orders,
                       int precision = 8)
        : Measurement{name, precision}, orders{orders} {}

    std::vector<double> orders;

    inline void info_runtime(std::ostream &output) const {
      Measurement::info_runtime(output);
      output << "  - Orders: ";
      io::print(output, orders, false, " ") << "\n";
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
    std::string filename = directories.dir_parameters + "/parameters_output_" +
                           parameter_set_name + ".dat";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    io::read(split_line, param_index, in_file + "Could not parse time units",
             time_units);
    if (!TimeUnits{}.contains(time_units))
      throw std::runtime_error{in_file + "Not supported"};
    time_unit_factor =
        ptof::time_unit_factor(time_units, params_transport, params_reaction);
    read_end_criterion(input, filename);
    read_measurement_spacing(input, filename, params_transport,
                             params_reaction);
    read_measurement_types(input, filename, geometry);
  }

  /** \brief Output generic information about object. */
  inline static void info(std::ostream &output) {
    output
        << "--------------------------------------------------------------"
           "\n"
           "Output parameters\n"
           "--------------------------------------------------------------"
           "\n"
           "- Time units for measurement times:\n"
           "  - diffusion\n"
           "    - Diffusion time units\n"
           "  - advection\n"
           "    - Advection time units\n"
           "  - reaction\n"
           "    - Reaction time units\n"
           "  - arbitrary\n"
           "    - Arbitary units (no rescaling)\n"
           "- End criterion to finish dynamics:\n"
           "  - time\n"
           "    - Specified time\n"
           "    - Pass on same line:\n"
           "      - End time\n"
           "  - time_max\n"
           "    - Maximum output time\n"
           "  - mass_below\n"
           "    - Total mass below or equal to value\n"
           "    - Pass on same line:\n"
           "      - Threshold mass\n"
           "  - mass_above\n"
           "    - Total mass above or equal to value\n"
           "    - Pass on same line:\n"
           "      - Threshold mass\n"
           "  - all_absorbed\n"
           "    - All particles absorbed\n"
           "  - one_absorbed\n"
           "    - One particle absorbed\n"
           "  - fraction_not_absorbed\n"
           "    - Fraction of particles not absorbed below or equal to value\n"
           "    - Pass on same line:\n"
           "      - Threshold fraction\n"
           "- Measurement spacing:\n"
           "  - linear\n"
           "    - Linear spacing, specified maximum time and number of\n"
           "      measurements\n"
           "    - Pass on same line:\n"
           "      - Minimum measurement time\n"
           "      - Maximum measurement time\n"
           "      - Number of measurements\n"
           "  - linear_step\n"
           "    - Linear spacing, specified time step (no maximum time)\n"
           "    - Pass on same line:\n"
           "      - Minimum measurement time\n"
           "      - Measurement time step\n"
           "  - log\n"
           "    - Log spacing, specified maximum time and number of\n"
           "      measurements\n"
           "    - Pass on same line:\n"
           "      - Minimum measurement time\n"
           "      - Maximum measurement time\n"
           "      - Number of measurements\n"
           "  - log_step\n"
           "    - Log spacing, specified time step factor (no maximum time)\n"
           "    - Pass on same line:\n"
           "      - Minimum measurement time\n"
           "      - Measurement time step factor\n"
           "- Measurement types:\n"
           "  (Note:\n"
           "    - Pass any number, one per line\n"
           "    - Precision [8] can be passed at end of line\n"
           "    - Masks are passed on construction of output handler; the\n"
           "      same masks are used for all measurements involving masks\n"
           "  - position\n"
           "    - Time, particle tags, particle positions, and particle\n"
           "      masses\n"
           "  - position_in_regions\n"
           "    - Time, particle tags, particle positions, and particle\n"
           "      masses within regions specified by masks\n"
           "  - position_mean\n"
           "    - Time and mean position\n"
           "  - position_second_moment\n"
           "    - Time and position second moment\n"
           "  - position_nth_moment\n"
           "    - Time and position nth moment\n"
           "  - position_moment\n"
           "    - Time and position moment with orders specified along each\n"
           "      dimension\n"
           "  - position_variance\n"
           "    - Time and position variance\n"
           "  - mass\n"
           "    - Time and total mass\n"
           "  - mass_in_regions\n"
           "    - Time and total mass within regions specified by masks\n"
           "  - velocity\n"
           "    - Time, particle tags, local velocities, and particle masses\n"
           "  - velocity_mean\n"
           "    - Time and mean of velocity field over particles\n"
           "  - velocity_gradient\n"
           "    - Time, particle tags, local velocity gradient, and particle\n"
           "      masses\n"
           "  - velocity_gradient_mean\n"
           "    - Time and mean of velocity gradient field over particles\n"
           "  - scalar_field\n"
           "    - Time, particle tags, local scalar field values, and\n"
           "      particle masses\n"
           "    - Pass on same line:\n"
           "      - Name of field to be read from OF file\n"
           "  - scalar_field_mean\n"
           "    - Time and mean of scalar field over particles\n"
           "    - Pass on same line:\n"
           "      - Name of field to be read from OF file\n"
           "  - vector_field\n"
           "    - Time, particle tags, local vector field value, and particle\n"
           "      masses\n"
           "    - Pass on same line:\n"
           "      - Name of field to be read from OF file\n"
           "  - vector_field_mean\n"
           "    - Time and mean of vector field over particles\n"
           "    - Pass on same line:\n"
           "      - Name of field to be read from OF file\n"
           "  - tensor_field\n"
           "    - Time, particle tags, local tensor field values, and\n"
           "      particle masses\n"
           "    - Pass on same line:\n"
           "      - Name of field to be read from OF file\n"
           "  - tensor_field_mean\n"
           "    - Time and mean of tensor field over particles\n"
           "    - Pass on same line:\n"
           "      - Name of field to be read from OF file\n"
           "  - position_periodic\n"
           "    - Time, particle tags, true positions accounting for\n"
           "      periodicity, and masses\n"
           "  - position_in_regions_periodic\n"
           "    - Time, particle tags, true positions accounting for\n"
           "      periodicity, and masses, in regions specified by masks\n"
           "  - position_mean_periodic\n"
           "    - Time and true mean position accounting for periodicity\n"
           "  - position_second_moment_periodic\n"
           "    - Time and true position second moment accounting for\n"
           "      periodicity\n"
           "  - position_nth_moment_periodic\n"
           "    - Time and true position nth moment accounting for\n"
           "      periodicity\n"
           "  - position_moment_periodic\n"
           "    - Time and true position moment with orders specified along\n"
           "      each dimension accountin for periodicity\n"
           "  - position_variance_periodic\n"
           "    - Time and true position variance accounting for periodicity\n"
           "  - first_crossing_time\n"
           "    - Crossing times, particle tags, and particle masses at end\n"
           "      of dynamics\n"
           "    - Pass on same line:\n"
           "      - Dimension\n"
           "      - Crossing position along dimension\n"
           "  - absorption_time\n"
           "    - Particle absorption times, particle tags and particle\n"
           "      masses at end of dynamics\n"
           "--------------------------------------------------------------"
           "\n";
  }

  /** \brief Read end criterion from input stream. */
  template <typename IStream>
  void read_end_criterion(IStream &input, std::string const &filename) {
    std::string in_file = std::string{"In file "} + filename + " : ";
    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    io::read(split_line, param_index, in_file + "Could not parse end criterion",
             end_criterion);
    std::string for_end_criterion =
        std::string{"End criterion "} + end_criterion + " : ";

    switch (EndCriterion::type(end_criterion)) {
    case EndCriterion::Type::time: {
      io::read(split_line, param_index,
               in_file + for_end_criterion + "Could not parse end time",
               end_value);
      end_value *= time_unit_factor;
      break;
    }
    case EndCriterion::Type::time_max:
      break;
    case EndCriterion::Type::mass_below:
    case EndCriterion::Type::mass_above: {
      io::read(split_line, param_index,
               in_file + for_end_criterion + "Could not parse mass threshold",
               end_value);
      break;
    }
    case EndCriterion::Type::all_absorbed:
      break;
    case EndCriterion::Type::one_absorbed:
      break;
    case EndCriterion::Type::fraction_not_absorbed: {
      io::read(split_line, param_index,
               in_file + for_end_criterion +
                   "Could not parse absorbed fraction threshold",
               end_value);
      break;
    }
    default:
      throw std::runtime_error{std::string{"End criterion "} + end_criterion +
                               " not supported"};
    }
  }

  /** \brief Read measurement spacing type from input stream. */
  template <typename IStream, typename TransportParameters,
            typename ReactionParameters>
  void read_measurement_spacing(IStream &input, std::string const &filename,
                                TransportParameters const &params_transport,
                                ReactionParameters const &params_reaction) {
    std::string in_file = std::string{"In file "} + filename + " : ";
    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse measurement spacing",
             measurement_spacing);
    std::string for_measurement_spacing =
        std::string{"Measurement spacing "} + measurement_spacing + " : ";

    io::read(split_line, param_index,
             in_file + for_measurement_spacing +
                 "Could not parse minimum measurement time",
             time_min);
    time_min *= time_unit_factor;
    switch (MeasurementSpacing::type(measurement_spacing)) {
    case MeasurementSpacing::Type::linear_step: {
      io::read(split_line, param_index, in_file + "Could not parse time step",
               time_increment);
      time_increment *= time_unit_factor;
      time_max = std::numeric_limits<double>::infinity();
      if (!(time_increment > 0.))
        throw std::runtime_error{in_file + for_measurement_spacing +
                                 "Non-positive time increment"};
      break;
    }
    case MeasurementSpacing::Type::log_step: {
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
    case MeasurementSpacing::Type::linear: {
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
    case MeasurementSpacing::Type::log: {
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
    default:
      throw std::runtime_error{in_file + for_measurement_spacing +
                               "Not supported"};
    }
  }

  /** \brief Read measurement type from input stream. */
  template <typename IStream, typename Geometry>
  void read_measurement_types(IStream &input, std::string const &filename,
                              Geometry const &geometry) {
    std::string in_file = std::string{"In file "} + filename + " : ";
    for (std::vector<std::string> split_line; io::split_line(input, split_line);
         split_line.clear()) {
      std::size_t param_index = 0;
      auto name =
          io::read<std::string>(split_line, param_index,
                                in_file + "Could not parse measurement type");
      std::string for_measurement_type =
          std::string{"Measurement type "} + name + " : ";
      if (!MeasurementList{}.contains(name))
        throw std::runtime_error{"Output type " + name + " not supported"};

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
        if (dim >= Geometry::dim)
          throw std::runtime_error{
              in_file + for_measurement_type +
              "Crossing dimension higher than simulation dimension"};
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

      if (param_index < split_line.size())
        io::read(split_line, param_index,
                 in_file + for_measurement_type + "Could not parse precision",
                 measurements.back()->precision);
    }
  }
};

/**
   \class Output_Cases PTOF/Output_Cases.h "PTOF/Output_Cases.h"
   \brief Output object to handle implemented output options.
*/
template <typename Subject, typename Geometry, typename Parameters>
class Output_Cases {
public:
  Parameters parameters; /**< Output parameters .*/

  /**
     \brief Constructor.
     \param subject CTRW object to measure.
     \param velocity_field Velocity field as a function of state.
     \param geometry Domain geometry info and utilities.
     \param directories Current case directory information.
     \param params_output Output parameters.
     \param identifier String to include in names of output files.
     \param masks Scalar fields reference wrappers.
     \param thresholds Thresholds for each mask, such that cells where a mask is
     above or equal to the threshold are considered.
  */
  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  Output_Cases(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters &&params_output, std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks = {},
               std::vector<double> thresholds = {})
      : parameters{std::forward<Parameters>(params_output)} {
    set_measurement_types(subject, geometry, directories, identifier,
                          velocity_field, masks, thresholds);
    set_end_criterion(subject);
    set_next_measurement_time();
  }

  /**
     \brief Constructor.
     \details Overload for initializer list arguments.
  */
  template <typename VelocityField = meta::Empty, typename Mask>
  Output_Cases(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters &&params_output, std::string const &identifier,
               std::initializer_list<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {})
      : Output_Cases(subject, velocity_field, geometry, directories,
                     std::forward<Parameters>(params_output), identifier,
                     std::vector<std::reference_wrapper<const Mask>>{masks},
                     std::vector<double>{thresholds}) {}

  /**
     \brief Constructor.
     \details Overload for initializer list arguments.
  */
  template <typename VelocityField = meta::Empty, typename Mask>
  Output_Cases(Subject const &subject, VelocityField &&velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters &&params_output, std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {})
      : Output_Cases(subject, velocity_field, geometry, directories,
                     std::forward<Parameters>(params_output), identifier, masks,
                     std::vector<double>{thresholds}) {}

  /**
     \brief Constructor.
     \details Overload for initializer list arguments.
  */
  template <typename VelocityField = meta::Empty, typename Mask>
  Output_Cases(Subject const &subject, VelocityField &&velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters &&params_output, std::string const &identifier,
               std::initializer_list<std::reference_wrapper<const Mask>> masks,
               std::vector<double> thresholds = {})
      : Output_Cases(subject, velocity_field, geometry, directories,
                     std::forward<Parameters>(params_output), identifier, masks,
                     std::vector<double>{thresholds}) {}

  /**
     \return \c true  if end simulation criterion is satisfied, \c false
     otherwise.
  */
  bool done(double time) const { return _end_criterion->operator()(time); }

  /** \brief Set time of next measurement. */
  void set_next_measurement_time() {
    switch (MeasurementSpacing::type(parameters.measurement_spacing)) {
    case MeasurementSpacing::Type::linear_step: {
      _next_measurement = std::make_unique<
          NextMeasurementTime_linear_step<std::remove_reference_t<Parameters>>>(
          parameters);
      break;
    }
    case MeasurementSpacing::Type::log_step: {
      _next_measurement = std::make_unique<
          NextMeasurementTime_log_step<std::remove_reference_t<Parameters>>>(
          parameters);
      break;
    }
    case MeasurementSpacing::Type::linear: {
      _next_measurement = std::make_unique<
          NextMeasurementTime_linear<std::remove_reference_t<Parameters>>>(
          parameters);
      break;
    }
    case MeasurementSpacing::Type::log: {
      _next_measurement = std::make_unique<
          NextMeasurementTime_log<std::remove_reference_t<Parameters>>>(
          parameters);
      break;
    }
    default:
      throw std::runtime_error{"Measurement spacing " +
                               parameters.measurement_spacing +
                               " not supported"};
    }
  }

  /** \return Time of next measurement. */
  double next_measurement_time() const { return _next_measurement->time(); }

  /**
     \brief Output requested measurements at given time and advance to next
     measurement.
  */
  void operator()(double time) {
    for (auto const &output : _output_time)
      output->operator()(time);
    _next_measurement->advance();
  }

  /** \brief Output current information. */
  void operator()() {
    for (auto const &output : _output)
      output->print();
    for (auto const &output : _output_time)
      output->print();
  }

  /** \brief Update internal state. */
  void update(double time, double time_of_change) {
    for (auto const &output : _output_time)
      output->update(time, time_of_change);
  }

  /** \brief Output information about current object. */
  inline void info_runtime(std::ostream &output) const {
    output << "------------------------------------------------------------"
              "--\n"
              "Output\n"
              "------------------------------------------------------------"
              "--\n"
              "- End criterion: "
           << parameters.end_criterion << "\n";
    if (EndCriterion::type(parameters.end_criterion) ==
            EndCriterion::Type::time ||
        EndCriterion::type(parameters.end_criterion) ==
            EndCriterion::Type::mass_below ||
        EndCriterion::type(parameters.end_criterion) ==
            EndCriterion::Type::mass_above ||
        EndCriterion::type(parameters.end_criterion) ==
            EndCriterion::Type::fraction_not_absorbed)
      output << "- End value: " << parameters.end_value << "\n"
             << "- Time units: " << parameters.time_units
             << "\n"
                "- Measurement spacing: "
             << parameters.measurement_spacing
             << "\n"
                "- Minimum measurement time: "
             << parameters.time_min
             << "\n"
                "- Maximum measurement time: ";
    if (EndCriterion::type(parameters.end_criterion) ==
            EndCriterion::Type::time &&
        (MeasurementSpacing::type(parameters.measurement_spacing) ==
             MeasurementSpacing::Type::linear_step ||
         MeasurementSpacing::type(parameters.measurement_spacing) ==
             MeasurementSpacing::Type::log_step))
      output << parameters.end_value << "\n";
    else
      output << parameters.time_max << "\n";
    output << "- Measurement types:\n";
    info_runtime_measurements(output);
    output << "------------------------------------------------------------"
              "--\n";
  }

private:
  /** \brief Output information about measurement types. */
  void info_runtime_measurements(std::ostream &output) const {
    if (parameters.measurements.empty()) {
      output << "\tNone\n";
      return;
    }
    for (auto const &measurement : parameters.measurements) {
      output << "  - ";
      measurement->info_runtime(output);
      output << "\n";
    }
  }

  /**
     \brief Set up output streams for requested output types.
     \param subject CTRW object to measure.
     \param directories Current case directory information.
     \param geometry Domain geometry info and utilities.
     \param identifier String to include in names of output files.
     \param velocity_field Velocity field as a function of state.
     \param masks Scalar field reference wrappers.
     \param thresholds Thresholds for each mask, such that cells where a mask is
     above or equal to the threshold are considered.
  */
  template <typename VelocityField = meta::Empty, typename Mask = meta::Empty>
  void set_measurement_types(
      Subject const &subject, Geometry const &geometry,
      Directories const &directories, std::string const &identifier,
      VelocityField &&velocity_field = {},
      std::vector<std::reference_wrapper<const Mask>> masks = {},
      std::vector<double> thresholds = {}) {
    for (auto const &measurement : parameters.measurements) {
      switch (MeasurementList::type(measurement->name)) {
      case MeasurementList::Type::position: {
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_position<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::position_in_regions: {
        if constexpr (!std::is_same_v<Mask, meta::Empty>)
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_in_regions<Subject, Geometry, Mask>>(
                  subject, geometry, directories, identifier, masks, thresholds,
                  measurement->precision));
        else
          throw std::runtime_error{std::string{"Measurement type "} +
                                   measurement->name +
                                   " : "
                                   "Region masks not provided"};
        break;
      }
      case MeasurementList::Type::position_mean: {
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_position_mean<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::position_second_moment: {
        _output_time.emplace_back(
            std::make_unique<
                MeasurerTime_position_second_moment<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::position_nth_moment: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Order;
        Measurement *measurement_derived = cast<Measurement>(measurement);
        if (measurement_derived->order == 0)
          throw std::runtime_error{std::string{"Measurement type "} +
                                   measurement->name +
                                   " : "
                                   "Moment order not provided"};
        _output_time.emplace_back(
            std::make_unique<
                MeasurerTime_position_nth_moment<Subject, Geometry>>(
                subject, measurement_derived->order, geometry, directories,
                identifier, measurement->precision));
        break;
      }
      case MeasurementList::Type::position_moment: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Orders;
        Measurement *measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_position_moment<Subject, Geometry>>(
                subject, measurement_derived->orders, geometry, directories,
                identifier, measurement->precision));
        break;
      }
      case MeasurementList::Type::position_variance: {
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_position_variance<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::mass: {
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_mass<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::mass_in_regions: {
        if constexpr (!std::is_same_v<Mask, meta::Empty>)
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_mass_in_regions<Subject, Geometry, Mask>>(
                  subject, geometry, directories, identifier, masks, thresholds,
                  measurement->precision));
        else
          throw std::runtime_error{std::string{"Measurement type "} +
                                   measurement->name +
                                   " : "
                                   "Region masks not provided"};
        break;
      }
      case MeasurementList::Type::velocity: {
        if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_vector_field<
                  Subject, Geometry, VelocityField const &>>(
                  subject, std::forward<VelocityField>(velocity_field),
                  geometry, directories, identifier, "U",
                  measurement->precision));
        else
          throw std::runtime_error{std::string{"Measurement type "} +
                                   measurement->name + " : " +
                                   "Velocity field not provided"};
        break;
      }
      case MeasurementList::Type::velocity_mean: {
        if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_vector_field_mean<
                  Subject, Geometry, VelocityField const &>>(
                  subject, std::forward<VelocityField>(velocity_field),
                  geometry, directories, identifier, "U",
                  measurement->precision));
        else
          throw std::runtime_error{std::string{"Measurement type "} +
                                   measurement->name + " : " +
                                   "Velocity field not provided"};
        break;
      }
      case MeasurementList::Type::velocity_gradient: {
        if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_tensor_field<Subject, Geometry>>(
                  subject, Foam::fvc::grad(velocity_field.field()), geometry,
                  directories, identifier, "gradU", measurement->precision));
        else
          throw std::runtime_error{std::string{"Measurement type "} +
                                   measurement->name + " : " +
                                   "Velocity field not provided"};
        break;
      }
      case MeasurementList::Type::velocity_gradient_mean: {
        if constexpr (!std::is_same_v<VelocityField, meta::Empty>)
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_tensor_field_mean<Subject, Geometry>>(
                  subject, Foam::fvc::grad(velocity_field.field()), geometry,
                  directories, identifier, "gradU", measurement->precision));
        else
          throw std::runtime_error{std::string{"Measurement type "} +
                                   measurement->name + " : " +
                                   "Velocity field not provided"};
        break;
      }
      case MeasurementList::Type::scalar_field: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement *measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_scalar_field<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::vector_field: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement *measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_vector_field<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::tensor_field: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement *measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_tensor_field<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::scalar_field_mean: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement *measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_scalar_field_mean<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::vector_field_mean: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement *measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_vector_field_mean<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::tensor_field_mean: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement *measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_tensor_field_mean<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::position_periodic: {
        if constexpr (meta::has_periodicity_v<
                          typename Subject::Particle::State>)
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_periodic<Subject, Geometry>>(
                  subject, geometry, directories, identifier,
                  measurement->precision));
        else
          throw std::runtime_error{
              std::string{"Measurement type "} + measurement->name +
              " : "
              "Particle state does not define periodicity"};
        break;
      }
      case MeasurementList::Type::position_in_regions_periodic: {
        if constexpr (meta::has_periodicity_v<
                          typename Subject::Particle::State> &&
                      !std::is_same_v<Mask, meta::Empty>)
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_position_in_regions_periodic<
                  Subject, Geometry, Mask>>(subject, geometry, directories,
                                            identifier, masks, thresholds,
                                            measurement->precision));
        else {
          if constexpr (!meta::has_periodicity_v<
                            typename Subject::Particle::State>)
            throw std::runtime_error{
                std::string{"Measurement type "} + measurement->name +
                " : "
                "Particle state does not define periodicity"};
          throw std::runtime_error{std::string{"Measurement type "} +
                                   measurement->name +
                                   " : "
                                   "Region masks not provided"};
        }
        break;
      }
      case MeasurementList::Type::position_mean_periodic: {
        if constexpr (meta::has_periodicity_v<
                          typename Subject::Particle::State>)
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_mean_periodic<Subject, Geometry>>(
                  subject, geometry, directories, identifier,
                  measurement->precision));
        else
          throw std::runtime_error{
              std::string{"Measurement type "} + measurement->name +
              " : "
              "Particle state does not define periodicity"};
        break;
      }
      case MeasurementList::Type::position_second_moment_periodic: {
        if constexpr (meta::has_periodicity_v<
                          typename Subject::Particle::State>)
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_position_second_moment_periodic<
                  Subject, Geometry>>(subject, geometry, directories,
                                      identifier, measurement->precision));
        else
          throw std::runtime_error{
              std::string{"Measurement type "} + measurement->name +
              " : "
              "Particle state does not define periodicity"};
        break;
      }
      case MeasurementList::Type::position_variance_periodic: {
        if constexpr (meta::has_periodicity_v<
                          typename Subject::Particle::State>)
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_variance_periodic<Subject, Geometry>>(
                  subject, geometry, directories, identifier,
                  measurement->precision));
        else
          throw std::runtime_error{
              std::string{"Measurement type "} + measurement->name +
              " : "
              "Particle state does not define periodicity"};
        break;
      }
      case MeasurementList::Type::position_nth_moment_periodic: {
        if constexpr (meta::has_periodicity_v<
                          typename Subject::Particle::State>) {
          using Measurement =
              typename std::remove_reference_t<Parameters>::Measurement_Order;
          Measurement *measurement_derived = cast<Measurement>(measurement);
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_nth_moment_periodic<Subject, Geometry>>(
                  subject, measurement_derived->order, geometry, directories,
                  identifier, measurement->precision));
        } else
          throw std::runtime_error{
              std::string{"Measurement type "} + measurement->name +
              " : "
              "Particle state does not define periodicity"};
        break;
      }
      case MeasurementList::Type::position_moment_periodic: {
        if constexpr (meta::has_periodicity_v<
                          typename Subject::Particle::State>) {
          using Measurement =
              typename std::remove_reference_t<Parameters>::Measurement_Orders;
          Measurement *measurement_derived = cast<Measurement>(measurement);
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_moment_periodic<Subject, Geometry>>(
                  subject, measurement_derived->orders, geometry, directories,
                  identifier, measurement->precision));
        } else
          throw std::runtime_error{
              std::string{"Measurement type "} + measurement->name +
              " : "
              "Particle state does not define periodicity"};
        break;
      }
      case MeasurementList::Type::first_crossing_time: {
        using Measurement = typename std::remove_reference_t<
            Parameters>::Measurement_Dim_Position;
        Measurement *measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<
                MeasurerTime_first_crossing_time<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->dim, measurement_derived->position,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::absorption_time: {
        _output.emplace_back(
            std::make_unique<Measurer_absorption_time<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      default:
        throw std::runtime_error{std::string{"Measurement type "} +
                                 measurement->name +
                                 " "
                                 "not supported"};
      }
    }
  }

  /** \brief Set end criterion. */
  void set_end_criterion(Subject const &subject) {
    // Note: The non-standard control flow structure of the if statements is
    // needed to deal with the quirks of if constexpr
    switch (EndCriterion::type(parameters.end_criterion)) {
    case EndCriterion::Type::time: {
      _end_criterion = std::make_unique<Criterion_time<Subject>>(
          subject, parameters.end_value);
      break;
    }
    case EndCriterion::Type::time_max: {
      if (parameters.time_max == std::numeric_limits<double>::infinity()) {
        throw std::runtime_error{
            std::string{"Infinite maximum time with end criterion "} +
            parameters.end_criterion};
      }
      _end_criterion = std::make_unique<Criterion_time<Subject>>(
          subject, parameters.time_max);
      break;
    }
    case EndCriterion::Type::mass_below: {
      _end_criterion = std::make_unique<Criterion_mass_below<Subject>>(
          subject, parameters.end_value);
      break;
    }
    case EndCriterion::Type::mass_above: {
      _end_criterion = std::make_unique<Criterion_mass_above<Subject>>(
          subject, parameters.end_value);
      break;
    }
    case EndCriterion::Type::all_absorbed: {
      _end_criterion =
          std::make_unique<Criterion_all_absorbed<Subject>>(subject);
      break;
    }
    case EndCriterion::Type::one_absorbed: {
      _end_criterion =
          std::make_unique<Criterion_one_absorbed<Subject>>(subject);
      break;
    }
    case EndCriterion::Type::fraction_not_absorbed: {
      _end_criterion =
          std::make_unique<Criterion_fraction_not_absorbed<Subject>>(
              subject, parameters.end_value);
      break;
    }
    default:
      throw std::runtime_error{std::string{"End criterion "} +
                               parameters.end_criterion + " not supported"};
    }
  }

  template <typename SpecificMeasurement>
  SpecificMeasurement *cast(
      std::unique_ptr<OutputParameters_Cases::Measurement> const &measurement) {
    return dynamic_cast<SpecificMeasurement *>(measurement.get());
  }

  std::unique_ptr<Criterion<Subject>>
      _end_criterion; /**< To check if end criterion is met. */
  std::unique_ptr<NextMeasurementTime<std::remove_reference_t<Parameters>>>
      _next_measurement; /**< Handle next measurement time. */
  std::vector<std::unique_ptr<MeasurerTime<Subject, Geometry>>>
      _output_time; /**< Handle each output type given time. */
  std::vector<std::unique_ptr<Measurer<Subject, Geometry>>>
      _output; /**< Handle each output type given nothing. */
};
template <typename Subject, typename Geometry, typename Parameters,
          typename VelocityField, typename Mask>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Directories const &, Parameters &&, std::string const &,
             std::vector<std::reference_wrapper<const Mask>>,
             std::vector<double>)
    -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename VelocityField, typename Mask>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Directories const &, Parameters &&, std::string const &,
             std::vector<std::reference_wrapper<const Mask>>)
    -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename VelocityField>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Directories const &, Parameters &&, std::string const &)
    -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename VelocityField, typename Mask>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Directories const &, Parameters &&, std::string const &,
             std::initializer_list<std::reference_wrapper<const Mask>>,
             std::initializer_list<double>)
    -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename VelocityField, typename Mask>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Directories const &, Parameters &&, std::string const &,
             std::initializer_list<std::reference_wrapper<const Mask>>)
    -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename VelocityField, typename Mask>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Directories const &, Parameters &&, std::string const &,
             std::vector<std::reference_wrapper<const Mask>>,
             std::initializer_list<double>)
    -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename VelocityField, typename Mask>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Directories const &, Parameters &&, std::string const &,
             std::initializer_list<std::reference_wrapper<const Mask>>,
             std::vector<double>)
    -> Output_Cases<Subject, Geometry, Parameters>;
} // namespace ptof

#endif /* PTOF_OUTPUT_CASES_H */
