/**
 \file PTOF/Output_Cases.h
 \author Tomás Aquino
 \date 07/03/2022
*/

#ifndef PTOF_OUTPUT_CASES_H
#define PTOF_OUTPUT_CASES_H

#include "CTRW/Meta.h"
#include "General/Meta.h"
#include "General/Useful.h"
#include "PTOF/BoundaryConditionList.h"
#include "PTOF/Criterion.h"
#include "PTOF/Directories.h"
#include "PTOF/Field.h"
#include "PTOF/MeasurementList.h"
#include "PTOF/MeasurementSpacing.h"
#include "PTOF/Measurer.h"
#include "PTOF/MeasurerTime.h"
#include "PTOF/NextMeasurement.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <fvcGrad.H>
#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace ptof {
/** \struct OutputParameters_Cases_Cases PTOF/Output_Cases.h
 * "PTOF/Output_Cases.h" \brief Parameters for output. */
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

    template <typename OStream> void info_runtime(OStream &output) const {
      output << name;
      output << "\tPrecision: " << precision;
    }
  };

  struct Measurement_field : public Measurement {
    Measurement_field(std::string name, std::string field_name,
                      int precision = 8)
        : Measurement{name, precision}, field_name{field_name} {}

    std::string field_name;

    template <typename OStream> void info_runtime(OStream &output) const {
      Measurement::info_runtime(output);
      output << "\tField: " << field_name;
    }
  };

  struct Measurement_dim_position : public Measurement {
    Measurement_dim_position(std::string name, std::size_t dim, double position,
                             int precision = 8)
        : Measurement{name, precision}, dim{dim}, position{position} {}

    std::size_t dim{dim};
    double position{position};

    template <typename OStream> void info(OStream &output) const {
      Measurement::info_runtime(output);
      output << "\tDimension: " << dim << "\t"
             << "\tPosition: " << position << "\t";
    }
  };

  struct Measurement_order : public Measurement {
    Measurement_order(std::string name, double order, int precision = 8)
        : Measurement{name, precision}, order{order} {}

    double order;

    template <typename OStream> void info_runtime(OStream &output) const {
      Measurement::info_runtime(output);
      output << "\tOrder: " << order;
    }
  };

  struct Measurement_orders : public Measurement {
    Measurement_orders(std::string name, std::vector<double> orders,
                       int precision = 8)
        : Measurement{name, precision}, orders{orders} {}

    std::vector<double> orders;

    template <typename OStream> void info_runtime(OStream &output) const {
      Measurement::info_runtime(output);
      output << "\tOrders: ";
      useful::print(output, orders, false, " ");
    }
  };

  std::vector<std::unique_ptr<Measurement>> measurements;

  /** Constructor.
   \param directories Current case directory information.
   \param name Name of parameter set.
   \param params_transport Transport parameters.
   \param params_reaction Reaction parameters.
   \param params_solvers Solver parameters. */
  template <typename Geometry, typename TransportParameters,
            typename ReactionParameters, typename SolverParameters>
  OutputParameters_Cases(Directories const &directories,
                         std::string const &name, Geometry const &geometry,
                         TransportParameters const &params_transport,
                         ReactionParameters const &params_reaction,
                         SolverParameters const &params_solvers) {
    auto input = useful::open_read(directories.dir_parameters +
                                   "/parameters_output_" + name + ".dat");
    read_time_units(input, params_transport, params_reaction);
    read_end_criterion(input);
    read_measurement_spacing(input, params_transport, params_reaction);
    read_measurement_types<Geometry>(input);
    input.close();
  }

  /** \brief Output generic information about object. */
  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------"
           "\n"
           "Output parameters\n"
           "--------------------------------------------------------------"
           "\n"
           "- Time units to rescale measurement times:\n"
           "\tdiffusion: Rescale by diffusion time\n"
           "\tadvection: Rescale by advection time\n"
           "\treaction: Rescale by reaction time\n"
           "\tarbitrary: Do not rescale\n"
           "- End criterion to finish dynamics:\n"
           "\ttime: Specified time\n"
           "\ttime_max: Maximum output time\n"
           "\tmass_below: Total mass below value\n"
           "\tmass_above: Total mass above value\n"
           "\tall_absorbed: All particles absorbed\n"
           "\tone_absorbed: One particle absorbed\n"
           "\tfraction_not_absorbed: Fraction of particles not absorbed\n"
           "\t                       below value\n"
           "- End value, if required by end criterion\n"
           "- Measurement spacing:\n"
           "\tstep: Linear spacing, specified time step (no maximum time)\n"
           "\tlinear: Linear spacing, specified maximum time and number\n"
           "\t        of measurements\n"
           "\tlog: Log spacing, specified maximum time and number of\n"
           "\t     measurements\n"
           "- Minimum measurement time\n"
           "- Measurement time increment, if required by measure spacing\n"
           "- Maximum measurement time, if required by measure spacing\n"
           "- Number of measurements, if required by measure spacing\n"
           "- Measurement types (any number, can pass precision [8] at the\n"
           "  end of line):\n"
           "\tposition: Time, particle tags, particle positions, and\n"
           "\t          particle masses\n"
           "\tposition_in_regions: Time, particle tags, particle\n"
           "\t                     positions, and particle masses within\n"
           "\t                     regions specified by masks\n"
           "\tposition_mean: Time and mean position\n"
           "\tposition_second_moment: Time and position second moment\n"
           "\tposition_nth_moment: Time and position nth moment\n"
           "\tposition__moment: Time and position moment with orders\n"
           "\t                  specified along each dimension\n"
           "\tposition_variance: Time and position variance\n"
           "\tmass: Time and total mass\n"
           "\tmass_in_regions: Time and total mass within regions "
           "\t                 specified by masks\n"
           "\tvelocity: Time, particle tags, and local velocities\n"
           "\tvelocity: Time and mean of velocity field over particles\n"
           "\tvelocity_gradient: Time, particle tags, and local velocity\n"
           "\tvelocity_gradient_mean: Time and mean of velocity gradient\n"
           "\t                        field over particles\n"
           "\tscalar_field: Time, particle tags, and local scalar field\n"
           "\t              values (specify field name on same line;\n"
           "\t              field to be read from OF file)\n"
           "\tvector_field: Time, particle tags, and local vector field\n"
           "\t              values (specify field name on same line;\n"
           "\t              field to be read from OF file)\n"
           "\ttensor_field: Time, particle tags, and local tensor field\n"
           "\t              values (specify field name on same line; "
           "\t              field to be read from OF file)\n"
           "\tscalar_field_mean: Time and mean of scalar field over\n"
           "\t                   particles (specify field name on same line;\n"
           "\t                   field to be read from OF file)\n"
           "\tvector_field_mean: Time and mean of vector field over\n"
           "\t                   particles (specify field name on same line;\n"
           "\t                   field to be read from OF file)\n"
           "\ttensor_field_mean: Time and mean of tensor field over\n"
           "\t                   particles (specify field name on same line;\n"
           "\t                   field to be read from OF file)\n"
           "\tposition_periodic: Time, particle tags, true positions\n"
           "\t                   accounting for periodicity, and masses\n"
           "\tposition_in_regions_periodic: Time, particle tags, true\n"
           "\t                              positions accounting for\n"
           "\t                              periodicity, and masses in\n"
           "\t                              regions specified by masks\n"
           "\tposition_mean_periodic: Time and true mean position\n"
           "\t                        accounting for periodicity\n"
           "\tposition_second_moment_periodic: Time and true position\n"
           "\t                                 second moment accounting for\n"
           "\t                                 periodicity\n"
           "\tposition_nth_moment_periodic: Time and true position nth\n"
           "\t                              moment accounting for\n"
           "\t                              periodicity\n"
           "\tposition_moment: Time and true position moment with orders\n"
           "\t                 specified along each dimension accounting\n"
           "\t                 for periodicity\n"
           "\tposition_variance_periodic: Time and true position variance\n"
           "\t                            accounting for periodicity\n"
           "\tfirst_crossing_time: Crossing times, particle tags, and\n"
           "\t                     particle masses at end of dynamics\n"
           "\t                     (specify dimension and crossing\n"
           "\t                     position along dimension on same line)\n"
           "\tabsorption_time: Particle absorption times, particle tags,\n"
           "\t                 and particle masses at end of dynamics\n"
           "--------------------------------------------------------------"
           "\n";
  }

  /** \brief Read end criterion from input stream. */
  template <typename IStream> void read_end_criterion(IStream &input) {
    useful::read_first_from_line(end_criterion, input);
    switch (EndCriterion::type(end_criterion)) {
    case EndCriterion::Type::time: {
      useful::read_first_from_line(end_value, input);
      end_value *= time_unit_factor;
      break;
    }
    case EndCriterion::Type::time_max:
      break;
    case EndCriterion::Type::mass_below: {
      useful::read_first_from_line(end_value, input);
      break;
    }
    case EndCriterion::Type::mass_above: {
      useful::read_first_from_line(end_value, input);
      break;
    }
    case EndCriterion::Type::all_absorbed:
      break;
    case EndCriterion::Type::one_absorbed:
      break;
    case EndCriterion::Type::fraction_not_absorbed: {
      useful::read_first_from_line(end_value, input);
      break;
    }
    default:
      throw std::runtime_error{std::string("End criterion ") + end_criterion +
                               " not supported"};
    }
  }

  /** \brief Read measure spacing time units from input stream. */
  template <typename IStream, typename TransportParameters,
            typename ReactionParameters>
  void read_time_units(IStream &input,
                       TransportParameters const &params_transport,
                       ReactionParameters const &params_reaction) {
    useful::read_first_from_line(time_units, input);
    switch (MeasurementSpacingUnits::type(time_units)) {
    case MeasurementSpacingUnits::Type::diffusion: {
      time_unit_factor = params_transport.diffusion_time;
      break;
    }
    case MeasurementSpacingUnits::Type::advection: {
      time_unit_factor = params_transport.advection_time;
      break;
    }
    case MeasurementSpacingUnits::Type::reaction: {
      time_unit_factor = params_reaction.reaction_time;
      break;
    }
    case MeasurementSpacingUnits::Type::arbitrary: {
      time_unit_factor = 1.;
      break;
    }
    default:
      throw std::runtime_error{std::string("Measurement spacing units ") +
                               time_units + " not supported"};
    }
  }

  /** \brief Read measure spacing type from input stream. */
  template <typename IStream, typename TransportParameters,
            typename ReactionParameters>
  void read_measurement_spacing(IStream &input,
                                TransportParameters const &params_transport,
                                ReactionParameters const &params_reaction) {
    useful::read_first_from_line(measurement_spacing, input);
    useful::read_first_from_line(time_min, input);
    time_min *= time_unit_factor;
    switch (MeasurementSpacing::type(measurement_spacing)) {
    case MeasurementSpacing::Type::step: {
      useful::read_first_from_line(time_increment, input);
      time_increment *= time_unit_factor;
      time_max = std::numeric_limits<double>::infinity();
      if (!(time_increment > 0.))
        throw std::runtime_error{
            std::string("Non-positive time increment for ") +
            measurement_spacing + " measurement spacing"};
      break;
    }
    case MeasurementSpacing::Type::linear: {
      useful::read_first_from_line(time_max, input);
      time_max *= time_unit_factor;
      std::size_t nr_measures;
      useful::read_first_from_line(nr_measures, input);
      if (nr_measures < 2)
        throw std::runtime_error{
            std::string("Minimum of two measurements required for ") +
            measurement_spacing + " measurement spacing"};
      time_increment = (time_max - time_min) / (nr_measures - 1);
      if (!(time_increment > 0.))
        throw std::runtime_error{
            std::string("Non-positive time increment for ") +
            measurement_spacing + " measurement spacing"};
      break;
    }
    case MeasurementSpacing::Type::log: {
      if (!(time_min > 0.))
        throw std::runtime_error{
            std::string("Non-positive minimum measurement time for ") +
            measurement_spacing + " measurement spacing"};
      useful::read_first_from_line(time_max, input);
      time_max *= time_unit_factor;
      std::size_t nr_measures;
      useful::read_first_from_line(nr_measures, input);
      if (nr_measures < 2)
        throw std::runtime_error{
            std::string("Minimum of two measurements required for ") +
            measurement_spacing + " measurement spacing"};
      time_increment = std::pow(time_max / time_min, 1. / (nr_measures - 1));
      if (!(time_increment > 1.))
        throw std::runtime_error{
            std::string("Non-positive time increment for ") +
            measurement_spacing + " measurement spacing"};
      break;
    }
    default:
      throw std::runtime_error{std::string("Measurement spacing ") +
                               measurement_spacing + " not supported"};
    }
  }

  /** \brief Read measure spacing type from input stream. */
  template <typename Geometry, typename IStream>
  void read_measurement_types(IStream &input) {
    std::vector<std::string> split_line;
    while (useful::split_line(input, split_line)) {
      if (split_line.empty())
        continue;
      std::size_t required_size = 1;
      if (split_line.size() < required_size)
        throw std::runtime_error{"Could not parse measurement type"};
      std::string const &name = split_line[0];
      if (!MeasurementList{}.contains(name))
        throw std::runtime_error{"Output type " + name + " not supported"};
      if (name == "scalar_field" || name == "vector_field" ||
          name == "tensor_field") {
        required_size = 2;
        if (split_line.size() < required_size)
          throw std::runtime_error{std::string("Measurement type ") + name +
                                   ": "
                                   "Could not parse field name"};
        measurements.emplace_back(
            std::make_unique<Measurement_field>(name, split_line[1]));
      } else if (name == "first_crossing_time") {
        required_size = 3;
        if (split_line.size() < required_size)
          throw std::runtime_error{
              std::string("Measurement type ") + name +
              ": "
              "Could not parse crossing dimension and position"};
        std::size_t dim = std::stoul(split_line[1]);
        if (dim >= Geometry::dim)
          throw std::runtime_error{
              std::string("Measurement type ") + name +
              ": "
              "Crossing dimension higher than simulation dimension"};
        measurements.emplace_back(std::make_unique<Measurement_dim_position>(
            name, dim, std::stod(split_line[2])));
      } else if (name == "position_nth_moment" ||
                 name == "position_nth_moment_periodic") {
        required_size = 2;
        if (split_line.size() < required_size)
          throw std::runtime_error{std::string("Measurement type ") + name +
                                   ": "
                                   "Could not parse moment order"};
        measurements.emplace_back(std::make_unique<Measurement_order>(
            name, std::stod(split_line[1])));
      } else if (name == "position_moment" ||
                 name == "position_moment_periodic") {
        required_size = Geometry::dim + 1;
        if (split_line.size() < required_size)
          throw std::runtime_error{
              std::string("Measurement type ") + name +
              ": "
              "Could not parse coordinate moment orders along each dimension"};
        std::vector<double> orders;
        orders.reserve(Geometry::dim);
        for (std::size_t ii = 1; ii < Geometry::dim + 1; ++ii)
          orders.push_back(std::stod(split_line[ii]));
        measurements.emplace_back(
            std::make_unique<Measurement_orders>(name, orders));
      } else {
        measurements.emplace_back(std::make_unique<Measurement>(name));
      }
      if (split_line.size() > required_size)
        measurements.back()->precision = std::stod(split_line.back());
    }
  }
};

/** \class Output_Cases PTOF/Output_Cases.h "PTOF/Output_Cases.h"
 * \brief Output object to handle implemented output options. */
template <typename Subject, typename Geometry> class Output_Cases {
public:
  using Parameters = OutputParameters_Cases; /**< Output parameters. */
  Parameters const &parameters;              /**< Output parameters .*/

  /** Constructor.
   \param subject CTRW object to measure.
   \param velocity_field Velocity field as a function of state.
   \param geometry Domain geometry info and utilities.
   \param directories Current case directory information.
   \param parameters Output parameters.
   \param identifier String to include in names of output files.
   \param masks Container of mask reference wrappers. Masks are scalar
   fields assigned values to mesh cells through operator[]. \param
   thresholds Vector of thresholds for each mask, such that cells where a
   mask is above the threshold are considered.
  */
  template <typename VelocityField, typename Mask = useful::Empty>
  Output_Cases(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters const &parameters, std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks = {},
               std::vector<double> thresholds = {})
      : parameters{parameters} {
    set_measurement_types(subject, geometry, directories, identifier,
                          velocity_field, masks, thresholds);
    set_end_criterion(subject);
    set_next_measurement_time();
  }

  template <typename VelocityField, typename Mask = useful::Empty>
  Output_Cases(Subject const &, VelocityField const &, Geometry const &,
               Directories const &, Parameters &&, std::string const &,
               std::vector<std::reference_wrapper<const Mask>> = {},
               std::vector<double> = {}) = delete;

  /** Constructor.
   \brief Overload for initializer list arguments.
  */
  template <typename VelocityField, typename Mask>
  Output_Cases(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters const &parameters, std::string const &identifier,
               std::initializer_list<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {})
      : Output_Cases(subject, velocity_field, geometry, directories, parameters,
                     identifier,
                     std::vector<std::reference_wrapper<const Mask>>{masks},
                     std::vector<double>{thresholds}) {}

  /** Constructor.
   \brief Overload for initializer list arguments.
  */
  template <typename VelocityField, typename Mask>
  Output_Cases(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters const &parameters, std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {})
      : Output_Cases(subject, velocity_field, geometry, directories, parameters,
                     identifier, masks, std::vector<double>{thresholds}) {}

  /** Constructor.
   \brief Overload without velocity field
  */
  template <typename Mask = useful::Empty>
  Output_Cases(Subject const &subject, Geometry const &geometry,
               Directories const &directories, Parameters const &parameters,
               std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks = {},
               std::vector<double> thresholds = {})
      : Output_Cases(subject, useful::Empty{}, geometry, directories,
                     parameters, identifier, masks, thresholds) {
    set_measurement_types(directories, identifier, useful::Empty{}, masks,
                          thresholds);
    set_end_criterion();
    set_next_measurement_time();
  }

  /** Constructor.
   \brief Overload for initializer list arguments.
  */
  template <typename VelocityField, typename Mask>
  Output_Cases(Subject const &subject, Geometry const &geometry,
               Directories const &directories, Parameters const &parameters,
               std::string const &identifier,
               std::initializer_list<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {})
      : Output_Cases(subject, geometry, directories, parameters, identifier,
                     std::vector<std::reference_wrapper<const Mask>>{masks},
                     std::vector<double>{thresholds}) {}

  /** Constructor.
   \brief Overload for initializer list arguments.
  */
  template <typename VelocityField, typename Mask>
  Output_Cases(Subject const &subject, Geometry const &geometry,
               Directories const &directories, Parameters const &parameters,
               std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {})
      : Output_Cases(subject, geometry, directories, parameters, identifier,
                     masks, std::vector<double>{thresholds}) {}

  /** Constructor.
   \brief Overload for initializer list arguments.
  */
  template <typename VelocityField, typename Mask>
  Output_Cases(Subject const &subject, Geometry const &geometry,
               Directories const &directories, Parameters const &parameters,
               std::string const &identifier,
               std::initializer_list<std::reference_wrapper<const Mask>> masks,
               std::vector<double> thresholds = {})
      : Output_Cases(subject, geometry, directories, parameters, identifier,
                     std::vector<std::reference_wrapper<const Mask>>{masks},
                     thresholds) {}

  /** \return \c true  if end simulation criterion is satisfied, \c false
   * otherwise. */
  bool done(double time) const { return _end_criterion->operator()(time); }

  /** \brief Set time of next measurement. */
  void set_next_measurement_time() {
    switch (MeasurementSpacing::type(parameters.measurement_spacing)) {
    case MeasurementSpacing::Type::step: {
      _next_measurement =
          std::make_unique<NextMeasurementTime_step<Parameters>>(parameters);
      break;
    }
    case MeasurementSpacing::Type::linear: {
      _next_measurement =
          std::make_unique<NextMeasurementTime_linear<Parameters>>(parameters);
      break;
    }
    case MeasurementSpacing::Type::log: {
      _next_measurement =
          std::make_unique<NextMeasurementTime_log<Parameters>>(parameters);
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

  /** \brief Output requested measurements at given time and advance to next
   * measurement. */
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
  template <typename OStream> void info_runtime(OStream &output) const {
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
        MeasurementSpacing::type(parameters.measurement_spacing) ==
            MeasurementSpacing::Type::step)
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
  template <typename OStream>
  void info_runtime_measurements(OStream &output) const {
    if (parameters.measurements.empty()) {
      output << "\tNone\n";
      return;
    }
    for (auto const &measurement : parameters.measurements) {
      output << "\t";
      measurement->info_runtime(output);
      output << "\n";
    }
  }

  /** \brief Set up output streams for requested output types.
   \param directories Current case directory information.
   \param identifier String to include in names of output files.
   \param masks Container of mask reference wrappers. Masks are scalar
   fields assigned values to mesh cells through operator[]. \param
   thresholds Vector of thresholds for each mask, such that cells where a
   mask is above the threshold are considered. */
  template <typename VelocityField = useful::Empty,
            typename Mask = useful::Empty>
  void set_measurement_types(
      Subject const &subject, Geometry const &geometry,
      Directories const &directories, std::string const &identifier,
      VelocityField const &velocity_field = {},
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
        if constexpr (!std::is_same_v<Mask, useful::Empty>)
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_in_regions<Subject, Geometry, Mask>>(
                  subject, geometry, directories, identifier, masks, thresholds,
                  measurement->precision));
        else
          throw std::runtime_error{std::string("Measurement type ") +
                                   measurement->name +
                                   ": "
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
        auto measurement_derived =
            dynamic_cast<OutputParameters_Cases::Measurement_order *>(
                measurement.get());
        if (measurement_derived->order == 0)
          throw std::runtime_error{std::string("Measurement type ") +
                                   measurement->name +
                                   ": "
                                   "Moment order not provided"};
        _output_time.emplace_back(
            std::make_unique<
                MeasurerTime_position_nth_moment<Subject, Geometry>>(
                subject, measurement_derived->order, geometry, directories,
                identifier, measurement->precision));
        break;
      }
      case MeasurementList::Type::position_moment: {
        auto measurement_derived =
            dynamic_cast<OutputParameters_Cases::Measurement_orders *>(
                measurement.get());
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
        if constexpr (!std::is_same_v<Mask, useful::Empty>)
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_mass_in_regions<Subject, Geometry, Mask>>(
                  subject, geometry, directories, identifier, masks, thresholds,
                  measurement->precision));
        else
          throw std::runtime_error{std::string("Measurement type ") +
                                   measurement->name +
                                   ": "
                                   "Region masks not provided"};
        break;
      }
      case MeasurementList::Type::velocity: {
        if constexpr (!std::is_same_v<VelocityField, useful::Empty>)
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_vector_field<
                  Subject, Geometry, VelocityField const &>>(
                  subject, velocity_field, geometry, directories, identifier,
                  "U", measurement->precision));
        else
          throw std::runtime_error{std::string("Measurement type ") +
                                   measurement->name + ": " +
                                   "Velocity field not provided"};
        break;
      }
      case MeasurementList::Type::velocity_mean: {
        if constexpr (!std::is_same_v<VelocityField, useful::Empty>)
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_vector_field_mean<
                  Subject, Geometry, VelocityField const &>>(
                  subject, velocity_field, geometry, directories, identifier,
                  "U", measurement->precision));
        else
          throw std::runtime_error{std::string("Measurement type ") +
                                   measurement->name + ": " +
                                   "Velocity field not provided"};
        break;
      }
      case MeasurementList::Type::velocity_gradient: {
        if constexpr (!std::is_same_v<VelocityField, useful::Empty>)
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_tensor_field<Subject, Geometry>>(
                  subject, Foam::fvc::grad(velocity_field.field()), geometry,
                  directories, identifier, "gradU", measurement->precision));
        else
          throw std::runtime_error{std::string("Measurement type ") +
                                   measurement->name + ": " +
                                   "Velocity field not provided"};
        break;
      }
      case MeasurementList::Type::velocity_gradient_mean: {
        if constexpr (!std::is_same_v<VelocityField, useful::Empty>)
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_tensor_field_mean<Subject, Geometry>>(
                  subject, Foam::fvc::grad(velocity_field.field()), geometry,
                  directories, identifier, "gradU", measurement->precision));
        else
          throw std::runtime_error{std::string("Measurement type ") +
                                   measurement->name + ": " +
                                   "Velocity field not provided"};
        break;
      }
      case MeasurementList::Type::scalar_field: {
        auto measurement_derived =
            dynamic_cast<OutputParameters_Cases::Measurement_field *>(
                measurement.get());
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_scalar_field<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::vector_field: {
        auto measurement_derived =
            dynamic_cast<OutputParameters_Cases::Measurement_field *>(
                measurement.get());
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_vector_field<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::tensor_field: {
        auto measurement_derived =
            dynamic_cast<OutputParameters_Cases::Measurement_field *>(
                measurement.get());
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_tensor_field<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::scalar_field_mean: {
        auto measurement_derived =
            dynamic_cast<OutputParameters_Cases::Measurement_field *>(
                measurement.get());
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_scalar_field_mean<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::vector_field_mean: {
        auto measurement_derived =
            dynamic_cast<OutputParameters_Cases::Measurement_field *>(
                measurement.get());
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_vector_field_mean<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived->field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::tensor_field_mean: {
        auto measurement_derived =
            dynamic_cast<OutputParameters_Cases::Measurement_field *>(
                measurement.get());
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
              std::string("Measurement type ") + measurement->name +
              ": "
              "Particle state does not define periodicity"};
        break;
      }
      case MeasurementList::Type::position_in_regions_periodic: {
        if constexpr (meta::has_periodicity_v<
                          typename Subject::Particle::State> &&
                      !std::is_same_v<Mask, useful::Empty>)
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_position_in_regions_periodic<
                  Subject, Geometry, Mask>>(subject, geometry, directories,
                                            identifier, masks, thresholds,
                                            measurement->precision));
        else {
          if constexpr (!meta::has_periodicity_v<
                            typename Subject::Particle::State>)
            throw std::runtime_error{
                std::string("Measurement type ") + measurement->name +
                ": "
                "Particle state does not define periodicity"};
          throw std::runtime_error{std::string("Measurement type ") +
                                   measurement->name +
                                   ": "
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
              std::string("Measurement type ") + measurement->name +
              ": "
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
              std::string("Measurement type ") + measurement->name +
              ": "
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
              std::string("Measurement type ") + measurement->name +
              ": "
              "Particle state does not define periodicity"};
        break;
      }
      case MeasurementList::Type::position_nth_moment_periodic: {
        if constexpr (meta::has_periodicity_v<
                          typename Subject::Particle::State>) {
          auto measurement_derived =
              dynamic_cast<OutputParameters_Cases::Measurement_order *>(
                  measurement.get());
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_nth_moment_periodic<Subject, Geometry>>(
                  subject, measurement_derived->order, geometry, directories,
                  identifier, measurement->precision));
        } else
          throw std::runtime_error{
              std::string("Measurement type ") + measurement->name +
              ": "
              "Particle state does not define periodicity"};
        break;
      }
      case MeasurementList::Type::position_moment_periodic: {
        if constexpr (meta::has_periodicity_v<
                          typename Subject::Particle::State>) {
          auto measurement_derived =
              dynamic_cast<OutputParameters_Cases::Measurement_orders *>(
                  measurement.get());
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_moment_periodic<Subject, Geometry>>(
                  subject, measurement_derived->orders, geometry, directories,
                  identifier, measurement->precision));
        } else
          throw std::runtime_error{
              std::string("Measurement type ") + measurement->name +
              ": "
              "Particle state does not define periodicity"};
        break;
      }
      case MeasurementList::Type::first_crossing_time: {
        auto measurement_derived =
            dynamic_cast<OutputParameters_Cases::Measurement_dim_position *>(
                measurement.get());
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
        throw std::runtime_error{std::string("Measurement type ") +
                                 measurement->name +
                                 " "
                                 "not supported"};
      }
    }
  }

  /** Set end criterion. */
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
            std::string("Infinite maximum time with end criterion ") +
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
      throw std::runtime_error{std::string("End criterion ") +
                               parameters.end_criterion + " not supported"};
    }
  }

  std::unique_ptr<Criterion<Subject>>
      _end_criterion; /**< To check if end criterion is met. */
  std::unique_ptr<NextMeasurementTime<Parameters>>
      _next_measurement; /**< Handle next measure time. */
  std::vector<std::unique_ptr<MeasurerTime<Subject, Geometry>>>
      _output_time; /**< Handle each output type given time. */
  std::vector<std::unique_ptr<Measurer<Subject, Geometry>>>
      _output; /**< Handle each output type given nothing. */
};
} // namespace ptof

#endif /* PTOF_OUTPUT_CASES_H */
