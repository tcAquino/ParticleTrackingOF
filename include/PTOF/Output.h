/**
 \file PTOF/Output.h
 \author Tomás Aquino
 \date 07/03/2022
*/

#ifndef PTOF_OUTPUT_H
#define PTOF_OUTPUT_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
#include <fvcGrad.H>
#include "CTRW/Meta.h"
#include "General/Useful.h"
#include "PTOF/Criteria.h"
#include "PTOF/Directories.h"
#include "PTOF/Field.h"
#include "PTOF/Measure.h"
#include "PTOF/Measurer.h"
#include "PTOF/MeasurerTime.h"

namespace ptof
{
  /** \struct OutputParameters PTOF/Output.h "PTOF/Output.h"
   *  \brief Parameters for output. */
  struct OutputParameters
  {
    std::string time_units;
    double time_unit_factor;
    std::string end_criterion;
    double end_value;
    std::string measure_spacing;
    double time_min;
    double time_max;
    double time_increment;
    struct Measurement
    {
      Measurement(std::string name, std::string field_name = {}, int precision = 8)
      : name{ name }
      , field_name{ field_name }
      , precision{ precision }
      {}
      
      Measurement(std::string name, int precision)
      : name{ name }
      , field_name{}
      , precision{ precision }
      {}
      
      std::string name;
      std::string field_name;
      int precision;
    };
    std::vector<Measurement> measurements;
    
    /** Constructor.
     \param directories Current case directory information.
     \param name Name of parameter set.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters. */
    template
    <typename Geometry,
    typename TransportParameters,
    typename ReactionParameters,
    typename SolverParameters>
    OutputParameters
    (Directories const& directories,
     std::string const& name,
     Geometry const& geometry,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers)
    {
      auto input = useful::open_read(directories.dir_parameters
                                     + "/parameters_output_"
                                     + name + ".dat");
      read_time_units(input, params_transport, params_reaction);
      read_end_criterion(input);
      read_measure_spacing(input,
                           params_transport, params_reaction);
      input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
      read_measurement_types(input);
      input.close();
    }
    
    /** \brief Output generic information about object. */
    template <typename OStream>
    static void info(OStream& output)
    {
      output <<
        "--------------------------------------------------\n"
        "Output parameters\n"
        "--------------------------------------------------\n"
        "- Time units to rescale measurement times:\n"
        "\tdiffusion: Rescale by diffusion time\n"
        "\tadvection: Rescale by reaction time\n"
        "\treaction: Rescale by reaction time\n"
        "\tarbitrary: Do not rescale\n"
        "- End criterion to finish dynamics:\n"
        "\ttime: Specified time\n"
        "\ttime_max: Maximum output time\n"
        "\tmass_below: Total mass below value\n"
        "\tmass_above: Total mass above value\n"
        "\tall_absorbed: All particles absorbed\n"
        "\tone_absorbed: One particle absorbed\n"
        "\tfraction_not_absorbed: Fraction of particles not absorbed below value\n"
        "- End value, if required by end criterion\n"
        "- Measure spacing:\n"
        "\tstep: Linear spacing, specified time step (no maximum time)\n"
        "\tlinear: Linear spacing, specified maximum time and number of measurements\n"
        "\tlog: Log spacing, specified maximum time and number of measurements\n"
        "- Minimum measurement time\n"
        "- Measurement time increment, if required by measure spacing\n"
        "- Maximum measurement time, if required by measure spacing\n"
        "- Number of measurements, if required by measure spacing\n"
        "- Measurement types (any number, can pass precision [8] at the end of line):\n"
        "\tposition: Time, particle tags, particle positions, and particle masses\n"
        "\tposition: Time, particle tags, particle positions, and particle masses within regions specified by masks\n"
        "\tposition_mean: Time and mean position\n"
        "\tposition_second_moment: Time and position second moment\n"
        "\tposition_variance: Time and position variance\n"
        "\tmass: Time and total mass\n"
        "\tmass_in_regions: Time and total mass within regions specified by masks\n"
        "\tvelocity: Time, particle tags, and local velocities\n"
        "\tvelocity_gradient: Time, particle tags, and local velocity gradients\n"
        "\tscalar_field: Time, particle tags, and local scalar field values (Specify field name on same line, to be read from OF file)\n"
        "\tvector_field: Time, particle tags, and local vector field values (Specify field name on same line, to be read from OF file)\n"
        "\ttensor_field: Time, particle tags, and local tensor field values (Specify field name on same line, to be read from OF file)\n"
        "\tposition_periodic: Time, particle tags, true positions accounting for periodicity, and masses\n"
        "\tposition_in_regions_periodic: Time, particle tags, true positions accounting for periodicity, and masses in regions speciefied by masks\n"
        "\tposition_mean_periodic: Time and true mean position accounting for periodicity\n"
        "\tposition_second_moment_periodic: Time and true position second moment accounting for periodicity\n"
        "\tposition_variance_periodic: Time and true position variance accounting for periodicity\n"
        "\tabsorption_time: Particle absorption times, particle tags, and particle masses at end of dynamics\n"
        "--------------------------------------------------\n";
    }
    
  private:
    /** \brief Read end criterion from input stream. */
    template <typename IStream>
    void read_end_criterion(IStream& input)
    {
      useful::read(input, end_criterion);
      switch (EndCriterion::type(end_criterion))
      {
        case EndCriterion::Type::time:
        {
          useful::read(input, end_value);
          end_value *= time_unit_factor;
          break;
        }
        case EndCriterion::Type::time_max:
          break;
        case EndCriterion::Type::mass_below:
        {
          useful::read(input, end_value);
          break;
        }
        case EndCriterion::Type::mass_above:
        {
          useful::read(input, end_value);
          break;
        }
        case EndCriterion::Type::all_absorbed:
          break;
        case EndCriterion::Type::one_absorbed:
          break;
        case EndCriterion::Type::fraction_not_absorbed:
        {
          useful::read(input, end_value);
          break;
        }
        default:
          throw std::runtime_error{
            std::string("End criterion ")
            + end_criterion
            + " not supported" };
      }
    }
    
    /** \brief Read measure spacing time units from input stream. */
    template <typename IStream,
    typename TransportParameters,
    typename ReactionParameters>
    void read_time_units
    (IStream& input,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction)
    {
      useful::read(input, time_units);
      switch (MeasureSpacingUnits::type(time_units))
      {
        case MeasureSpacingUnits::Type::diffusion:
        {
          time_unit_factor = params_transport.diffusion_time;
          break;
        }
        case MeasureSpacingUnits::Type::advection:
        {
          time_unit_factor = params_transport.advection_time;
          break;
        }
        case MeasureSpacingUnits::Type::reaction:
        {
          time_unit_factor = params_reaction.reaction_time;
          break;
        }
        case MeasureSpacingUnits::Type::arbitrary:
        {
          time_unit_factor = 1.;
          break;
        }
        default:
          throw std::runtime_error{
            std::string("Measure spacing units ")
            + time_units
            + " not supported" };
      }
    }
    
    /** \brief Read measure spacing type from input stream. */
    template <typename IStream,
    typename TransportParameters,
    typename ReactionParameters>
    void read_measure_spacing
    (IStream& input,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction)
    {
      useful::read(input, measure_spacing);
      useful::read(input, time_min);
      time_min *= time_unit_factor;
      switch (MeasureSpacing::type(measure_spacing))
      {
        case MeasureSpacing::Type::step:
        {
          useful::read(input, time_increment);
          time_increment *= time_unit_factor;
          time_max = std::numeric_limits<double>::infinity();
          break;
        }
        case MeasureSpacing::Type::linear:
        {
          useful::read(input, time_max);
          time_max *= time_unit_factor;
          std::size_t nr_measures;
          useful::read(input, nr_measures);
          time_increment = (time_max - time_min)/(nr_measures - 1);
          break;
        }
        case MeasureSpacing::Type::log:
        {
          useful::read(input, time_max);
          time_max *= time_unit_factor;
          std::size_t nr_measures;
          useful::read(input, nr_measures);
          time_increment = std::pow(time_max/time_min, 1./(nr_measures - 1));
          break;
        }
        default:
          throw std::runtime_error{
            std::string("Measure spacing ")
            + measure_spacing
            + " not supported" };
      }
    }
    
    /** \brief Read measure spacing type from input stream. */
    template <typename IStream>
    void read_measurement_types
    (IStream& input)
    {
      std::string line;
      while (std::getline(input, line))
      {
        std::cout << line << std::endl;
        auto split_line = useful::split(useful::remove_carriage_return(line));
        std::size_t required_size = 1;
        if (split_line.size() < required_size)
          throw std::runtime_error{
            "Could not parse measurement types" };
        std::string const& name = split_line[0];
        measurements.emplace_back(name);
        if (name == "scalar_field" ||
            name == "vector_field" ||
            name == "tensor_field")
        {
          required_size = 2;
          if (split_line.size() < required_size)
            throw std::runtime_error{
              std::string("Measurement type ") + name + ": "
                "Field name not provided" };
          measurements.back().field_name = split_line[1];
        }
        if (split_line.size() > required_size)
          measurements.back().precision = std::stod(split_line.back());
      }
    }
  };
  
  /** \class Output_Cases PTOF/Output.h "PTOF/Output.h"
   * \brief Output object to handle implemented output options. */
  template
  <typename Subject, typename Geometry>
  class Output_Cases
  {
  public:
    using Parameters = OutputParameters;    /**< Output parameters. */
    
    /** Constructor.
     \param subject CTRW object to measure.
     \param velocity_field Velocity field as a function of state.
     \param geometry Domain geometry info and utilities.
     \param directories Current case directory information.
     \param parameters Output parameters.
     \param identifier String to include in names of output files.
     \param masks Container of mask reference wrappers. Masks are scalar fields assigned values to mesh cells through operator[].
     \param thresholds Vector of thresholds for each mask, such that cells where a mask is above the threshold are considered.
    */
    template <typename VelocityField, typename Mask = useful::Empty>
    Output_Cases
    (Subject const& subject,
     VelocityField const& velocity_field,
     Geometry const& geometry,
     Directories const& directories,
     Parameters parameters,
     std::string const& identifier,
     std::vector<std::reference_wrapper<const Mask>> masks = {},
     std::vector<double> thresholds = {})
    : parameters{ parameters }
    {
      set_measure_types(subject, geometry,
                        directories, identifier,
                        velocity_field,
                        masks, thresholds);
      set_end_criterion(subject);
      set_next_measure_time();
    }
    
    /** Constructor.
     \brief Overload for initializer list arguments.
    */
    template <typename VelocityField, typename Mask>
    Output_Cases
    (Subject const& subject,
     VelocityField const& velocity_field,
     Geometry const& geometry,
     Directories const& directories,
     Parameters parameters,
     std::string const& identifier,
     std::initializer_list<std::reference_wrapper<const Mask>> masks,
     std::initializer_list<double> thresholds = {})
    : Output_Cases(subject, velocity_field, geometry, directories,
                   parameters, identifier,
                   std::vector<std::reference_wrapper<const Mask>>{ masks },
                   std::vector<double>{ thresholds })
    {}
    
    /** Constructor.
     \brief Overload for initializer list arguments.
    */
    template <typename VelocityField, typename Mask>
    Output_Cases
    (Subject const& subject,
     VelocityField const& velocity_field,
     Geometry const& geometry,
     Directories const& directories,
     Parameters parameters,
     std::string const& identifier,
     std::vector<std::reference_wrapper<const Mask>> masks,
     std::initializer_list<double> thresholds = {})
    : Output_Cases(subject, velocity_field, geometry, directories,
                   parameters, identifier,
                   masks, std::vector<double>{ thresholds })
    {}
    
    /** Constructor.
     \brief Overload without velocity field
    */
    template <typename Mask = useful::Empty>
    Output_Cases
    (Subject const& subject,
     Geometry const& geometry,
     Directories const& directories,
     Parameters parameters,
     std::string const& identifier,
     std::vector<std::reference_wrapper<const Mask>> masks = {},
     std::vector<double> thresholds = {})
    : Output_Cases(subject, useful::Empty{}, geometry, directories,
      parameters, identifier,
      masks, thresholds)
    {
      set_measure_types(directories, identifier,
                        useful::Empty{},
                        masks, thresholds);
      set_end_criterion();
      set_next_measure_time();
    }
    
    /** Constructor.
     \brief Overload for initializer list arguments.
    */
    template <typename VelocityField, typename Mask>
    Output_Cases
    (Subject const& subject,
     Geometry const& geometry,
     Directories const& directories,
     Parameters parameters,
     std::string const& identifier,
     std::initializer_list<std::reference_wrapper<const Mask>> masks,
     std::initializer_list<double> thresholds = {})
    : Output_Cases(subject, geometry, directories,
                   parameters, identifier,
                   std::vector<std::reference_wrapper<const Mask>>{ masks },
                   std::vector<double>{ thresholds })
    {}
    
    /** Constructor.
     \brief Overload for initializer list arguments.
    */
    template <typename VelocityField, typename Mask>
    Output_Cases
    (Subject const& subject,
     Geometry const& geometry,
     Directories const& directories,
     Parameters parameters,
     std::string const& identifier,
     std::vector<std::reference_wrapper<const Mask>> masks,
     std::initializer_list<double> thresholds = {})
    : Output_Cases(subject, geometry, directories,
                   parameters, identifier,
                   masks,
                   std::vector<double>{ thresholds })
    {}
    
    /** Constructor.
     \brief Overload for initializer list arguments.
    */
    template <typename VelocityField, typename Mask>
    Output_Cases
    (Subject const& subject,
     Geometry const& geometry,
     Directories const& directories,
     Parameters parameters,
     std::string const& identifier,
     std::initializer_list<std::reference_wrapper<const Mask>> masks,
     std::vector<double> thresholds = {})
    : Output_Cases(subject, geometry, directories,
                   parameters, identifier,
                   std::vector<std::reference_wrapper<const Mask>>{ masks },
                   thresholds)
    {}
    
    /** \return \c true  if end simulation criterion is satisfied, \c false otherwise. */
    bool done(double time) const
    {
      return _end_criterion->operator()(time);
    }

    /** \brief Set time of next measurement. */
    void set_next_measure_time()
    {
      switch (MeasureSpacing::type(parameters.measure_spacing))
      {
        case MeasureSpacing::Type::step:
        {
          _next_measure = std::make_unique<
            NextMeasureTime_step>(parameters);
          break;
        }
        case MeasureSpacing::Type::linear:
        {
          _next_measure = std::make_unique<
            NextMeasureTime_linear>(parameters);
          break;
        }
        case MeasureSpacing::Type::log:
        {
          _next_measure = std::make_unique<
            NextMeasureTime_log>(parameters);
          break;
        }
        default:
          throw std::runtime_error{
            "Measure spacing "
            + parameters.measure_spacing
            + " not supported" };
      }
    }
    
    /** \return Time of next measurement. */
    double next_measure_time() const
    { return _next_measure->time(); }
        
    /** \brief Output requested measurements at given time and advance to next measurement. */
    void operator()(double time)
    {
      for (auto const& output : _output_time)
        output->operator()(time);
      _next_measure->advance();
    }
    
    /** \brief Output current information. */
    void operator()()
    {
      for (auto const& output : _output)
        output->operator()();
    }
    
    /** \brief Update internal state. */
    void update(double time, double time_of_change)
    {
      for (auto const& output : _output_time)
        output->update(time, time_of_change);
    }
    
    const Parameters parameters;  /**< Output parameters .*/
    
    /** \brief Output information about current object. */
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Output\n"
        "--------------------------------------------------\n"
        "- End criterion: " << parameters.end_criterion << "\n";
      if (EndCriterion::type(parameters.end_criterion) == EndCriterion::Type::time
          || EndCriterion::type(parameters.end_criterion) == EndCriterion::Type::mass_below
          || EndCriterion::type(parameters.end_criterion) == EndCriterion::Type::mass_above
          || EndCriterion::type(parameters.end_criterion) == EndCriterion::Type::fraction_not_absorbed)
        output << "- End value: " << parameters.end_value << "\n";
      output <<
        "- Time units: " << parameters.time_units << "\n"
        "- Measure spacing: " << parameters.measure_spacing << "\n"
        "- Minimum measurement time: " << parameters.time_min << "\n"
        "- Maximum measurement time: ";
      if (EndCriterion::type(parameters.end_criterion) == EndCriterion::Type::time
          && MeasureSpacing::type(parameters.measure_spacing) == MeasureSpacing::Type::step)
        output << parameters.end_value << "\n";
      else
        output << parameters.time_max << "\n";
      output << "- Measurement types:\n";
      info_runtime_measurements(output, "\t");
      output <<
        "--------------------------------------------------\n";
    }
        
  private:
    /** \brief Output information about measurement types. */
    template <typename OStream>
    void info_runtime_measurements(OStream& output, std::string const& prefix = "") const
    {
      if (parameters.measurements.empty())
      {
        output << "None\n";
        return;
      }
      int width_name = int(std::max_element(parameters.measurements.begin(),
                                            parameters.measurements.end(),
                                            [](auto const& aa, auto const& bb)
                                            { return aa.name.length() < bb.name.length(); })
                           ->name.length()) + 1;
      int width_field = int(std::max_element(parameters.measurements.begin(),
                                             parameters.measurements.end(),
                                             [](auto const& aa, auto const& bb)
                                             { return aa.field_name.length() < bb.field_name.length(); })
                            ->field_name.length()) + 1;
      for (auto const& measurement : parameters.measurements)
      {
        output << prefix
               << std::left << std::setw(width_name) << measurement.name;
        if (!measurement.field_name.empty())
          output << std::left << std::setw(width_field)
                 << measurement.field_name;
        output << "\n";
      }
    }
    
    /** \brief Set up output streams for requested output types.
     \param directories Current case directory information.
     \param identifier String to include in names of output files.
     \param masks Container of mask reference wrappers. Masks are scalar fields assigned values to mesh cells through operator[].
     \param thresholds Vector of thresholds for each mask, such that cells where a mask is above the threshold are considered. */
    template <typename VelocityField = useful::Empty, typename Mask = useful::Empty>
    void set_measure_types
    (Subject const& subject,
     Geometry const& geometry,
     Directories const& directories,
     std::string const& identifier,
     VelocityField const& velocity_field = {},
     std::vector<std::reference_wrapper<const Mask>> masks = {},
     std::vector<double> thresholds = {})
    {
      for (auto const& measurement : parameters.measurements)
      {
        switch (Measure::type(measurement.name))
        {
          case Measure::Type::position:
          {
            _output_time.emplace_back
            (std::make_unique<MeasurerTime_position<Subject, Geometry>>
             (subject, geometry,
              directories, identifier,
              measurement.precision));
            break;
          }
          case Measure::Type::position_in_regions:
          {
            if constexpr (!std::is_same_v<Mask, useful::Empty>)
              _output_time.emplace_back
              (std::make_unique<MeasurerTime_position_in_regions<Subject, Geometry, Mask>>
               (subject, geometry,
                directories, identifier,
                masks, thresholds,
                measurement.precision));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + measurement.name + ": "
                "Region masks not specified" };
            break;
          }
          case Measure::Type::position_mean:
          {
            _output_time.emplace_back
            (std::make_unique<MeasurerTime_position_mean<Subject, Geometry>>
             (subject, geometry,
              directories, identifier,
              measurement.precision));
            break;
          }
          case Measure::Type::position_second_moment:
          {
            _output_time.emplace_back
            (std::make_unique<MeasurerTime_position_second_moment<Subject, Geometry>>
             (subject, geometry,
              directories, identifier,
              measurement.precision));
            break;
          }
          case Measure::Type::position_variance:
          {
            _output_time.emplace_back
            (std::make_unique<MeasurerTime_position_variance<Subject, Geometry>>
             (subject, geometry,
              directories, identifier,
              measurement.precision));
            break;
          }
          case Measure::Type::mass:
          {
            _output_time.emplace_back
            (std::make_unique<MeasurerTime_mass<Subject, Geometry>>
             (subject, geometry,
              directories, identifier,
              measurement.precision));
            break;
          }
          case Measure::Type::mass_in_regions:
          {
            if constexpr (!std::is_same_v<Mask, useful::Empty>)
              _output_time.emplace_back
              (std::make_unique<MeasurerTime_mass_in_regions<Subject, Geometry, Mask>>
               (subject, geometry,
                directories, identifier,
                masks, thresholds,
                measurement.precision));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + measurement.name + ": "
                "Region masks not specified" };
            break;
          }
          case Measure::Type::velocity:
          {
            if constexpr (!std::is_same_v<VelocityField, useful::Empty>)
              _output_time.emplace_back
              (std::make_unique<MeasurerTime_vector_field<Subject, Geometry, VelocityField const&>>
               (subject, velocity_field, geometry,
                directories, identifier,
                "U",
                measurement.precision));
            else
              throw std::runtime_error{
                "Measurement type " + Measure::type(measurement.name) + ": "
                + "Velocity field not provided" };
            break;
          }
          case Measure::Type::velocity_gradient:
          {
            if constexpr (!std::is_same_v<VelocityField, useful::Empty>)
              _output_time.emplace_back
              (std::make_unique<MeasurerTime_tensor_field<Subject, Geometry>>
               (subject,
                Foam::fvc::grad(velocity_field.field()),
                geometry,
                directories, identifier,
                "gradU",
                measurement.precision));
            else
              throw std::runtime_error{
                "Measurement type " + measurement.name + ": "
                + "Velocity field not provided" };
            break;
          }
          case Measure::Type::scalar_field:
          {
            _output_time.emplace_back
            (std::make_unique<MeasurerTime_scalar_field<Subject, Geometry>>
             (subject, geometry,
              directories, identifier,
              measurement.field_name,
              measurement.precision));
            break;
          }
          case Measure::Type::vector_field:
          {
            _output_time.emplace_back
            (std::make_unique<MeasurerTime_vector_field<Subject, Geometry>>
             (subject, geometry,
              directories, identifier,
              measurement.field_name,
              measurement.precision));
            break;
          }
          case Measure::Type::tensor_field:
          {
            _output_time.emplace_back
            (std::make_unique<MeasurerTime_tensor_field<Subject, Geometry>>
             (subject, geometry,
              directories, identifier,
              measurement.field_name,
              measurement.precision));
            break;
          }
          case Measure::Type::position_periodic:
          {
            if constexpr (meta::has_periodicity_v<typename Subject::Particle::State>)
              _output_time.emplace_back
              (std::make_unique<MeasurerTime_position_periodic<Subject, Geometry>>
               (subject, geometry,
                directories, identifier,
                measurement.precision));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + measurement.name + ": "
                "Particle state must define periodicity" };
            break;
          }
          case Measure::Type::position_in_regions_periodic:
          {
            if constexpr (meta::has_periodicity_v<typename Subject::Particle::State>
                          && !std::is_same_v<Mask, useful::Empty>)
              _output_time.emplace_back
              (std::make_unique<MeasurerTime_position_in_regions_periodic<Subject, Geometry, Mask>>
               (subject, geometry,
                directories, identifier,
                masks, thresholds,
                measurement.precision));
            else
            {
              if constexpr (!meta::has_periodicity_v<typename Subject::Particle::State>)
                throw std::runtime_error{
                  std::string("Measurement type ") + measurement.name + ": "
                  "Particle state must define periodicity" };
              throw std::runtime_error{
                std::string("Measurement type ") + measurement.name + ": "
                            "Region masks not specified" };
            }
            break;
          }
          case Measure::Type::position_mean_periodic:
          {
            if constexpr (meta::has_periodicity_v<typename Subject::Particle::State>)
              _output_time.emplace_back
              (std::make_unique<MeasurerTime_position_mean_periodic<Subject, Geometry>>
               (subject, geometry,
                directories, identifier,
                measurement.precision));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + measurement.name + ": "
                "Particle state must define periodicity" };
            break;
          }
          case Measure::Type::position_second_moment_periodic:
          {
            if constexpr (meta::has_periodicity_v<typename Subject::Particle::State>)
              _output_time.emplace_back
              (std::make_unique<MeasurerTime_position_second_moment_periodic<Subject, Geometry>>
               (subject, geometry,
                directories, identifier,
                measurement.precision));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + measurement.name + ": "
                "Particle state must define periodicity" };
            break;
          }
          case Measure::Type::position_variance_periodic:
          {
            if constexpr (meta::has_periodicity_v<typename Subject::Particle::State>)
              _output_time.emplace_back
              (std::make_unique<MeasurerTime_position_variance_periodic<Subject, Geometry>>
               (subject, geometry,
                directories, identifier,
                measurement.precision));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + measurement.name + ": "
                "Particle state must define periodicity" };
            break;
          }
          case Measure::Type::absorption_time:
          {
            _output.emplace_back
            (std::make_unique<Measurer_absorption_time<Subject, Geometry>>
             (subject, geometry,
              directories, identifier,
              measurement.precision));
            break;
          }
          default:
            throw std::runtime_error{
              std::string("Measurement type ") + measurement.name + " "
              "not supported" };
        }
      }
    }
    
    /** Set end criterion. */
    void set_end_criterion(Subject const& subject)
    {
      switch (EndCriterion::type(parameters.end_criterion))
      {
        case EndCriterion::Type::time:
        {
          _end_criterion = std::make_unique<
            Criterion_time<Subject>>(subject, parameters.end_value);
          break;
        }
        case EndCriterion::Type::time_max:
        {
          _end_criterion = std::make_unique<
            Criterion_time<Subject>>(subject, parameters.time_max);
          break;
        }
        case EndCriterion::Type::mass_below:
        {
          _end_criterion = std::make_unique<
            Criterion_mass_below<Subject>>(subject, parameters.end_value);
          break;
        }
        case EndCriterion::Type::mass_above:
        {
          _end_criterion = std::make_unique<
            Criterion_mass_above<Subject>>(subject, parameters.end_value);
          break;
        }
        case EndCriterion::Type::all_absorbed:
        {
          _end_criterion = std::make_unique<
            Criterion_all_absorbed<Subject>>(subject);
          break;
        }
        case EndCriterion::Type::one_absorbed:
        {
          _end_criterion = std::make_unique<
            Criterion_one_absorbed<Subject>>(subject);
          break;
        }
        case EndCriterion::Type::fraction_not_absorbed:
        {
          _end_criterion = std::make_unique<
            Criterion_fraction_not_absorbed<Subject>>(subject, parameters.end_value);
          break;
        }
        default:
          throw std::runtime_error{
            std::string("End criterion ")
            + parameters.end_criterion
            + " not supported" };
      }
    }
    
    /** \struct NextMeasureTime PTOF/Output.h "PTOF/Output.h"
     *  \brief Polymorphic object to handle setting up the next measurement time
     \note To be used only through derived classes. */
    struct NextMeasureTime
    {
      /** Destructor. */
      virtual ~NextMeasureTime()
      {}
      
      /** \return Next measurement time. */
      double time() const
      { return _next_time; }
      
      /** \brief Set up for next measurement.  */
      void advance()
      {
        ++_next_measure;
        set_next_time();
      }
      
    protected:
      /** Constructor. */
      NextMeasureTime
      (Parameters const& parameters)
      : _parameters{ parameters }
      , _next_time{ parameters.time_min }
      {}
      
      Parameters const& _parameters;         /**< Output parameters.         */
      std::size_t _next_measure{ 0 };        /**< Index of next measurement. */
      double _next_time;                     /**< Time of next measurement.  */
      
      /** Set the time of the next measurment. */
      virtual void set_next_time() = 0;
    };
    
    /** \struct NextMeasureTime_step PTOF/Output.h "PTOF/Output.h"
    *  \brief NextMeasureTime object for constant time step. */
    struct NextMeasureTime_step final : NextMeasureTime
    {
      NextMeasureTime_step
      (Parameters const& parameters)
      : NextMeasureTime{ parameters }
      {}
      
      void set_next_time() override
      {
        NextMeasureTime::_next_time = NextMeasureTime::_parameters.time_min
          + NextMeasureTime::_next_measure
          * NextMeasureTime::_parameters.time_increment;
      }
    };
    
    /** \struct NextMeasureTime_linear PTOF/Output.h "PTOF/Output.h"
      \brief NextMeasureTime object for linearly-spaced times. */
    struct NextMeasureTime_linear final : NextMeasureTime
    {
      NextMeasureTime_linear
      (Parameters const& parameters)
      : NextMeasureTime{ parameters }
      {}
      
      void set_next_time() override
      {
        double time = NextMeasureTime::_parameters.time_min
          + NextMeasureTime::_next_measure
          * NextMeasureTime::_parameters.time_increment;
        NextMeasureTime::_next_time > NextMeasureTime::_parameters.time_max
          ? std::numeric_limits<double>::infinity()
          : time;
      }
    };
    
    /** \struct Output_Cases::NextMeasureTime_log PTOF/Output.h "PTOF/Output.h"
      \brief NextMeasureTime object for log-spaced times. */
    struct NextMeasureTime_log final : NextMeasureTime
    {
      NextMeasureTime_log
      (Parameters const& parameters)
      : NextMeasureTime{ parameters }
      {}
      
      void set_next_time() override
      {
         double time = NextMeasureTime::_parameters.time_min
           * std::pow(NextMeasureTime::_parameters.time_increment,
                      NextMeasureTime::_next_measure);
         NextMeasureTime::_next_time > NextMeasureTime::_parameters.time_max
           ? std::numeric_limits<double>::infinity()
           : time;
      }
    };
    
    std::unique_ptr<Criterion<Subject>> _end_criterion;       /**< To check if end criterion is met. */
    std::unique_ptr<NextMeasureTime> _next_measure;           /**< Handle next measure time. */
    std::vector<std::unique_ptr<MeasurerTime<Subject, Geometry>>>
      _output_time;                                           /**< Handle each output type given time. */
    std::vector<std::unique_ptr<Measurer<Subject, Geometry>>>
      _output;                                                /**< Handle each output type given nothing. */
  };
}
        
#endif /* PTOF_OUTPUT_H */
