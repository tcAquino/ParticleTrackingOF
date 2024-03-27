/**
 \file PTOF/Output.h
 \author Tomás Aquino
 \date 07/03/2022
*/

#ifndef PTOF_OUTPUT_H
#define PTOF_OUTPUT_H

#include <algorithm>
#include <array>
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
#include <Tensor.H>
#include <zero.H>
#include "CTRW/StateGetter.h"
#include "CTRW/Meta.h"
#include "General/Useful.h"
#include "PTOF/Directories.h"
#include "PTOF/Field.h"
#include "PTOF/Meta.h"
#include "PTOF/Useful.h"

namespace ptof
{
  /** \struct EndCriterion PTOF/Output.h "PTOF/Output.h"
   * \brief Keep track of names and types of end criteria. */
  struct EndCriterion
  {
    /** \enum Type
     *  \brief Implemented types. */
    enum class Type
    {
      time,                 /**< Specified time. */
      time_max,             /**< Maximum output time. */
      mass_below,           /**< Total mass below value. */
      mass_above,           /**< Total mass above value. */
      all_absorbed,         /**< All particles absorbed. */
      one_absorbed,         /**< One particle absorbed. */
      fraction_not_absorbed /**< Particle fraction not absorbed. */
    };
    
    /** \return Type from name. */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** \return Name from type. */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** \return \c true if name exists, \c false otherwise. */
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    /** Map of names to types. */
    inline static const
    std::unordered_map<std::string, Type> name_to_type
    {
      { "time", Type::time },
      { "time_max", Type::time_max },
      { "mass_below", Type::mass_below },
      { "mass_above", Type::mass_above },
      { "all_absorbed", Type::all_absorbed },
      { "one_absorbed", Type::one_absorbed },
      { "fraction_not_absorbed", Type::fraction_not_absorbed }
    };
    
    /** Map of types to names. */
    inline static const
    std::unordered_map<Type, std::string> type_to_name
    {
      { Type::time, "time" },
      { Type::time_max, "time_max" },
      { Type::mass_below, "mass_below" },
      { Type::mass_above, "mass_above" },
      { Type::all_absorbed, "all_absorbed" },
      { Type::one_absorbed, "one_absorbed" },
      { Type::fraction_not_absorbed, "fraction_not_absorbed" }
    };
  };
  
  /** \struct MeasureSpacing PTOF/Output.h "PTOF/Output.h"
   *  \brief Keep track of names and types of measure spacing */
  struct MeasureSpacing
  {
     /** \enum Type
      *  \brief Implemented types. */
    enum class Type
    {
      linear,  /**< Fixed number of linearly-spaced measurements, starting and ending at specified  times */ 
      log,     /**< Fixed number of log-spaced measurements, starting and ending at specified times.      */ 
      step     /**< At fixed time intervals, starting at specified time.                                  */ 
    };
    
    /** \return Type from name. */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** \return Name from type. */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** \return \c true if name exists, \c false otherwise. */
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    /** Map of names to types. */
    inline static const
    std::unordered_map<std::string, Type> name_to_type
    {
      { "linear", Type::linear },
      { "log", Type::log },
      { "step", Type::step }
    };
    
    /** Map of types to names. */
    inline static const
    std::unordered_map<Type, std::string> type_to_name
    {
      { Type::linear, "linear" },
      { Type::log, "log" },
      { Type::step, "step" }
    };
  };
  
  /** \struct MeasureSpacingUnits PTOF/Output.h "PTOF/Output.h"
   * \brief Keep track of names and types of measure spacing units (scaling) */
  struct MeasureSpacingUnits
  {
    /** \enum Type
     *  \brief Implemented types. */  
    enum class Type
    {
      diffusion, /**< Units of diffusion scale */
      advection, /**< Units of advection scale */
      reaction,  /**<  Units of reaction scale*/
      arbitrary  /**< Arbitrary units*/
    };
    
    /** \return Type from name. */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** \return Name from type. */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** \return \c true if criterion name exists, \c false otherwise. */
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    /** Map of names to types. */
    inline static const
    std::unordered_map<std::string, Type>
    name_to_type
    {
      { "diffusion", Type::diffusion },
      { "advection", Type::advection },
      { "reaction", Type::reaction },
      { "arbitrary", Type::arbitrary }
    };
    
    /** Map of types to names. */
    inline static const
    std::unordered_map<Type, std::string>
    type_to_name
    {
      { Type::diffusion, "diffusion" },
      { Type::advection, "advection" },
      { Type::reaction, "reaction" },
      { Type::arbitrary, "arbitrary" }
    };
  };
  
  /** \struct Measure PTOF/Output.h "PTOF/Output.h"
   * \brief Keep track of measurement types */
  struct Measure
  {
    /** \enum Type
     *  \brief Implemented types. */  
    enum class Type
    {
      position,                         /**< Time, particle tags, positions, and masses. */
      position_in_regions,              /**< Time, particle tags, positions, and masses in regions specified by masks. */
      position_mean,                    /**< Time and mean position. */
      position_second_moment,           /**< Time and position second moment. */
      position_variance,                /**< Time and position variance. */
      mass,                             /**< Time and total mass. */
      mass_in_regions,                  /**< Time and total mass in regions speciefied by masks. */
      velocity,                         /**< Time, particle tags, and local velocities. */
      velocity_gradient,                /**< Time, particle tags, and local velocity gradients. */
      scalar_field,                     /**< Time, particle tags, and local values of scalar field. */
      vector_field,                     /**< Time, particle tags, and local values of vector field. */
      tensor_field,                     /**< Time, particle tags and local values of tensor field. */
      position_periodic,                /**< Time, particle tags, true positions accounting for periodicity, and masses. */
      position_in_regions_periodic,     /**< Time, particle tags, true positions accounting for periodicity, and masses in regions specified by masks. */
      position_mean_periodic,           /**< Time and true mean position accounting for periodicity. */
      position_second_moment_periodic,  /**< Time and true position second moment accounting for periodicity. */
      position_variance_periodic,       /**< Time and true position variance accounting for periodicity. */
      absorption_time                   /**< Particle absorption times, tags, and masses at end of dynamics. */
    };
    
    /** \return Type from name. */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** \return Name from type. */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** \return \c true if criterion name exists, \c false otherwise. */
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    /** Map of names to types. */
    inline static const
    std::unordered_map<std::string, Type>name_to_type
    {
      { "position", Type::position },
      { "position_in_regions", Type::position_in_regions },
      { "position_mean", Type::position_mean },
      { "position_second_moment", Type::position_second_moment },
      { "position_variance", Type::position_variance },
      { "mass", Type::mass },
      { "mass_in_regions", Type::mass_in_regions },
      { "velocity", Type::velocity },
      { "velocity_gradient", Type::velocity_gradient },
      { "scalar_field", Type::scalar_field },
      { "vector_field", Type::scalar_field },
      { "tensor_field", Type::scalar_field },
      { "position_periodic", Type::position_periodic },
      { "position_in_regions_periodic", Type::position_in_regions_periodic },
      { "position_mean_periodic", Type::position_mean_periodic },
      { "position_second_moment_periodic", Type::position_second_moment_periodic },
      { "position_variance_periodic", Type::position_variance_periodic },
      { "absorption_time", Type::absorption_time }
    };
    
    /** Map of types to names. */
    inline static const
    std::unordered_map<Type, std::string>type_to_name
    {
      { Type::position, "position" },
      { Type::position_in_regions, "position_in_regions" },
      { Type::position_mean, "position_mean" },
      { Type::position_second_moment, "position_second_moment" },
      { Type::position_variance, "position_variance" },
      { Type::mass, "mass" },
      { Type::mass_in_regions, "mass_in_regions" },
      { Type::velocity, "velocity" },
      { Type::velocity_gradient, "velocity_gradient" },
      { Type::position_periodic, "position_periodic" },
      { Type::position_in_regions_periodic, "position_in_regions_periodic" },
      { Type::position_mean_periodic, "position_mean_periodic" },
      { Type::position_second_moment_periodic, "position_second_moment_periodic" },
      { Type::position_variance_periodic, "position_variance_periodic" },
      { Type::absorption_time, "absorption_time" }
    };
  };
  
  /** \class Criterion PTOF/Output.h "PTOF/Output.h"
   * \brief  Polymorphic functor to implement criteria for ending measurements.
   \note To be used only through derived classes. */
  template <typename Subject>
  struct Criterion
  {
    virtual bool operator()
    (Subject const& subject, double time) const;
    
    virtual ~Criterion()
    {}
    
  protected:
    Criterion()
    {}
  };
  
  /** \class Criterion_time PTOF/Output.h "PTOF/Output.h
   *  \brief Check if time is greater than value. */
  template <typename Subject>
  struct Criterion_time final : Criterion<Subject>
  {
    Criterion_time(double end_value)
    : end_value{ end_value }
    {}
    
    bool operator()
    (Subject const& subject, double time) const override
    {
      return time >= end_value;
    }
      
    double end_value;
  };
  
  /** \class Criterion_mass_below PTOF/Output.h "PTOF/Output.h"
   *  \brief Check if mass is less than or equal to value. */
  template <typename Subject>
  struct Criterion_mass_below final : Criterion<Subject>
  {
    Criterion_mass_below(double end_value)
    : end_value{ end_value }
    {}
    
    bool operator()
    (Subject const& subject, double time) const override
    {
      return mass(subject, time) <= end_value;
    }
      
    double end_value;
  };
  
  /** \class Criterion_mass_above PTOF/Output.h "PTOF/Output.h"
   *  \brief Check if mass is greater than or equal to value. */
  template <typename Subject>
  struct Criterion_mass_above final : Criterion<Subject>
  {
    Criterion_mass_above(double end_value)
    : end_value{ end_value }
    {}
    
    bool operator()
    (Subject const& subject, double time) const override
    {
      return mass(subject, time) >= end_value;
    }
      
    double end_value;
  };
  
  /** \class Criterion_all_absorbed PTOF/Output.h "PTOF/Output.h"
   * \brief Check if all particles have been absorbed. */
  template <typename Subject>
  struct Criterion_all_absorbed final : Criterion<Subject>
  {
    Criterion_all_absorbed()
    {}
    
    bool operator()
    (Subject const& subject, double time) const override
    {
      return nr_absorbed(subject, time) == subject.size();
    }
  };
  
  /** \class Criterion_one_absorbed PTOF/Output.h "PTOF/Output.h"
   *  \brief Check if at least one particle has been absorbed. */
  template <typename Subject>
  struct Criterion_one_absorbed final : Criterion<Subject>
  {
    Criterion_one_absorbed()
    {}
    
    bool operator()(Subject const& subject, double time) const
    {
      return nr_absorbed(subject, time) > 0;
    }
  };
  
  /** \class Criterion_fraction_not_absorbed PTOF/Output.h "PTOF/Output.h"
   *  \brief Check if the fraction of particles that have not been absorbed is at most a certain value. */
  template <typename Subject>
  struct Criterion_fraction_not_absorbed final : Criterion<Subject>
  {
    Criterion_fraction_not_absorbed(double end_value)
    : end_value{ end_value }
    {}
    
    bool operator()
    (Subject const& subject, double time) const override
    {
      return subject.size() - nr_absorbed(subject, time) <= std::size_t(end_value*subject.size());
    }
      
    double end_value;
  };
  
  /** \struct OutputParameters_Cases PTOF/Output.h "PTOF/Output.h"
   *  \brief Output parameters to handle the output types in Output_Cases. */
  struct OutputParameters_Cases
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
      Measurement(std::string name, std::string parameters = {})
      : name{ name }
      , parameters{ parameters }
      {}
      
      std::string name;
      std::string parameters;
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
    OutputParameters_Cases
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
      read_output_types(input);
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
        "- Measurement types (any number):\n"
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
    void read_output_types
    (IStream& input)
    {
      std::string line;
      while (std::getline(input, line))
      {
        auto split_line = useful::split(useful::remove_carriage_return(line));
        if (split_line.empty())
          throw std::runtime_error{
            "Could not parse measurement types" };
        std::string const& name = split_line[0];
        measurements.emplace_back(name);
        if (name == "scalar_field" ||
            name == "vector_field" ||
            name == "tensor_field")
        {
          if (split_line.size() != 2)
            throw std::runtime_error{
              std::string("Measurement type ") + name + ": "
                "Field name not provided" };
          measurements.back().parameters = split_line[1];
        }
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
    using Parameters = OutputParameters_Cases;    /**< Output parameters. */
    /** Interpolation of scalar fields for outputs that need it. */
    using InterpolatedScalarField = ptof::ScalarField_LinearInterpolation_OF
      <Foam::volScalarField, typename Geometry::Locator const&, CheckOptions::Check>;
    /** Interpolation of vector fields for outputs that need it. */
    using InterpolatedVectorField = ptof::VectorField_LinearInterpolation_OF
      <Foam::volVectorField, typename Geometry::Locator const&, CheckOptions::Check>;
    
    /** Constructor.
     \param subject CTRW object to measure.
     \param velocity_field Velocity field as a function of state.
     \param geometry Domain geometry info and utilities.
     \param directories Current case directory information.
     \param parameters Output parameters.
     \param identifier String to include in names of output files.
     \param masks Container of mask reference wrappers. Masks are scalar fields assigned values to mesh cells through operator[].
     \param thresholds Vector of thresholds for each mask, such that cells where a mask is above the threshold are considered.
     \param precision Number of digits after decimal point in output, in scientific notation.
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
     std::vector<double> thresholds = {},
     int precision = 8)
    : parameters{ parameters }
    , _subject{ subject }
    , _geometry{ geometry }
    {
      set_measure_types(directories, identifier,
                        velocity_field,
                        masks, thresholds,
                        precision);
      set_end_criterion();
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
     std::initializer_list<double> thresholds = {},
     int precision = 8)
    : Output_Cases(subject, velocity_field, geometry, directories,
                   parameters, identifier,
                   std::vector<std::reference_wrapper<const Mask>>(masks),
                   std::vector<double>(thresholds),
                   precision)
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
     std::initializer_list<double> thresholds = {},
     int precision = 8)
    : Output_Cases(subject, velocity_field, geometry, directories,
                   parameters, identifier,
                   masks,
                   std::vector<double>(thresholds),
                   precision)
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
     std::vector<double> thresholds = {},
     int precision = 8)
    : parameters{ parameters }
    , _subject{ subject }
    , _geometry{ geometry }
    {
      set_measure_types(directories, identifier,
                        useful::Empty{},
                        masks, thresholds,
                        precision);
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
     std::initializer_list<double> thresholds = {},
     int precision = 8)
    : Output_Cases(subject, geometry, directories,
                   parameters, identifier,
                   std::vector<std::reference_wrapper<const Mask>>(masks),
                   std::vector<double>(thresholds),
                   precision)
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
     std::initializer_list<double> thresholds = {},
     int precision = 8)
    : Output_Cases(subject, geometry, directories,
                   parameters, identifier,
                   masks,
                   std::vector<double>(thresholds),
                   precision)
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
     std::vector<double> thresholds = {},
     int precision = 8)
    : Output_Cases(subject, geometry, directories,
                   parameters, identifier,
                   std::vector<std::reference_wrapper<const Mask>>(masks),
                   thresholds,
                   precision)
    {}
    
    /** \return \c true  if end simulation criterion is satisfied, \c false otherwise. */
    bool done(double time) const
    {
      return _end_criterion->operator()(_subject, time);
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
    {
      return _next_measure->time();
    }
        
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
      output << "\n";
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
                                             { return aa.parameters.length() < bb.parameters.length(); })
                            ->parameters.length()) + 1;
      for (auto const& measurement : parameters.measurements)
      {
        output << prefix
               << std::left << std::setw(width_name) << measurement.name;
        if (!measurement.parameters.empty())
          output << std::left << std::setw(width_field)
                 << measurement.parameters;
        output << "\n";
      }
    }
    
    /** \brief Set up output streams for requested output types.
     \param directories Current case directory information.
     \param identifier String to include in names of output files.
     \param masks Container of mask reference wrappers. Masks are scalar fields assigned values to mesh cells through operator[].
     \param thresholds Vector of thresholds for each mask, such that cells where a mask is above the threshold are considered.
     \param precision Number of digits after decimal point in output, in scientific notation. */
    template <typename VelocityField = useful::Empty, typename Mask = useful::Empty>
    void set_measure_types
    (Directories const& directories,
     std::string const& identifier,
     VelocityField const& velocity_field = {},
     std::vector<std::reference_wrapper<const Mask>> masks = {},
     std::vector<double> thresholds = {},
     int precision = 8)
    {
      for (auto const& measurement : parameters.measurements)
      {
        switch (Measure::type(measurement.name))
        {
          case Measure::Type::position:
          {
            _output_time.emplace_back
            (std::make_unique<Output_position>
             (_subject, _geometry,
               directories, identifier, parameters,
               precision));
            break;
          }
          case Measure::Type::position_in_regions:
          {
            if constexpr (!std::is_same_v<Mask, useful::Empty>)
              _output_time.emplace_back
              (std::make_unique<Output_position_in_regions<Mask>>
               (_subject, _geometry,
                 directories, identifier, parameters,
                 masks, thresholds,
                 precision));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + measurement.name + ": "
                "Region masks not specified" };
            break;
          }
          case Measure::Type::position_mean:
          {
            _output_time.emplace_back
            (std::make_unique<Output_position_mean>
             (_subject, _geometry,
              directories, identifier, parameters,
              precision));
            break;
          }
          case Measure::Type::position_second_moment:
          {
            _output_time.emplace_back
            (std::make_unique<Output_position_second_moment>
             (_subject, _geometry,
              directories, identifier, parameters,
              precision));
            break;
          }
          case Measure::Type::position_variance:
          {
            _output_time.emplace_back
            (std::make_unique<Output_position_variance>
             (_subject, _geometry,
              directories, identifier, parameters,
              precision));
            break;
          }
          case Measure::Type::mass:
          {
            _output_time.emplace_back
            (std::make_unique<Output_mass>
             (_subject, _geometry,
              directories, identifier, parameters,
              precision));
            break;
          }
          case Measure::Type::mass_in_regions:
          {
            if constexpr (!std::is_same_v<Mask, useful::Empty>)
              _output_time.emplace_back
              (std::make_unique<Output_mass_in_regions<Mask>>
               (_subject, _geometry,
                directories, identifier, parameters,
                masks, thresholds,
                precision));
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
              (std::make_unique<Output_vector_field<VelocityField const&>>
               (_subject, velocity_field, _geometry,
                directories, identifier, parameters,
                "U",
                precision));
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
              (std::make_unique<Output_tensor_field<>>
               (_subject,
                Foam::fvc::grad(velocity_field.field()),
                _geometry,
                directories, identifier, parameters,
                "gradU",
                precision));
            else
              throw std::runtime_error{
                "Measurement type " + measurement.name + ": "
                + "Velocity field not provided" };
            break;
          }
          case Measure::Type::scalar_field:
          {
            _output_time.emplace_back
            (std::make_unique<Output_scalar_field<>>
             (_subject, _geometry,
              directories, identifier, parameters,
              measurement.parameters,
              precision));
            break;
          }
          case Measure::Type::vector_field:
          {
            _output_time.emplace_back
            (std::make_unique<Output_vector_field<>>
             (_subject, _geometry,
              directories, identifier, parameters,
              measurement.parameters,
              precision));
            break;
          }
          case Measure::Type::tensor_field:
          {
            _output_time.emplace_back
            (std::make_unique<Output_tensor_field<>>
             (_subject, _geometry,
              directories, identifier, parameters,
              measurement.parameters,
              precision));
            break;
          }
          case Measure::Type::position_periodic:
          {
            if constexpr (meta::has_periodicity_v<typename Subject::Particle::State>)
              _output_time.emplace_back
              (std::make_unique<Output_position_periodic>
               (_subject, _geometry,
                directories, identifier, parameters,
                precision));
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
              (std::make_unique<Output_position_in_regions_periodic<Mask>>
               (_subject, _geometry,
                 directories, identifier, parameters,
                 masks, thresholds,
                 precision));
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
              (std::make_unique<Output_position_mean_periodic>
               (_subject, _geometry,
                directories, identifier, parameters,
                precision));
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
              (std::make_unique<Output_position_second_moment_periodic>
               (_subject, _geometry,
                directories, identifier, parameters,
                precision));
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
              (std::make_unique<Output_position_variance_periodic>
               (_subject, _geometry,
                directories, identifier, parameters,
                precision));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + measurement.name + ": "
                "Particle state must define periodicity" };
            break;
          }
          case Measure::Type::absorption_time:
          {
            _output.emplace_back
            (std::make_unique<Output_absorption_time>
             (_subject, _geometry,
              directories, identifier, parameters,
              precision));
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
    void set_end_criterion()
    {
      switch (EndCriterion::type(parameters.end_criterion))
      {
        case EndCriterion::Type::time:
        {
          _end_criterion = std::make_unique<
            Criterion_time<Subject>>(parameters.end_value);
          break;
        }
        case EndCriterion::Type::time_max:
        {
          _end_criterion = std::make_unique<
            Criterion_time<Subject>>(parameters.time_max);
          break;
        }
        case EndCriterion::Type::mass_below:
        {
          _end_criterion = std::make_unique<
            Criterion_mass_below<Subject>>(parameters.end_value);
          break;
        }
        case EndCriterion::Type::mass_above:
        {
          _end_criterion = std::make_unique<
            Criterion_mass_above<Subject>>(parameters.end_value);
          break;
        }
        case EndCriterion::Type::all_absorbed:
        {
          _end_criterion = std::make_unique<
            Criterion_all_absorbed<Subject>>();
          break;
        }
        case EndCriterion::Type::one_absorbed:
        {
          _end_criterion = std::make_unique<
            Criterion_one_absorbed<Subject>>();
          break;
        }
        case EndCriterion::Type::fraction_not_absorbed:
        {
          _end_criterion = std::make_unique<
            Criterion_fraction_not_absorbed<Subject>>(parameters.end_value);
          break;
        }
        default:
          throw std::runtime_error{
            std::string("End criterion ")
            + parameters.end_criterion
            + " not supported" };
      }
    }
    
    /** \struct Output_Cases::NextMeasureTime PTOF/Output.h "PTOF/Output.h"
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
    
    /** \struct Output_Cases::NextMeasureTime_step PTOF/Output.h "PTOF/Output.h"
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
    
    /** \struct Output_Cases::OutputTime PTOF/Output.h "PTOF/Output.h"
     \brief Polymorphic output handler for outputting time and some quantity.
     \note To be used only through derived classes. */
    struct OutputTime
    {
      /** Destructor. */
      virtual ~OutputTime()
      { _output.close(); }
      
      /** \brief Make measurement and output, given current time \c time. */
      virtual void operator()(double time) = 0;
      
      /** \brief Update internal state. */
      virtual void update(double time, double time_of_change)
      {}
      
    protected:
      /** \brief Set up output streams for requested output types.
      \param subject CTRW object to measure.
      \param geometry Domain geometry info and utilities.
      \param directories Current case directory information.
      \param output_name Name of output type.
      \param identifier String to include in names of output files.
      \param parameters Output parameters.
      \param precision Number of digits after decimal point in output, in scientific notation. */
      OutputTime
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : _subject{ subject }
      , _geometry{ geometry }
      , _locator{ geometry.locator }
      , _output{ open_write(directories, output_name, identifier) }
      {
        _output << std::setprecision(precision)
                << std::scientific;
      }
      
      Subject const& _subject;                      /**< CTRW object to measure. */
      Geometry const& _geometry;                    /**< Domain geometry info and utilities. */
      typename Geometry::Locator const& _locator;   /**< Object to locate positions in mesh. */
      std::ofstream _output;                        /**< Output stream. */
        
      /** \brief Open output file for a given output type.
      \param directories Current case directory information.
      \param output_name Name of output type.
      \param identifier String to include in names of output files.
      \return Output stream. */
      std::ofstream open_write
      (Directories const& directories,
       std::string const& output_name,
       std::string const& identifier)
      {
        return
          useful::open_write(directories.dir_output
          + "/Data_" + output_name
          + (identifier.empty() ? "" : "_" + identifier)
          + ".dat") ;
      }
    };
        
    /** \struct Output_Cases::Output PTOF/Output.h "PTOF/Output.h"
     \brief Polymorphic output handler for outputting some quantity.
     \note To be used only through derived classes. */
    struct Output
    {
      virtual ~Output()
      { _output.close(); }
      
      /** \brief Make measurement and output. */
      virtual void operator()() = 0;
      
    protected:
      /** Constructor.
      \param subject CTRW object to measure.
      \param geometry Domain geometry info and utilities.
      \param directories Current case directory information.
      \param output_name Name of output type.
      \param identifier String to include in names of output files.
      \param parameters Output parameters.
      \param precision Number of digits after decimal point in output, in scientific notation. */
      Output
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : _subject{ subject }
      , _geometry{ geometry }
      , _locator{ geometry.locator }
      , _output{ open_write(directories, output_name, identifier) }
      {
        _output << std::setprecision(precision)
                << std::scientific;
      }
      
      Subject const& _subject;                      /**< CTRW object to measure. */
      Geometry const& _geometry;                    /**< Domain geometry info and utilities. */
      typename Geometry::Locator const& _locator;   /**< Object to locate positions in mesh. */
      std::ofstream _output;                        /**< Output stream. */
      
      /**
      \param directories Current case directory information.
      \param output_name Name of output type.
      \param identifier String to include in names of output files.
      \return Output stream. */
      std::ofstream open_write
      (Directories const& directories,
       std::string const& output_name,
       std::string const& identifier)
      {
        return
          useful::open_write(directories.dir_output
          + "/Data_" + output_name
          + (identifier.empty() ? "" : "_" + identifier)
          + ".dat");
      }
    };
    
    /** \struct Output_Cases::Output_position PTOF/Output.h "PTOF/Output.h"
     *  \brief Output time, tags, positions, and masses. */
    struct Output_position final : OutputTime
    {
      Output_position
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "position", identifier, parameters,
          precision }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(12, int(1 + std::string{ "Tag" }.length())),
        std::max(9 + precision, int(2 + std::string{ "Position_" }.length())),
        std::max(9 + precision, int(1 + std::string{ "Mass" }.length())) }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time"
                            << std::setw(_column_widths[1]) << "Tag";
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
          OutputTime::_output << std::setw(_column_widths[2])
                              << "Position_" + std::to_string(dd);
        OutputTime::_output << std::setw(_column_widths[3]) << "Mass";
          for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
        OutputTime::_output << std::setw(_column_widths[1]) << "Tag"
                            << std::setw(_column_widths[2]) << "..."
                            << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed
              && part.state_old().time <= time)
          {
            OutputTime::_output << std::setw(_column_widths[1]) << state.tag;
            useful::print(OutputTime::_output, state.position, _column_widths[2]);
            OutputTime::_output << _column_widths[3] << state.mass;
          }
        }
        OutputTime::_output << "\n";
      }
      
    private:
      std::array<int, 4> _column_widths;
    };
        
    /** \class Output_Cases::Output_position_in_regions PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time, tags, positions, and masses within regions specified by masks. */
    template <typename Mask>
    struct Output_position_in_regions final : OutputTime
    {
      Output_position_in_regions
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       std::vector<std::reference_wrapper<const Mask>> masks,
       std::vector<double> thresholds = {},
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "position_in_regions", identifier, parameters,
          precision }
      , _masks{ masks }
      , _thresholds{ thresholds }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(12, int(1 + std::string{ "Tag" }.length())),
        std::max(9 + precision, int(2 + std::string{ "Position_" }.length())),
        std::max(9 + precision, int(1 + std::string{ "Mass" }.length())),
        std::max(9, int(4 + std::string{ "In_region_" }.length())) }
      {
        _thresholds.resize(masks.size(), 0.);
        OutputTime::_output << std::setw(_column_widths[0]) << "Time"
                            << std::setw(_column_widths[1]) << "Tag";
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
          OutputTime::_output << std::setw(_column_widths[2])
                              << "Position_" + std::to_string(dd);
        OutputTime::_output << std::setw(_column_widths[3]) << "Mass";
        for (std::size_t ii = 0; ii < _masks.size(); ++ii)
          OutputTime::_output << std::setw(_column_widths[4])
                              << "In_region_" + std::to_string(ii);
        OutputTime::_output << std::setw(_column_widths[1]) << "Tag"
                            << std::setw(_column_widths[2]) << "..."
                            << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          auto cell = OutputTime::_locator(state);
          OutputTime::_output << std::setw(_column_widths[1]) << state.tag;
          if (!state.info.absorbed
              && part.state_old().time <= time
              && !outside(cell))
          {
            std::vector<int> in_region(_masks.size(), 0);
            for (std::size_t ii = 0; ii < _masks.size(); ++ii)
              if (_masks[ii].get()[cell] > _thresholds[ii])
                in_region[ii] = 1;
            if (std::any_of(in_region.begin(), in_region.end(),
                            [](int ii){ return ii > 0; }))
            {
              useful::print(OutputTime::_output, state.position, _column_widths[2]);
              OutputTime::_output << std::setw(_column_widths[3]) << state.mass;
              useful::print(OutputTime::_output, in_region, _column_widths[4]);
            }
          }
        }
        OutputTime::_output << "\n";
      }
      
    private:
      std::vector<std::reference_wrapper<const Mask>> _masks;
      std::vector<double> _thresholds;
      std::array<int, 5> _column_widths;
    };
    
    /** \struct Output_Cases::Output_position_mean PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and mean position (weighted by mass). */
    struct Output_position_mean final : OutputTime
    {
      Output_position_mean
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "position_mean", identifier, parameters,
          precision }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(12, int(2 + std::string{ "Position_mean_" }.length())) }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time";
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
          OutputTime::_output << std::setw(_column_widths[1])
                              << "Position_mean_" + std::to_string(dd);
        OutputTime::_output << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        useful::print(OutputTime::_output,
                      position_mean(OutputTime::_subject, time),
                      _column_widths[1]);
        OutputTime::_output << "\n";
      }
      
    private:
      std::array<int, 2> _column_widths;
    };
        
    /** \struct Output_Cases::Output_position_second_moment PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and second moment of position (weighted by mass). */
    struct Output_position_second_moment final : OutputTime
    {
      Output_position_second_moment
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "position_second_moment", identifier, parameters,
          precision }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(9 + precision, int(2 + std::string{ "Position_second_moment_" }.length())) }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time";
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
          OutputTime::_output << std::setw(_column_widths[1])
                              << "Position_second_moment_" + std::to_string(dd);
        OutputTime::_output << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        useful::print(OutputTime::_output,
                      position_second_moment(OutputTime::_subject, time),
                      _column_widths[1]);
        OutputTime::_output << "\n";
      }
      
    private:
      std::array<int, 2> _column_widths;
    };
        
    /** \struct Output_Cases::Output_position_variance PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and position variance (weighted by mass). */
    struct Output_position_variance final : OutputTime
    {
      Output_position_variance
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "position_variance", identifier, parameters,
          precision }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(9 + precision, int(2 + std::string{ "Position_variance_" }.length())) }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time";
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
          OutputTime::_output << std::setw(_column_widths[1])
                              << "Position_variance_" + std::to_string(dd);
        OutputTime::_output << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        useful::print(OutputTime::_output,
                      position_variance(OutputTime::_subject, time),
                      _column_widths[1]);
        OutputTime::_output << "\n";
      }
      
    private:
      std::array<int, 2> _column_widths;
    };
        
    /** \struct Output_Cases::Output_mass PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and total mass. */
    struct Output_mass final : OutputTime
    {
      Output_mass
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "mass", identifier, parameters,
          precision }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(9 + precision, int(1 + std::string{ "Mass" }.length())) }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time"
                            << std::setw(_column_widths[1]) << "Mass"
                            << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time
                            << std::setw(_column_widths[1]) << mass(OutputTime::_subject, time)
                            << "\n";
      }
      
    private:
      std::array<int, 2> _column_widths;
    };
        
    /** \class Output_Cases::Output_mass_in_regions PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and total mass in regions specified by maks. */
    template <typename Mask>
    struct Output_mass_in_regions final : OutputTime
    {
      Output_mass_in_regions
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       std::vector<std::reference_wrapper<const Mask>> masks,
       std::vector<double> thresholds = {},
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "mass_in_regions", identifier, parameters,
          precision }
      , _masks{ masks }
      , _thresholds{ thresholds }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(9 + precision, int(4 + std::string{ "Mass_region_" }.length())) }
      {
        _thresholds.resize(masks.size(), 0.);
        OutputTime::_output << std::setw(_column_widths[0]) << "Time";
        for (std::size_t ii = 0; ii < masks.size(); ++ii)
          OutputTime::_output << std::setw(_column_widths[4])
                              << "Mass_region_" + std::to_string(ii);
        OutputTime::_output << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        useful::print(OutputTime::_output,
                      mass(OutputTime::_subject, time,
                           OutputTime::_locator,
                           _masks, _thresholds),
                      _column_widths[1]);
        OutputTime::_output << "\n";
      }
      
    private:
      std::vector<std::reference_wrapper<const Mask>> _masks;
      std::vector<double> _thresholds;
      std::array<int, 2> _column_widths;
    };
        
        /**
     \struct Output_Cases::Output_scalar_field PTOF/Output.h "PTOF/Output.h"
     \brief Output time, tags, and scalar field values.
    */
    template <typename Field = InterpolatedScalarField>
    struct Output_scalar_field final : OutputTime
    {
      Output_scalar_field
      (Subject const& subject,
       Field&& field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       std::string const& field_name,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          std::string{ "scalar_field_" } + field_name,
          identifier, parameters,
          precision }
      , _field_name{ field_name }
      , _field{ std::forward<Field>(field) }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(12, int(1 + std::string{ "Tag" }.length())),
        std::max(9 + precision, int(2 + (field_name + "_").length())) }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time"
                            << std::setw(_column_widths[1]) << "Tag"
                            << std::setw(_column_widths[2]) << field_name
                            << std::setw(_column_widths[1]) << "Tag"
                            << std::setw(_column_widths[2]) << "..."
                            << "\n";
      }
      
      Output_scalar_field
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       std::string const& field_name,
       int precision = 8)
      : Output_scalar_field{ subject,
          typename Output_Cases::InterpolatedScalarField{
            { Foam::IOobject{ field_name,
              geometry.mesh().time().timeName(),
              geometry.mesh(), Foam::IOobject::MUST_READ,
              Foam::IOobject::NO_WRITE },
            geometry.mesh() },
            geometry.locator },
          geometry, directories,
          identifier, parameters,
          field_name,
          precision }
      {}
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed
              && part.state_old().time <= time)
            OutputTime::_output << std::setw(_column_widths[1]) << state.tag
                                << std::setw(_column_widths[2]) << _field(state);
        }
        OutputTime::_output << "\n";
      }
      
      void update(double time, double time_of_change) override
      {
        if constexpr (!std::is_reference_v<Field>)
          if (time >= time_of_change)
            _field.set({
              Foam::IOobject{ _field_name,
                OutputTime::_geometry.mesh().time().timeName(),
                OutputTime::_geometry.mesh(), Foam::IOobject::MUST_READ,
                Foam::IOobject::NO_WRITE },
              OutputTime::_geometry.mesh() });
      }
      
    private:
      std::string _field_name;
      Field _field;
      std::array<int, 3> _column_widths;
    };
    template <typename Field>
    Output_scalar_field
    (Subject const&, Field&&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&,
     int) ->
    Output_scalar_field<Field>;
    template <typename Field>
    Output_scalar_field
    (Subject const&, Field&&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&) ->
    Output_scalar_field<Field>;
    Output_scalar_field
    (Subject const&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&,
     int) ->
    Output_scalar_field<typename Output_Cases::InterpolatedScalarField>;
    Output_scalar_field
    (Subject const&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&) ->
    Output_scalar_field<typename Output_Cases::InterpolatedScalarField>;
        
     /** \class Output_Cases::Output_vector_field PTOF/Output.h "PTOF/Output.h"
      *  \brief  Output time, tags, and vector field values. */
    template <typename Field = InterpolatedVectorField>
    struct Output_vector_field final : OutputTime
    {
      Output_vector_field
      (Subject const& subject,
       Field&& field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       std::string const& field_name,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          std::string{ "vector_field_" } + field_name,
          identifier, parameters,
          precision }
      , _field{ std::forward<Field>(field) }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(12, int(1 + std::string{ "Tag" }.length())),
        std::max(9 + precision, int(2 + (field_name + "_").length())) }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time"
                            << std::setw(_column_widths[1]) << "Tag";
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
          OutputTime::_output << std::setw(_column_widths[2])
                              << field_name + "_" + std::to_string(dd);
        OutputTime::_output << std::setw(_column_widths[1]) << "Tag"
                            << std::setw(_column_widths[2]) << "..."
                            << "\n";
      }
      
      Output_vector_field
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       std::string const& field_name,
       int precision = 8)
      : Output_vector_field{ subject,
          typename Output_Cases::InterpolatedVectorField{
            { Foam::IOobject{ field_name,
              geometry.mesh().time().timeName(),
              geometry.mesh(), Foam::IOobject::MUST_READ,
              Foam::IOobject::NO_WRITE },
            geometry.mesh() },
            geometry.locator },
          geometry, directories,
          identifier, parameters,
          field_name,
          precision }
      {}

      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed
              && part.state_old().time <= time)
          {
            OutputTime::_output << std::setw(_column_widths[1]) << state.tag;
            useful::print(OutputTime::_output,
                          _field(state),
                          _column_widths[2]);
          }
        }
        OutputTime::_output << "\n";
      }
      
      void update(double time, double time_of_change) override
      {
        if constexpr (!std::is_reference_v<Field>)
          if (time >= time_of_change)
            _field.set({
              Foam::IOobject{ _field_name,
                OutputTime::_geometry.mesh().time().timeName(),
                OutputTime::_geometry.mesh(), Foam::IOobject::MUST_READ,
                Foam::IOobject::NO_WRITE },
              OutputTime::_geometry.mesh() });
      }
       
     private:
      std::string _field_name;
      Field _field;
      std::array<int, 3> _column_widths;
    };
    template <typename Field>
    Output_vector_field
    (Subject const&, Field&&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&,
     int) ->
    Output_vector_field<Field>;
    template <typename Field>
    Output_vector_field
    (Subject const&, Field&&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&) ->
    Output_vector_field<Field>;
    Output_vector_field
    (Subject const&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&,
     int) ->
    Output_vector_field<typename Output_Cases::InterpolatedVectorField>;
    Output_vector_field
    (Subject const&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&) ->
    Output_vector_field<typename Output_Cases::InterpolatedVectorField>;
        
     /** \struct Output_Cases::Output_tensor_field PTOF/Output.h "PTOF/Output.h"
      *  \brief  Output time, tag, and tensor field values. */
    template <typename Field = Foam::volTensorField>
    struct Output_tensor_field final : OutputTime
    {
      Output_tensor_field
      (Subject const& subject,
       Field&& field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       std::string const& field_name,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          std::string{ "tensor_field_" } + field_name,
          identifier, parameters,
          precision }
      , _field{ std::forward<Field>(field) }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(12, int(1 + std::string{ "Tag" }.length())),
        std::max(9 + precision, int(3 + (field_name + "_").length())) }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time"
                            << std::setw(_column_widths[1]) << "Tag";
        for (std::size_t dd1 = 0; dd1 < Geometry::dim; ++dd1)
          for (std::size_t dd2 = 0; dd2 < Geometry::dim; ++dd2)
            OutputTime::_output << std::setw(_column_widths[2])
                                << field_name + "_" + std::to_string(dd1) + std::to_string(dd2);
        OutputTime::_output << std::setw(_column_widths[1]) << "Tag"
                            << std::setw(_column_widths[2]) << "..."
                            << "\n";
      }
      
      Output_tensor_field
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       std::string const& field_name,
       int precision = 8)
      : Output_tensor_field{ subject,
          Foam::volTensorField{ Foam::IOobject{ field_name,
            geometry.mesh().time().timeName(),
            geometry.mesh(), Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE },
            geometry.mesh() },
          geometry, directories,
          identifier, parameters,
          field_name,
          precision }
      {}

      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed
              && part.state_old().time <= time)
          {
            OutputTime::_output << std::setw(_column_widths[1]) << state.tag;
            auto cell = OutputTime::_locator(state);
            if (outside(cell))
              for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
                useful::print(OutputTime::_output,
                              std::vector<double>(Geometry::dim, 0.),
                              _column_widths[2]);
            else
              for (std::size_t dd1 = 0; dd1 < Geometry::dim; ++dd1)
                for (std::size_t dd2 = 0; dd2 < Geometry::dim; ++dd2)
                  OutputTime::_output << std::setw(_column_widths[2])
                                      << _field[cell].row(dd1)[dd2];
          }
        }
        OutputTime::_output << "\n";
      }

      void update(double time, double time_of_change) override
      {
        if constexpr (!std::is_reference_v<Field>)
          if (time >= time_of_change)
            _field = {
              Foam::IOobject{ _field_name,
                OutputTime::_geometry.mesh().time().timeName(),
                OutputTime::_geometry.mesh(), Foam::IOobject::MUST_READ,
                Foam::IOobject::NO_WRITE },
              OutputTime::_geometry.mesh() };
      }
    
    private:
      std::string _field_name;
      Field _field;
      std::array<int, 3> _column_widths;
    };
    template <typename Field>
    Output_tensor_field
    (Subject const&, Field&&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&,
     int) ->
    Output_tensor_field<Field>;
    template <typename Field>
    Output_tensor_field
    (Subject const&, Field&&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&) ->
    Output_tensor_field<Field>;
    Output_tensor_field
    (Subject const&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&,
     int) ->
    Output_tensor_field<Foam::volTensorField>;
    Output_tensor_field
    (Subject const&,
     Geometry const&,
     Directories const&,
     std::string const&,
     Parameters const&,
     std::string const&) ->
    Output_tensor_field<Foam::volTensorField>;
        
    /** \class Output_Cases::Output_position_periodic PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time, tags, positions, and masses, with positions accounting for periodicity. */
    struct Output_position_periodic final : OutputTime
    {
      Output_position_periodic
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "position_periodic", identifier, parameters,
          precision }
      , _getter_position{ geometry.boundary_periodic }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(12, int(1 + std::string{ "Tag" }.length())),
        std::max(9 + precision, int(3 + std::string{ "Position_" }.length())),
        std::max(9 + precision, int(1 + std::string{ "Mass" }.length()))
      }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time"
                            << std::setw(_column_widths[1]) << "Tag";
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
          OutputTime::_output << std::setw(_column_widths[2])
                              << "Position_" + std::to_string(dd);
        OutputTime::_output << std::setw(_column_widths[3]) << "Mass"
                            << std::setw(_column_widths[1]) << "Tag"
                            << std::setw(_column_widths[2]) << "..."
                            << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed
              && part.state_old().time <= time)
          {
            OutputTime::_output << std::setw(_column_widths[1]) << state.tag;
            useful::print(OutputTime::_output, _getter_position(state),
                          _column_widths[2]);
            OutputTime::_output << std::setw(_column_widths[3]) << state.mass;
          }
        }
        OutputTime::_output << "\n";
      }
      
    private:
      using Boundary = std::decay_t<decltype(std::declval<Geometry>().boundary_periodic)>;
      ctrw::Get_position_periodic<Boundary const&> _getter_position;
      std::array<int, 4> _column_widths;
    };
        
    /** \class Output_Cases::Output_position_in_regions_periodic PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time, tags, positions, and masses within regions specified by masks. */
    template <typename Mask>
    struct Output_position_in_regions_periodic final : OutputTime
    {
      Output_position_in_regions_periodic
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       std::vector<std::reference_wrapper<const Mask>> masks,
       std::vector<double> thresholds = {},
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "position_in_regions_periodic", identifier, parameters,
          precision }
      , _getter_position{ geometry.boundary_periodic }
      , _masks{ masks }
      , _thresholds{ thresholds }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(12, int(1 + std::string{ "Tag" }.length())),
        std::max(9 + precision, int(2 + std::string{ "Position_" }.length())),
        std::max(9 + precision, int(1 + std::string{ "Mass" }.length())),
        std::max(2, int(4 + std::string{ "In_region_" }.length())) }
      {
        _thresholds.resize(masks.size(), 0.);
        OutputTime::_output << std::setw(_column_widths[0]) << "Time"
                            << std::setw(_column_widths[1]) << "Tag";
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
          OutputTime::_output << std::setw(_column_widths[2])
                              << "Position_" + std::to_string(dd);
        OutputTime::_output << std::setw(_column_widths[3]) << "Mass";
        for (std::size_t ii = 0; ii < _masks.size(); ++ii)
          OutputTime::_output << std::setw(_column_widths[4])
                              << "In_region_" + std::to_string(ii);
        OutputTime::_output << std::setw(_column_widths[1]) << "Tag"
                            << std::setw(_column_widths[2]) << "..."
                            << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          auto cell = OutputTime::_locator(state);
          OutputTime::_output << std::setw(_column_widths[1]) << state.tag;
          if (!state.info.absorbed
              && part.state_old().time <= time
              && !outside(cell))
          {
            std::vector<int> in_region(_masks.size(), 0);
            for (std::size_t ii = 0; ii < _masks.size(); ++ii)
              if (_masks[ii].get()[cell] > _thresholds[ii])
                in_region[ii] = 1;
            if (std::any_of(in_region.begin(), in_region.end(),
                            [](int ii){ return ii > 0; }))
            {
              useful::print(OutputTime::_output, _getter_position(state), _column_widths[2]);
              OutputTime::_output << std::setw(_column_widths[3]) << state.mass;
              useful::print(OutputTime::_output, in_region, _column_widths[4]);
            }
          }
        }
        OutputTime::_output << "\n";
      }
      
    private:
      using Boundary = std::decay_t<decltype(std::declval<Geometry>().boundary_periodic)>;
      ctrw::Get_position_periodic<Boundary const&> _getter_position;
      std::vector<std::reference_wrapper<const Mask>> _masks;
      std::vector<double> _thresholds;
      std::array<int, 5> _column_widths;
    };
    
    /** \class Output_Cases::Output_position_mean_periodic PTOF/Output.h "PTOF/Output.h"
     *  \brief Output time and mean position (weighted by mass), with positions accounting for periodicity. */
    struct Output_position_mean_periodic final
    : OutputTime
    {
      Output_position_mean_periodic
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "position_mean_periodic", identifier, parameters,
          precision }
      , getter_position{ geometry.boundary_periodic }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(12, int(2 + std::string{ "Position_mean_" }.length())) }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time";
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
          OutputTime::_output << std::setw(_column_widths[1])
                              << "Position_mean_" + std::to_string(dd);
        OutputTime::_output << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        useful::print(OutputTime::_output,
                      position_mean(OutputTime::_subject, time,
                                    getter_position),
                      _column_widths[1]);
        OutputTime::_output << "\n";
      }
      
    private:
      using Boundary = std::decay_t<decltype(std::declval<Geometry>().boundary_periodic)>;
      ctrw::Get_position_periodic<Boundary const&> getter_position;
      std::array<int, 2> _column_widths;
    };
        
    /** \class Output_Cases::Output_position_second_moment_periodic PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and second moment of position (weighted by mass), with position accounting for periodicity. */
    struct Output_position_second_moment_periodic final
    : OutputTime
    {
      Output_position_second_moment_periodic
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "position_second_moment_periodic", identifier, parameters,
          precision }
      , getter_position{ geometry.boundary_periodic }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(9 + precision, int(2 + std::string{ "Position_second_moment_" }.length())) }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time";
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
          OutputTime::_output << std::setw(_column_widths[1])
                              << "Position_second_moment_" + std::to_string(dd);
        OutputTime::_output << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        useful::print(OutputTime::_output,
                      position_second_moment(OutputTime::_subject, time,
                                             getter_position),
                      _column_widths[1]);
        OutputTime::_output << "\n";
      }
      
      private:
        using Boundary = std::decay_t<decltype(std::declval<Geometry>().boundary_periodic)>;
        ctrw::Get_position_periodic<Boundary const&> getter_position;
        std::array<int, 2> _column_widths;
    };
        
    /** \class Output_Cases::Output_position_variance_periodic PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and position variance (weighted by mass), with position accounting for periodicity. */
    struct Output_position_variance_periodic final
    : OutputTime
    {
      Output_position_variance_periodic
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : OutputTime{ subject,
          geometry, directories,
          "position_variance_periodic", identifier, parameters,
          precision }
      , getter_position{ geometry.boundary_periodic }
      , _column_widths{
        std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
        std::max(9 + precision, int(2 + std::string{ "Position_variance_" }.length())) }
      {
        OutputTime::_output << std::setw(_column_widths[0]) << "Time";
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
          OutputTime::_output << std::setw(_column_widths[1])
                              << "Position_variance_" + std::to_string(dd);
        OutputTime::_output << "\n";
      }
      
      void operator()(double time) override
      {
        OutputTime::_output << std::setw(_column_widths[0]) << time;
        useful::print(OutputTime::_output,
                      position_variance(OutputTime::_subject, time,
                                        getter_position),
                      _column_widths[1]);
        OutputTime::_output << "\n";
      }
      
      private:
        using Boundary = std::decay_t<decltype(std::declval<Geometry>().boundary_periodic)>;
        ctrw::Get_position_periodic<Boundary const&> getter_position;
        std::array<int, 2> _column_widths;
    };
        
    /** \struct Output_Cases::Output_absorption_time PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output absorption times, tags, masses of absorbed particles. */
    struct Output_absorption_time final : Output
    {
      Output_absorption_time
      (Subject const& subject,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8)
      : Output{ subject,
        geometry, directories,
        "absorption_time", identifier, parameters,
        precision }
      , _column_widths{
          std::max(9 + precision, int(1 + std::string{ "Time" }.length())),
          std::max(12, int(1 + std::string{ "Tag" }.length())),
          std::max(9 + precision, int(1 + std::string{ "Mass" }.length())) }
      {
        Output::_output << std::setw(_column_widths[0]) << "Time"
                        << std::setw(_column_widths[1]) << "Tag"
                        << std::setw(_column_widths[2]) << "Mass"
                        << "\n";
      }
      
      void operator()() override
      {
        for (auto const& part : Output::_subject.particles())
        {
          auto const& state = part.state_new();
          if (state.info.absorbed)
            Output::_output << std::setw(_column_widths[0]) << state.time
                            << std::setw(_column_widths[1]) << state.tag
                            << std::setw(_column_widths[2]) << state.mass
                            << "\n";
        }
      }
      
      private:
        std::array<int, 3> _column_widths;
    };
    
    Subject const& _subject;                                /**< Subject to measure. */
    Geometry const& _geometry;                              /**< Domain geometry info and utilities.*/
    std::unique_ptr<Criterion<Subject>> _end_criterion;     /**< To check if end criterion is met. */
    std::unique_ptr<NextMeasureTime> _next_measure;         /**< Handle next measure time. */
    std::vector<std::unique_ptr<OutputTime>> _output_time;  /**< Handle each output type given time. */
    std::vector<std::unique_ptr<Output>> _output;           /**< Handle each output type given nothing. */
  };
}
        
#endif /* PTOF_OUTPUT_H */
