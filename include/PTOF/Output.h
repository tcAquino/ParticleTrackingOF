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
#include <exception>
#include <fstream>
#include <memory>
#include <iomanip>
#include <string>
#include <type_traits>
#include <vector>
#include <zero.H>
#include "CTRW/StateGetter.h"
#include "General/Useful.h"
#include "PTOF/Directories.h"
#include "PTOF/Field.h"
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
      time,               /**< Specified time. */
      time_max,           /**< Maximum output time. */
      mass_below,         /**< Total mass below value. */
      mass_above,         /**< Total mass above value. */
      all_absorbed,       /**< All particles absorbed. */
      one_absorbed,       /**< One particle absorbed. */
      fraction_absorbed   /**< A particle fraction absorbed. */
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
      { "fraction_absorbed", Type::fraction_absorbed }
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
      { Type::fraction_absorbed, "fraction_absorbed" }
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
      position,                         /**< Time and particle tags, positions, and masses.                                 */
      position_in_regions,              /**< Time and particle tags, positions, and masses in regions speciefied by masks. */
      position_mean,                    /**< Time and mean position.                                                        */ 
      position_second_moment,           /**< Time and position second moment.                                               */ 
      position_variance,                /**< Time and position variance.                                                    */ 
      mass,                             /**< Time and total mass.                                                           */
      mass_in_regions,                  /**< Time and total mass in regions speciefied by masks.                                            */
      velocity,                         /**< Time and particle tags and local velocities.                                   */ 
      velocity_gradient,                /**< Time and particle tags and local velocity gradient components.                 */ 
      pressure,                         /**< Time and particle tags and local pressures.                                    */ 
      position_periodic,                /**< Time and particle tags, true positions accounting for periodicity, and masses. */ 
      position_mean_periodic,           /**< Time and true mean position accounting for periodicity.                        */ 
      position_second_moment_periodic,  /**< Time and true position second moment accounting for periodicity.               */ 
      position_variance_periodic,       /**< Time and true position variance accounting for periodicity.                    */ 
      absorption_time                   /**< Particle absorption times, tags, and masses at end of dynamics.                */ 
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
      { "pressure", Type::pressure },
      { "position_periodic", Type::position_periodic },
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
      { Type::pressure, "pressure" },
      { Type::position_periodic, "position_periodic" },
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
   *  \brief  Check if mass is greater than or equal to value. */
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
   * \brief  Check if all particles have been absorbed. */
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
  
  /** \class Criterion_fraction_absorbed PTOF/Output.h "PTOF/Output.h"
   *  \brief Check if at a least a given fraction of particles have been absorbed. */
  template <typename Subject>
  struct Criterion_fraction_absorbed final : Criterion<Subject>
  {
    Criterion_fraction_absorbed(double end_value)
    : end_value{ end_value }
    {}
    
    bool operator()
    (Subject const& subject, double time) const override
    {
      return nr_absorbed(subject, time) >= end_value*subject.size();
    }
      
    double end_value;
  };
  
  /** \class Output_Cases PTOF/Output.h "PTOF/Output.h"
   * \brief Output object to handle implemented output options. */
  template
  <typename Subject, typename VelocityField, typename Geometry>
  class Output_Cases
  {
  public:
    /** \struct Output_Cases::Parameters PTOF/Output.h "PTOF/Output.h"
     *  \brief Output parameters. */
    struct Parameters
    {
      double velocity_rescaling;
      std::size_t run_nr;
      std::string time_units;
      double time_unit_factor;
      std::string end_criterion;
      double end_value;
      std::string measure_spacing;
      double time_min;
      double time_max;
      double time_increment;
      std::vector<std::string> measure_names;
      
      /** Constructor.
       \param directories Current case directory information.
       \param name Name of parameter set.
       \param params_transport Transport parameters.
       \param params_reaction Reaction parameters.
       \param params_solvers Solver parameters.
       \param velocity_rescaling Factor by which velocity field was rescaled, used to rescale, e.g., velocity gradients. */
      template
      <typename TransportParameters,
      typename ReactionParameters,
      typename SolverParameters>
      Parameters
      (Directories const& directories,
       std::string const& name,
       TransportParameters const& params_transport,
       ReactionParameters const& params_reaction,
       SolverParameters const& params_solvers,
       double velocity_rescaling)
      : velocity_rescaling{ velocity_rescaling }
      {
        auto input = useful::open_read(directories.dir_parameters
                                       + "/parameters_output_"
                                       + name + ".dat");
        useful::read(input, run_nr);
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
          "- Run number (nonnegative integer to index output)\n"
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
          "\tfraction_absorbed: Fraction of particles absorbed\n"
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
          case EndCriterion::Type::fraction_absorbed:
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
        while (1)
        {
          std::string measurement;
          input >> measurement;
          if (!input.fail())
            measure_names.push_back(measurement);
          else
            break;
        }
      }
    };
    
    /** Constructor.
     \param subject CTRW object to measure.
     \param velocity_field Velocity field as a function of state.
     \param geometry Domain geometry info and utilities.
     \param directories Current case directory information.
     \param parameters Output parameters.
     \param identifier String to include in names of output files.
     \param masks Vector of pointers to Masks. Masks are scalar fields assigning values to mesh cell indices through operator[].
     \param thresholds Vector of thresholds for each mask, such that cells where a mask is above the threshold are considered.
     \param precision Number of digits after decimal point in output, in scientific notation.
     \param delimiter Delimiter between output values on same line.
    */
    template <typename Mask = useful::Empty>
    Output_Cases
    (Subject const& subject,
     VelocityField const& velocity_field,
     Geometry const& geometry,
     Directories const& directories,
     Parameters parameters,
     std::string const& identifier,
     std::vector<Mask const*> masks = {},
     std::vector<double> thresholds = {},
     int precision = 8, std::string delimiter = "\t")
    : parameters{ parameters }
    , _subject{ subject }
    , _velocity_field{ velocity_field }
    , _geometry{ geometry }
    {
      set_measure_types(directories, identifier,
                        masks, thresholds,
                        precision, delimiter);
      set_end_criterion();
      set_next_measure_time();
    }
    
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
      for (auto const& measure : _output_time)
        measure->operator()(time);
      _next_measure->advance();
    }
    
    /** \brief Output current information. */
    void operator()()
    {
      for (auto const& measure : _output)
        measure->operator()();
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
        "- Run number: " << parameters.run_nr << "\n"
        "- End criterion: " << parameters.end_criterion << "\n"
        "- Time units: " << parameters.time_units << "\n"
        "- Measure spacing: " << parameters.measure_spacing << "\n"
        "- Minimum measurement time: " << parameters.time_min << "\n"
        "- Maximum measurement time: " << std::to_string(parameters.time_max) << "\n"
        "- Measurement types:";
      useful::print(output, parameters.measure_names, true, "\n\t");
      output << "\n";
      output <<
        "--------------------------------------------------\n";
    }
        
  private:
    /** \brief Set up output streams for requested output types.
     \param directories Current case directory information.
     \param identifier String to include in names of output files.
     \param masks Vector of pointers to Masks. Masks are scalar fields assigning values to mesh cell indices through operator[].
     \param thresholds Vector of thresholds for each mask, such that cells where a mask is above the threshold are considered.
     \param precision Number of digits after decimal point in output, in scientific notation.
     \param delimiter Delimiter between output values on same line. */
    template <typename Mask = useful::Empty>
    void set_measure_types
    (Directories const& directories,
     std::string const& identifier,
     std::vector<Mask const*> masks = {},
     std::vector<double> thresholds = {},
     int precision = 8,
     std::string delimiter = "\t")
    {
      for (auto const& name : parameters.measure_names)
      {
        switch (Measure::type(name))
        {
          case Measure::Type::position:
          {
            _output_time.emplace_back
            (std::make_unique<Output_position>
             (_subject, _velocity_field, _geometry,
               directories, name, identifier, parameters,
               precision, delimiter));
            break;
          }
          case Measure::Type::position_in_regions:
          {
            _output_time.emplace_back
            (std::make_unique<Output_position_in_regions<Mask>>
             (_subject, _velocity_field, _geometry,
               directories, name, identifier, parameters,
               masks, thresholds,
               precision, delimiter));
            break;
          }
          case Measure::Type::position_mean:
          {
            _output_time.emplace_back
            (std::make_unique<Output_position_mean>
             (_subject, _velocity_field, _geometry,
              directories, name, identifier, parameters,
              precision, delimiter));
            break;
          }
          case Measure::Type::position_second_moment:
          {
            _output_time.emplace_back
            (std::make_unique<Output_position_second_moment>
             (_subject, _velocity_field, _geometry,
              directories, name, identifier, parameters,
              precision, delimiter));
            break;
          }
          case Measure::Type::position_variance:
          {
            _output_time.emplace_back
            (std::make_unique<Output_position_variance>
             (_subject, _velocity_field, _geometry,
              directories, name, identifier, parameters,
              precision, delimiter));
            break;
          }
          case Measure::Type::mass:
          {
            _output_time.emplace_back
            (std::make_unique<Output_mass>
             (_subject, _velocity_field, _geometry,
              directories, name, identifier, parameters,
              precision, delimiter));
            break;
          }
          case Measure::Type::mass_in_regions:
          {
            _output_time.emplace_back
            (std::make_unique<Output_mass_in_regions<Mask>>
             (_subject, _velocity_field, _geometry,
              directories, name, identifier, parameters,
              masks, thresholds,
              precision, delimiter));
            break;
          }
          case Measure::Type::velocity:
          {
            _output_time.emplace_back
            (std::make_unique<Output_velocity>
             (_subject, _velocity_field, _geometry,
              directories, name, identifier, parameters,
              precision, delimiter));
            break;
          }
          case Measure::Type::velocity_gradient:
          {
            _output_time.emplace_back
            (std::make_unique<Output_velocity_gradient>
             (_subject, _velocity_field, _geometry,
              directories, name, identifier, parameters,
              precision, delimiter));
            break;
          }
          case Measure::Type::pressure:
          {
            _output_time.emplace_back
            (std::make_unique<Output_pressure>
             (_subject, _velocity_field, _geometry,
              directories, name, identifier, parameters,
              precision, delimiter));
            break;
          }
          case Measure::Type::position_periodic:
          {
            if constexpr (Has_periodicity<
                          typename Subject::Particle::State>::value
                          == true)
              _output_time.emplace_back
              (std::make_unique<Output_position_periodic<
               decltype(_geometry.boundary_periodic)>>
               (_subject, _velocity_field, _geometry,
                directories, name, identifier, parameters,
                precision, delimiter));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + name + ": "
                "Particle state must define periodicity" };
            break;
          }
          case Measure::Type::position_mean_periodic:
          {
            if constexpr (Has_periodicity<
                          typename Subject::Particle::State>::value
                          == true)
              _output_time.emplace_back
              (std::make_unique<Output_position_mean_periodic<
               decltype(_geometry.boundary_periodic)>>
               (_subject, _velocity_field, _geometry,
                directories, name, identifier, parameters,
                precision, delimiter));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + name + ": "
                "Particle state must define periodicity" };
            break;
          }
          case Measure::Type::position_second_moment_periodic:
          {
            if constexpr (Has_periodicity<
                          typename Subject::Particle::State>::value
                          == true)
              _output_time.emplace_back
              (std::make_unique<Output_position_second_moment_periodic<
               decltype(_geometry.boundary_periodic)>>
               (_subject, _velocity_field, _geometry,
                directories, name, identifier, parameters,
                precision, delimiter));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + name + ": "
                "Particle state must define periodicity" };
            break;
          }
          case Measure::Type::position_variance_periodic:
          {
            if constexpr (Has_periodicity<
                          typename Subject::Particle::State>::value
                          == true)
              _output_time.emplace_back
              (std::make_unique<Output_position_variance_periodic<
               decltype(_geometry.boundary_periodic)>>
               (_subject, _velocity_field, _geometry,
                directories, name, identifier, parameters,
                precision, delimiter));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + name + ": "
                "Particle state must define periodicity" };
            break;
          }
          case Measure::Type::absorption_time:
          {
            _output.emplace_back
            (std::make_unique<Output_absorption_time>
             (_subject, _velocity_field, _geometry,
              directories, name, identifier, parameters,
              precision, delimiter));
            break;
          }
          default:
            throw std::runtime_error{
              std::string("Measurement type ") + name + " "
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
        case EndCriterion::Type::fraction_absorbed:
        {
          _end_criterion = std::make_unique<
            Criterion_fraction_absorbed<Subject>>(parameters.end_value);
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
      
    protected:
      /** \brief Set up output streams for requested output types.
      \param subject CTRW object to measure.
      \param velocity_field Velocity field as a function of state.
      \param geometry Domain geometry info and utilities.
      \param directories Current case directory information.
      \param output_name Name of output type.
      \param identifier String to include in names of output files.
      \param parameters Output parameters.
      \param precision Number of digits after decimal point in output, in scientific notation.
      \param delimiter Delimiter between output values on same line. */
      OutputTime
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : _subject{ subject }
      , _velocity_field{ velocity_field }
      , _geometry{ geometry }
      , _locator{ geometry.locator }
      , _output{ open_write(directories, output_name, identifier, parameters.run_nr) }
      , _delimiter{ delimiter }
      {
        _output << std::setprecision(precision)
                << std::scientific;
      }
      
      Subject const& _subject;                      /**< CTRW object to measure. */
      VelocityField const& _velocity_field;         /**< Velocity field as a function of state. */
      Geometry const& _geometry;                    /**< Domain geometry info and utilities. */
      typename Geometry::Locator const& _locator;   /**< Object to locate positions in mesh. */
      std::ofstream _output;                        /**< Output stream. */
      std::string _delimiter;                       /**< Delimiter between output values on same line. */
        
      /** \brief Open output file for a given output type, tagged with run number identifier.
       \param directories Current case directory information.
       \param output_name Name of output type.
       \param identifier String to include in names of output files.
       \param run_nr Unsigned integer to index output file.
       \return Output stream. */
      std::ofstream open_write
      (Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       std::size_t run_nr)
      {
        return
          useful::open_write(directories.dir_output
          + "/Data_" + output_name
          + (identifier.empty() ? "" : "_" + identifier)
          + "_RUN_" + std::to_string(run_nr)
          + ".dat") ;
      }
        
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
      \param velocity_field Velocity field as a function of state.
      \param geometry Domain geometry info and utilities.
      \param directories Current case directory information.
      \param output_name Name of output type.
      \param identifier String to include in names of output files.
      \param parameters Output parameters.
      \param precision Number of digits after decimal point in output, in scientific notation.
      \param delimiter Delimiter between output values on same line. */
      Output
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : _subject{ subject }
      , _velocity_field{ velocity_field }
      , _geometry{ geometry }
      , _locator{ geometry.locator }
      , _output{ open_write(directories, output_name, identifier, parameters.run_nr) }
      , _delimiter{ delimiter }
      {
        _output << std::setprecision(precision)
                << std::scientific;
      }
      
      Subject const& _subject;                      /**< CTRW object to measure. */
      VelocityField const& _velocity_field;         /**< Velocity field as a function of state. */
      Geometry const& _geometry;                    /**< Domain geometry info and utilities. */
      typename Geometry::Locator const& _locator;   /**< Object to locate positions in mesh. */
      std::ofstream _output;                        /**< Output stream. */
      std::string _delimiter;                       /**< Delimiter between output values on same line. */
      
      /** \brief Open output file for a given output type, tagged with run number identifier.
      \param directories Current case directory information.
      \param output_name Name of output type.
      \param identifier String to include in names of output files.
      \param run_nr Unsigned integer to index output file.
      \return Output stream. */
      std::ofstream open_write
      (Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       std::size_t run_nr)
      {
        return
          useful::open_write(directories.dir_output
          + "/Data_" + output_name
          + (identifier.empty() ? "" : "_" + identifier)
          + "_RUN_" + std::to_string(run_nr)
          + ".dat");
      }
      
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
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::_output << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed
              && part.state_old().time <= time)
          {
            OutputTime::_output << OutputTime::_delimiter << state.tag;
            useful::print(OutputTime::_output, state.position,
                          true, OutputTime::_delimiter);
            OutputTime::_output << OutputTime::_delimiter << state.mass;
          }
        }
        OutputTime::_output << "\n";
      }
    };
        
    /** \class Output_Cases::Output_position_in_regions PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time, tags, positions, and masses within regions specified by masks. */
    template <typename Mask>
    struct Output_position_in_regions final : OutputTime
    {
      Output_position_in_regions
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       std::vector<Mask const*> masks,
       std::vector<double> thresholds = {},
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      , masks{ masks }
      , thresholds{ thresholds }
      {}
      
      void operator()(double time) override
      {
        OutputTime::_output << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          if constexpr (!std::is_same_v<Mask, useful::Empty>)
          {
            auto const& state = part.state_new();
            auto cell = OutputTime::_locator(state);
            for (std::size_t ii = 0; ii < masks.size(); ++ii)
            {
              if (!state.info.absorbed
                  && part.state_old().time <= time
                  && cell >= 0
                  && (*masks[ii])[cell] > thresholds[ii])
              {
                OutputTime::_output << OutputTime::_delimiter << state.tag;
                useful::print(OutputTime::_output, state.position,
                              true, OutputTime::_delimiter);
                OutputTime::_output << OutputTime::_delimiter << state.mass;
              }
            }
          }
        }
        OutputTime::_output << "\n";
      }
      
      std::vector<Mask const*> masks;
      std::vector<double> thresholds;
    };
    
    /** \struct Output_Cases::Output_position_mean PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and mean position (weighted by mass). */
    struct Output_position_mean final : OutputTime
    {
      Output_position_mean
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::_output << time;
        useful::print(OutputTime::_output,
                      position_mean(OutputTime::_subject, time),
                      true, OutputTime::_delimiter);
        OutputTime::_output << "\n";
      }
    };
        
    /** \struct Output_Cases::Output_position_second_moment PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and second moment of position (weighted by mass). */
    struct Output_position_second_moment final : OutputTime
    {
      Output_position_second_moment
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::_output << time;
        useful::print(OutputTime::_output,
                      position_second_moment(OutputTime::_subject, time),
                      true, OutputTime::_delimiter);
        OutputTime::_output << "\n";
      }
    };
        
    /** \struct Output_Cases::Output_position_variance PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and position variance (weighted by mass). */
    struct Output_position_variance final : OutputTime
    {
      Output_position_variance
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::_output << time;
        useful::print(OutputTime::_output,
                      position_variance(OutputTime::_subject, time),
                      true, OutputTime::_delimiter);
        OutputTime::_output << "\n";
      }
    };
        
    /** \struct Output_Cases::Output_mass PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and total mass. */
    struct Output_mass final : OutputTime
    {
      Output_mass
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::_output << time << OutputTime::_delimiter
                            << mass(OutputTime::_subject, time) << "\n";
      }
    };
        
    /** \class Output_Cases::Output_mass_in_regions PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and total mass in regions specified by maks. */
    template <typename Mask>
    struct Output_mass_in_regions final : OutputTime
    {
      Output_mass_in_regions
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       std::vector<Mask const*> masks,
       std::vector<double> thresholds = {},
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      , masks{ masks }
      , thresholds{ thresholds }
      {
        while(thresholds.size() < masks.size())
          thresholds.push_back(0.);
      }
      
      void operator()(double time) override
      {
        if constexpr (!std::is_same_v<Mask, useful::Empty>)
        {
          auto masses = mass(OutputTime::_subject, time,
                             OutputTime::_locator,
                             masks, thresholds);
          OutputTime::_output << time;
          useful::print(OutputTime::_output, masses, true, OutputTime::_delimiter);
          OutputTime::_output << "\n";
        }
      }
      
      std::vector<Mask const*> masks;
      std::vector<double> thresholds;
    };
        
     /** \struct Output_Cases::Output_velocity PTOF/Output.h "PTOF/Output.h"
      *  \brief  Output time, tags, velocities. */
    struct Output_velocity final : OutputTime
    {
      Output_velocity
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      {}

      void operator()(double time) override
      {
        OutputTime::_output << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed
              && part.state_old().time <= time)
          {
            OutputTime::_output << OutputTime::_delimiter << state.tag;
            useful::print(OutputTime::_output,
                          OutputTime::_velocity_field(state),
                          true, OutputTime::_delimiter);
          }
        }
        OutputTime::_output << "\n";
      }
    };
        
     /** \struct Output_Cases::Output_velocity_gradient PTOF/Output.h "PTOF/Output.h"
      *  \brief  Output time and velocity gradient components.
      *  \details
      * Component output order is xx, xy, xz, yx, yy, yz, zx, zy, zz.
      * Spurious components are omitted if geometry dimension is < 3. */
    struct Output_velocity_gradient final : OutputTime
    {
      Output_velocity_gradient
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      , velocity_gradient{
        Foam::IOobject{
          "gradU", geometry.mesh().time().timeName(),
          geometry.mesh(), Foam::IOobject::MUST_READ,
          Foam::IOobject::NO_WRITE },
          geometry.mesh() }
      {
        if (parameters.velocity_rescaling != 1.)
          velocity_gradient *= parameters.velocity_rescaling;
      }

      void operator()(double time) override
      {
        OutputTime::_output << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed
              && part.state_old().time <= time)
          {
            OutputTime::_output << OutputTime::_delimiter << state.tag;
            auto cell = OutputTime::_locator(state);
            if (cell < 0)
              useful::print(OutputTime::_output,
                            std::vector<double>(Geometry::dim, 0.),
                            true, OutputTime::_delimiter);
            else
            {
              auto const& grad = velocity_gradient[cell];
              if constexpr (Geometry::dim == 1)
              {
                OutputTime::_output << OutputTime::_delimiter << grad.xx();
              }
              if constexpr (Geometry::dim == 2)
              {
                OutputTime::_output << OutputTime::_delimiter << grad.xx();
                OutputTime::_output << OutputTime::_delimiter << grad.xy();
                OutputTime::_output << OutputTime::_delimiter << grad.yx();
                OutputTime::_output << OutputTime::_delimiter << grad.yy();
              }
              if constexpr (Geometry::dim == 3)
              {
                OutputTime::_output << OutputTime::_delimiter << grad.xx();
                OutputTime::_output << OutputTime::_delimiter << grad.xy();
                OutputTime::_output << OutputTime::_delimiter << grad.xz();
                OutputTime::_output << OutputTime::_delimiter << grad.yx();
                OutputTime::_output << OutputTime::_delimiter << grad.yy();
                OutputTime::_output << OutputTime::_delimiter << grad.yz();
                OutputTime::_output << OutputTime::_delimiter << grad.zx();
                OutputTime::_output << OutputTime::_delimiter << grad.zy();
                OutputTime::_output << OutputTime::_delimiter << grad.zz();
              }
            }
          }
        }
        OutputTime::_output << "\n";
      }

      Foam::volTensorField velocity_gradient;
    };
        
    /**
     \struct Output_Cases::Output_pressure PTOF/Output.h "PTOF/Output.h"
     \brief Output time, tags, and pressures.
    */
    struct Output_pressure final : OutputTime
    {
      Output_pressure
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      , pressure{
        { Foam::IOobject{
          "U", geometry.mesh().time().timeName(),
          geometry.mesh(), Foam::IOobject::MUST_READ,
          Foam::IOobject::NO_WRITE },
          geometry.mesh() },
        geometry.locator,
        CheckOptions::Check{} }
      {
        if (parameters.velocity_rescaling != 1.)
          pressure.rescale(parameters.velocity_rescaling);
      }

      void operator()(double time) override
      {
        OutputTime::_output << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed
              && part.state_old().time <= time)
          {
            OutputTime::_output << OutputTime::_delimiter << state.tag;
            useful::print(OutputTime::_output, pressure(state),
                          true, OutputTime::_delimiter);
          }
        }
        OutputTime::_output << "\n";
      }
      
      ptof::ScalarField_LinearInterpolation_OF
      <Foam::volScalarField, typename Geometry::Locator const&, 1, 0> pressure;
    };
        
    /** \class Output_Cases::Output_position_periodic PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time, tags, positions, and masses, with positions accounting for periodicity. */
    template <typename Boundary>
    struct Output_position_periodic final : OutputTime
    {
      Output_position_periodic
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      , getter_position{ geometry.boundary_periodic }
      {}
      
      void operator()(double time) override
      {
        OutputTime::_output << time;
        for (auto const& part : OutputTime::_subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed
              && part.state_old().time <= time)
          {
            OutputTime::_output << OutputTime::_delimiter << state.tag;
            useful::print(OutputTime::_output, getter_position(state),
                          true, OutputTime::_delimiter);
            OutputTime::_output << OutputTime::_delimiter << state.mass;
          }
        }
        OutputTime::_output << "\n";
      }
      
      ctrw::Get_position_periodic<Boundary const&> getter_position;
    };
    
    /** \class Output_Cases::Output_position_mean_periodic PTOF/Output.h "PTOF/Output.h"
     *  \brief Output time and mean position (weighted by mass), with positions accounting for periodicity. */
    template <typename Boundary>
    struct Output_position_mean_periodic final
    : OutputTime
    {
      Output_position_mean_periodic
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      , getter_position{ geometry.boundary_periodic }
      {}
      
      void operator()(double time) override
      {
        OutputTime::_output << time;
        useful::print(OutputTime::_output,
                      position_mean(OutputTime::_subject, time,
                                    getter_position),
                      true, OutputTime::_delimiter);
        OutputTime::_output << "\n";
      }
      
      ctrw::Get_position_periodic<Boundary const&> getter_position;
    };
        
    /** \class Output_Cases::Output_position_second_moment_periodic PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and second moment of position (weighted by mass), with position accounting for periodicity. */
    template <typename Boundary>
    struct Output_position_second_moment_periodic final
    : OutputTime
    {
      Output_position_second_moment_periodic
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      , getter_position{ geometry.boundary_periodic }
      {}
      
      void operator()(double time) override
      {
        OutputTime::_output << time;
        useful::print(OutputTime::_output,
                      position_second_moment(OutputTime::_subject, time,
                                             getter_position),
                      true, OutputTime::_delimiter);
        OutputTime::_output << "\n";
      }
      
      ctrw::Get_position_periodic<Boundary const&> getter_position;
    };
        
    /** \class Output_Cases::Output_position_variance_periodic PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and position variance (weighted by mass), with position accounting for periodicity. */
    template <typename Boundary>
    struct Output_position_variance_periodic final
    : OutputTime
    {
      Output_position_variance_periodic
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
          geometry, directories,
          output_name, identifier, parameters,
          precision, delimiter }
      , getter_position{ geometry.boundary_periodic }
      {}
      
      void operator()(double time) override
      {
        OutputTime::_output << time;
        useful::print(OutputTime::_output,
                      position_variance(OutputTime::_subject, time,
                                        getter_position),
                      true, OutputTime::_delimiter);
        OutputTime::_output << "\n";
      }
      
      ctrw::Get_position_periodic<Boundary const&> getter_position;
    };
        
    /** \struct Output_Cases::Output_absorption_time PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output absorption times, tags, masses of absorbed particles. */
    struct Output_absorption_time final : Output
    {
      Output_absorption_time
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& parameters,
       int precision = 8,
       std::string delimiter = "\t")
      : Output{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, parameters,
        precision, delimiter }
      {}
      
      void operator()() override
      {
        for (auto const& part : Output::_subject.particles())
        {
          auto const& state = part.state_new();
          if (state.info.absorbed)
          {
            Output::_output << state.time << Output::_delimiter;
            Output::_output << state.tag << Output::_delimiter;
            Output::_output << state.mass << "\n";
          }
        }
      }
    };
    
    Subject const& _subject;                                /**< Subject to measure.                   */
    VelocityField const& _velocity_field;                   /**< Compute velocity from particle state. */
    Geometry const& _geometry;                              /**< Mesh and other domain information.    */
    std::unique_ptr<Criterion<Subject>> _end_criterion;     /**< To check if end criterion is met.     */
    std::unique_ptr<NextMeasureTime> _next_measure;         /**< Handle next measure time.             */
    std::vector<std::unique_ptr<OutputTime>> _output_time;  /**< Handle each output type given time.   */
    std::vector<std::unique_ptr<Output>> _output;           /**< Handle each output type given nothing.*/
  };
}
        
#endif /* PTOF_OUTPUT_H */
