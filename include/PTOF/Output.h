/**
* \file PTOF/Output.h
* \author Tomás Aquino
* \date 07/03/2022
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
      time,               /**< Specified time               */ 
      time_max,           /**< Maximum output time          */ 
      mass_below,         /**< Total mass below value       */ 
      mass_above,         /**< Total mass above value       */ 
      all_absorbed,       /**< All particles absorbed       */ 
      one_absorbed,       /**< One particle absorbed        */ 
      fraction_absorbed   /**< A particle fraction absorbed */ 
    };
    
    /** Type from name. */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** Name from type. */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** Check if name exists. */
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    /** Map names to types. */
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
    
    /** Map types to names. */
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
    
    /** Type from name. */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** Name from type. */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** Check if name exists. */
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    /** Map names to types. */
    inline static const
    std::unordered_map<std::string, Type> name_to_type
    {
      { "linear", Type::linear },
      { "log", Type::log },
      { "step", Type::step }
    };
    
    /** Map types to names. */
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
      diffusion, /**<  */ 
      advection, /**<  */ 
      reaction,  /**<  */ 
      arbitrary  /**<  */ 
    };
    
    /** Type from name. */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** Name from type. */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** Check if name exists. */
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    /** Map names to types. */
    inline static const
    std::unordered_map<std::string, Type>
    name_to_type
    {
      { "diffusion", Type::diffusion },
      { "advection", Type::advection },
      { "reaction", Type::reaction },
      { "arbitrary", Type::arbitrary }
    };
    
    /** Map types to names. */
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
      position_mean,                    /**< Time and mean position.                                                        */ 
      position_second_moment,           /**< Time and position second moment.                                               */ 
      position_variance,                /**< Time and position variance.                                                    */ 
      mass,                             /**< Time and total mass.                                                           */ 
      velocity,                         /**< Time and particle tags and local velocities.                                   */ 
      velocity_gradient,                /**< Time and particle tags and local velocity gradient components.                 */ 
      pressure,                         /**< Time and particle tags and local pressures.                                    */ 
      position_periodic,                /**< Time and particle tags, true positions accounting for periodicity, and masses. */ 
      position_mean_periodic,           /**< Time and true mean position accounting for periodicity.                        */ 
      position_second_moment_periodic,  /**< Time and true position second moment accounting for periodicity.               */ 
      position_variance_periodic,       /**< Time and true position variance accounting for periodicity.                    */ 
      absorption_time                   /**< Particle absorption times, tags, and masses at end of dynamics.                */ 
    };
    
    /** Type from name. */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** Name from type. */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** Check if name exists. */
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    /** Map names to types. */
    inline static const
    std::unordered_map<std::string, Type>name_to_type
    {
      { "position", Type::position },
      { "position_mean", Type::position_mean },
      { "position_second_moment", Type::position_second_moment },
      { "position_variance", Type::position_variance },
      { "mass", Type::mass },
      { "velocity", Type::velocity },
      { "velocity_gradient", Type::velocity_gradient },
      { "pressure", Type::pressure },
      { "position_periodic", Type::position_periodic },
      { "position_mean_periodic", Type::position_mean_periodic },
      { "position_second_moment_periodic", Type::position_second_moment_periodic },
      { "position_variance_periodic", Type::position_variance_periodic },
      { "absorption_time", Type::absorption_time }
    };
    
    /** Map types to names. */
    inline static const
    std::unordered_map<Type, std::string>type_to_name
    {
      { Type::position, "position" },
      { Type::position_mean, "position_mean" },
      { Type::position_second_moment, "position_second_moment" },
      { Type::position_variance, "position_variance" },
      { Type::mass, "mass" },
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
  
  /** \struct Criterion PTOF/Output.h "PTOF/Output.h"
   * \brief  Polymorphic criteria. */
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
  
  /** \struct Criterion_time final PTOF/Output.h "PTOF/Output.h
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
  
  /** \struct Criterion_mass_below final PTOF/Output.h "PTOF/Output.h"
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
  
  /** \struct Criterion_mass_above final PTOF/Output.h "PTOF/Output.h"
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
  
  /** \struct Criterion_all_absorbed final PTOF/Output.h "PTOF/Output.h"
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
  
  /** \struct Criterion_one_absorbed final PTOF/Output.h "PTOF/Output.h"
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
  
  /** \struct Criterion_fraction_absorbed final PTOF/Output.h "PTOF/Output.h"
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
    /** \struct Parameters PTOF/Output.h "PTOF/Output.h"
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
      
      /** Get Output parameters from:
       * - directories info; 
       * - parameter set name;
       * - transport parameters;
       * - reaction parameters;
       * - solver parameters. */
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
                                       + "/output.dat");
        useful::read(input, run_nr);
        read_time_units(input, params_transport, params_reaction);
        read_end_criterion(input);
        read_measure_spacing(input,
                             params_transport, params_reaction);
        read_output_types(input);
        input.close();
      }
      
      /** Output generic information about object. */
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
          "- Maximum measurement time, if required by measure spacing\n"
          "- Number of measurements, if required by measure spacing\n"
          "- Measurement types (any number):\n"
          "\tposition: Time, particle positions, and particle masses\n"
          "\tposition_mean: Time and mean position\n"
          "\tposition_second_moment: Time and position second moment\n"
          "\tposition_variance: Time and position variance\n"
          "\tmass: Time and total mass\n"
          "\tabsorption_time: Particle absorption times and masses at end of dynamics\n"
          "--------------------------------------------------\n";
      }
      
    private:
      /** Read end criterion from input stream. */
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
      
      /** Read measure spacing time units from input stream. */
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
      
      /** Read measure spacing type from input stream. */
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
      
      /** Read measure spacing type from input stream. */
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
    
    /** Construct given directory information, output parameters,
     * and identifier for output filenames.
     * Optionaly set output precision and delimiter character. */
    Output_Cases
    (Subject const& subject,
     VelocityField const& velocity_field,
     Geometry const& geometry,
     Directories const& directories,
     Parameters params,
     std::string const& identifier,
     int precision = 8, std::string delimiter = "\t")
    : params{ params }
    , subject{ subject }
    , velocity_field{ velocity_field }
    , geometry{ geometry }
    {
      set_measure_types(geometry, directories, identifier,
                        precision, delimiter);
      set_end_criterion();
      set_next_measure_time();
    }
    
    /** Check if end simulation criterion is satisfied. */
    bool done(double time) const
    {
      return end_criterion->operator()(subject, time);
    }

    /** Set time of next measurement. */
    void set_next_measure_time()
    {
      switch (MeasureSpacing::type(params.measure_spacing))
      {
        case MeasureSpacing::Type::step:
        {
          next_measure = std::make_unique<
            NextMeasureTime_step>(params);
          break;
        }
        case MeasureSpacing::Type::linear:
        {
          next_measure = std::make_unique<
            NextMeasureTime_linear>(params);
          break;
        }
        case MeasureSpacing::Type::log:
        {
          next_measure = std::make_unique<
            NextMeasureTime_log>(params);
          break;
        }
        default:
          throw std::runtime_error{
            "Measure spacing "
            + params.measure_spacing
            + " not supported" };
      }
    }
    
    /** Get time of next measurement. */
    double next_measure_time() const
    {
      return next_measure->time();
    }
        
    /** Output requested measurements at given time advance to next measurement. */
    void operator()(double time)
    {
      for (auto const& measure : output_time)
        measure->operator()(time);
      next_measure->advance();
    }
    
    /** Output current information. */
    void operator()()
    {
      for (auto const& measure : output)
        measure->operator()();
    }
    
    const Parameters params;  /**< Output parameters .*/
    
    /** Output information about current object. */
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Output\n"
        "--------------------------------------------------\n"
        "- Run number: " << params.run_nr << "\n"
        "- End criterion: " << params.end_criterion << "\n"
        "- Time units: " << params.time_units << "\n"
        "- Measure spacing: " << params.measure_spacing << "\n"
        "- Minimum measurement time: " << params.time_min << "\n"
        "- Maximum measurement time: " << std::to_string(params.time_max) << "\n"
        "- Measurement types:";
      useful::print(output, params.measure_names, 1, "\n\t");
      output << "\n";
      output <<
        "--------------------------------------------------\n";
    }
        
  private:
    /** Open output file for a given output type. */
    std::ofstream open_write
    (Directories const& directories,
     std::string const& output_name,
     std::string const& identifier)
    {
      return
        useful::open_write(directories.dir_output
        + "/" + output_name
        + ".dat") ;
    }
  
    /** Set up output streams for requested output types. */
    void set_measure_types
    (Geometry const& geometry,
     Directories const& directories,
     std::string const& identifier,
     int precision = 8,
     std::string delimiter = "\t")
    {
      for (auto const& name : params.measure_names)
      {
        switch (Measure::type(name))
        {
          case Measure::Type::position:
          {
            output_time.emplace_back
            (std::make_unique<Output_position>
              (subject, velocity_field, geometry,
               directories, name, identifier, params));
            break;
          }
          case Measure::Type::position_mean:
          {
            output_time.emplace_back
            (std::make_unique<Output_position_mean>
             (subject, velocity_field, geometry,
              directories, name, identifier, params));
            break;
          }
          case Measure::Type::position_second_moment:
          {
            output_time.emplace_back
            (std::make_unique<Output_position_second_moment>
             (subject, velocity_field, geometry,
              directories, name, identifier, params));
            break;
          }
          case Measure::Type::position_variance:
          {
            output_time.emplace_back
            (std::make_unique<Output_position_variance>
             (subject, velocity_field, geometry,
              directories, name, identifier, params));
            break;
          }
          case Measure::Type::mass:
          {
            output_time.emplace_back
            (std::make_unique<Output_mass>
             (subject, velocity_field, geometry,
              directories, name, identifier, params));
            break;
          }
          case Measure::Type::velocity:
          {
            output_time.emplace_back
            (std::make_unique<Output_velocity>
             (subject, velocity_field, geometry,
              directories, name, identifier, params));
            break;
          }
          case Measure::Type::velocity_gradient:
          {
            output_time.emplace_back
            (std::make_unique<Output_velocity_gradient>
             (subject, velocity_field, geometry,
              directories, name, identifier, params));
            break;
          }
          case Measure::Type::pressure:
          {
            output_time.emplace_back
            (std::make_unique<Output_pressure>
             (subject, velocity_field, geometry,
              directories, name, identifier, params));
            break;
          }
          case Measure::Type::position_periodic:
          {
            if constexpr (Has_periodicity<
                          typename Subject::Particle::State>::value
                          == true)
              output_time.emplace_back
              (std::make_unique<Output_position_periodic<
               decltype(geometry.boundary_periodic)>>
               (subject, velocity_field, geometry,
                directories, name, identifier, params));
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
              output_time.emplace_back
              (std::make_unique<Output_position_mean_periodic<
               decltype(geometry.boundary_periodic)>>
               (subject, velocity_field, geometry,
                directories, name, identifier, params));
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
              output_time.emplace_back
              (std::make_unique<Output_position_second_moment_periodic<
               decltype(geometry.boundary_periodic)>>
               (subject, velocity_field, geometry,
                directories, name, identifier, params));
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
              output_time.emplace_back
              (std::make_unique<Output_position_variance_periodic<
               decltype(geometry.boundary_periodic)>>
               (subject, velocity_field, geometry,
                directories, name, identifier, params));
            else
              throw std::runtime_error{
                std::string("Measurement type ") + name + ": "
                "Particle state must define periodicity" };
            break;
          }
          case Measure::Type::absorption_time:
          {
            output.emplace_back
            (std::make_unique<Output_absorption_time>
             (subject, velocity_field, geometry,
              directories, name, identifier, params));
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
      switch (EndCriterion::type(params.end_criterion))
      {
        case EndCriterion::Type::time:
        {
          end_criterion = std::make_unique<
            Criterion_time<Subject>>(params.end_value);
          break;
        }
        case EndCriterion::Type::time_max:
        {
          end_criterion = std::make_unique<
            Criterion_time<Subject>>(params.time_max);
          break;
        }
        case EndCriterion::Type::mass_below:
        {
          end_criterion = std::make_unique<
            Criterion_mass_below<Subject>>(params.end_value);
          break;
        }
        case EndCriterion::Type::mass_above:
        {
          end_criterion = std::make_unique<
            Criterion_mass_above<Subject>>(params.end_value);
          break;
        }
        case EndCriterion::Type::all_absorbed:
        {
          end_criterion = std::make_unique<
            Criterion_all_absorbed<Subject>>();
          break;
        }
        case EndCriterion::Type::one_absorbed:
        {
          end_criterion = std::make_unique<
            Criterion_one_absorbed<Subject>>();
          break;
        }
        case EndCriterion::Type::fraction_absorbed:
        {
          end_criterion = std::make_unique<
            Criterion_fraction_absorbed<Subject>>(params.end_value);
          break;
        }
        default:
          throw std::runtime_error{
            std::string("End criterion ")
            + params.end_criterion
            + " not supported" };
      }
    }
    
    /** \struct NextMeasureTime PTOF/Output.h "PTOF/Output.h"
     *  \brief Polymorphic object to handle next measurement value */
    struct NextMeasureTime
    {
      virtual ~NextMeasureTime()
      {}
      
      /** Get next measurement time. */
      double time() const
      { return next_time; }
      
      /** Prepare for next.  */
      void advance()
      {
        ++next_measure;
        set_next_time();
      }
      
    protected:
      /** Base constructor. */
      NextMeasureTime
      (Parameters const& params)
      : params{ params }
      , next_time{ params.time_min }
      {}
      
      Parameters const& params;      /**< Output parameters.         */
      std::size_t next_measure{ 0 }; /**< Index of next measurement. */
      double next_time;              /**< Time of next measurement.  */
      
      virtual void set_next_time() = 0;
    };
    
    struct NextMeasureTime_step final : NextMeasureTime
    {
      NextMeasureTime_step
      (Parameters const& params)
      : NextMeasureTime{ params }
      {}
      
      void set_next_time() override
      {
        NextMeasureTime::next_time = NextMeasureTime::params.time_min
          + NextMeasureTime::next_measure
          * NextMeasureTime::params.time_increment;
      }
    };
    
    struct NextMeasureTime_linear final : NextMeasureTime
    {
      NextMeasureTime_linear
      (Parameters const& params)
      : NextMeasureTime{ params }
      {}
      
      void set_next_time() override
      {
        double time = NextMeasureTime::params.time_min
          + NextMeasureTime::next_measure
          * NextMeasureTime::params.time_increment;
        NextMeasureTime::next_time > NextMeasureTime::params.time_max
          ? std::numeric_limits<double>::infinity()
          : time;
      }
    };
    
    struct NextMeasureTime_log final : NextMeasureTime
    {
      NextMeasureTime_log
      (Parameters const& params)
      : NextMeasureTime{ params }
      {}
      
      void set_next_time() override
      {
         double time = NextMeasureTime::params.time_min
           * std::pow(NextMeasureTime::params.time_increment,
                      NextMeasureTime::next_measure);
         NextMeasureTime::next_time > NextMeasureTime::params.time_max
           ? std::numeric_limits<double>::infinity()
           : time;
      }
    };
    
    /** \struct OutputTime PTOF/Output.h "PTOF/Output.h"
     *  \brief Polymorphic output handler for outputting time and some quantity. */
    struct OutputTime
    {
      virtual ~OutputTime()
      { output.close(); }
      virtual void operator()(double time) = 0;
      
    protected:
      OutputTime
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : subject{ subject }
      , velocity_field{ velocity_field }
      , geometry{ geometry }
      , output{ useful::open_write(directories.dir_output
        + "/" + output_name
        + ".dat") } 
      , delimiter{ delimiter }
      {
        output << std::setprecision(precision)
               << std::scientific;
      }
      
      Subject const& subject;
      VelocityField const& velocity_field;
      Geometry const& geometry;
      std::ofstream output;
      std::string delimiter;
    };
        
    /** \struct Output PTOF/Output.h "PTOF/Output.h"
     * \brief Polymorphic output handler for outputting some quantity. */
    struct Output
    {
      virtual ~Output()
      { output.close(); }
      virtual void operator()() = 0;
      
    protected:
      Output
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : subject{ subject }
      , velocity_field{ velocity_field }
      , geometry{ geometry }
      , output{ useful::open_write(directories.dir_output
        + "/" + output_name
        + ".dat") }
      , delimiter{ delimiter }
      {
        output << std::setprecision(precision)
               << std::scientific;
      }
      
      Subject const& subject;
      VelocityField const& velocity_field;
      Geometry const& geometry;
      std::ofstream output;
      std::string delimiter;
    };
    
    /** \struct Output_position final PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time, tags, positions, and masses. */
    struct Output_position final : OutputTime
    {
      Output_position
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time;
        for (auto const& part : OutputTime::subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed)
          {
            OutputTime::output << OutputTime::delimiter << state.tag;
            useful::print(OutputTime::output, state.position,
                          1, OutputTime::delimiter);
            OutputTime::output << OutputTime::delimiter << state.mass;
          }
        }
        OutputTime::output << "\n";
      }
    };
    
    /** \struct Output_position_mean final PTOF/Output.h "PTOF/Output.h"
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
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time;
        useful::print(OutputTime::output,
                      position_mean(OutputTime::subject, time),
                      1, OutputTime::delimiter);
        OutputTime::output << "\n";
      }
    };
        
    /** \struct Output_position_second_moment final PTOF/Output.h "PTOF/Output.h"
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
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time;
        useful::print(OutputTime::output,
                      position_second_moment(OutputTime::subject, time),
                      1, OutputTime::delimiter);
        OutputTime::output << "\n";
      }
    };
        
    /** \struct Output_position_variance final PTOF/Output.h "PTOF/Output.h"
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
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time;
        useful::print(OutputTime::output,
                      position_variance(OutputTime::subject, time),
                      1, OutputTime::delimiter);
        OutputTime::output << "\n";
      }
    };
        
    /** \struct Output_mass final PTOF/Output.h "PTOF/Output.h"
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
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time << OutputTime::delimiter
                           << mass(OutputTime::subject, time) << "\n";
      }
    };
        
     /** \struct Output_velocity final PTOF/Output.h "PTOF/Output.h"
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
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      {}

      void operator()(double time) override
      {
        OutputTime::output << time;
        for (auto const& part : OutputTime::subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed)
          {
            OutputTime::output << OutputTime::delimiter << state.tag;
            useful::print(OutputTime::output,
                          OutputTime::velocity_field(state),
                          1, OutputTime::delimiter);
          }
        }
        OutputTime::output << "\n";
      }
    };
        
     /** \struct Output_velocity_gradient final PTOF/Output.h "PTOF/Output.h"
      *  \brief  Output time and velocity gradient components. 
      *  
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
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      , velocity_gradient{
        Foam::IOobject{
          "gradU", geometry.mesh.time().timeName(),
          geometry.mesh, Foam::IOobject::MUST_READ,
          Foam::IOobject::NO_WRITE },
          geometry.mesh }
      {
        if (params.velocity_rescaling != 1.)
          velocity_gradient *= params.velocity_rescaling;
      }

      void operator()(double time) override
      {
        OutputTime::output << time;
        for (auto const& part : OutputTime::subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed)
          {
            OutputTime::output << OutputTime::delimiter << state.tag;
            auto cell = OutputTime::geometry.locator(part.state_new());
            if (cell == -1)
              useful::print(OutputTime::output,
                            std::vector<double>(Geometry::dim, 0.),
                            1, OutputTime::delimiter);
            else
            {
              auto const& grad = velocity_gradient[cell];
              if constexpr (Geometry::dim == 1)
              {
                OutputTime::output << OutputTime::delimiter << grad.xx();
              }
              if constexpr (Geometry::dim == 2)
              {
                OutputTime::output << OutputTime::delimiter << grad.xx();
                OutputTime::output << OutputTime::delimiter << grad.xy();
                OutputTime::output << OutputTime::delimiter << grad.yx();
                OutputTime::output << OutputTime::delimiter << grad.yy();
              }
              if constexpr (Geometry::dim == 3)
              {
                OutputTime::output << OutputTime::delimiter << grad.xx();
                OutputTime::output << OutputTime::delimiter << grad.xy();
                OutputTime::output << OutputTime::delimiter << grad.xz();
                OutputTime::output << OutputTime::delimiter << grad.yx();
                OutputTime::output << OutputTime::delimiter << grad.yy();
                OutputTime::output << OutputTime::delimiter << grad.yz();
                OutputTime::output << OutputTime::delimiter << grad.zx();
                OutputTime::output << OutputTime::delimiter << grad.zy();
                OutputTime::output << OutputTime::delimiter << grad.zz();
              }
            }
          }
        }
        OutputTime::output << "\n";
      }

      Foam::volTensorField velocity_gradient;
    };
        
    /** \struct Output_pressure final PTOF/Output.h "PTOF/Output.h"
      *  \brief Output time, tags, and pressures. */
    struct Output_pressure final : OutputTime
    {
      Output_pressure
      (Subject const& subject,
       VelocityField const& velocity_field,
       Geometry const& geometry,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      , pressure{
        { Foam::IOobject{
          "U", geometry.mesh.time().timeName(),
          geometry.mesh, Foam::IOobject::MUST_READ,
          Foam::IOobject::NO_WRITE },
          geometry.mesh },
        geometry.locator,
        FieldOptions::Check{} }
      {
        if (params.velocity_rescaling != 1.)
          pressure.rescale(params.velocity_rescaling);
      }

      void operator()(double time) override
      {
        OutputTime::output << time;
        for (auto const& part : OutputTime::subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed)
          {
            OutputTime::output << OutputTime::delimiter << state.tag;
            useful::print(OutputTime::output, pressure(state),
                          1, OutputTime::delimiter);
          }
        }
        OutputTime::output << "\n";
      }

      ptof::ScalarField_LinearInterpolation_OF
      <Foam::volScalarField, ptof::Locator_Cell const&, 1, 1> pressure;
    };
        
    /** \struct Output_position_periodic final PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time, tags, positions, and masses
     *  with position accounting for periodicity. */
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
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      , getter_position{ geometry.boundary_periodic }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time;
        for (auto const& part : OutputTime::subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed)
          {
            OutputTime::output << OutputTime::delimiter << state.tag;
            useful::print(OutputTime::output, getter_position(state),
                          1, OutputTime::delimiter);
            OutputTime::output << OutputTime::delimiter << state.mass;
          }
        }
        OutputTime::output << "\n";
      }
      
      ctrw::Get_position_periodic<Boundary const&> getter_position;
    };
    
    /** \struct Output_position_mean_periodic final PTOF/Output.h "PTOF/Output.h"
     *  \brief Output time and mean position (weighted by mass)
     * with position accounting for periodicity. */
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
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      , getter_position{ geometry.boundary_periodic }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time;
        useful::print(OutputTime::output,
                      position_mean(OutputTime::subject, time,
                                    getter_position),
                      1, OutputTime::delimiter);
        OutputTime::output << "\n";
      }
      
      ctrw::Get_position_periodic<Boundary const&> getter_position;
    };
        
    /** \struct Output_position_second_moment_periodic final PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and second moment of position (weighted by mass)
     * with position accounting for periodicity. */
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
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      , getter_position{ geometry.boundary_periodic }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time;
        useful::print(OutputTime::output,
                      position_second_moment(OutputTime::subject, time,
                                             getter_position),
                      1, OutputTime::delimiter);
        OutputTime::output << "\n";
      }
      
      ctrw::Get_position_periodic<Boundary const&> getter_position;
    };
        
    /** \struct Output_position_variance_periodic final PTOF/Output.h "PTOF/Output.h"
     *  \brief  Output time and position variance (weighted by mass)
     * with position accounting for periodicity. */
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
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      , getter_position{ geometry.boundary_periodic }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time;
        useful::print(OutputTime::output,
                      position_variance(OutputTime::subject, time,
                                        getter_position),
                      1, OutputTime::delimiter);
        OutputTime::output << "\n";
      }
      
      ctrw::Get_position_periodic<Boundary const&> getter_position;
    };
        
    /** \struct Output_absorption_time final PTOF/Output.h "PTOF/Output.h"
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
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : Output{ subject, velocity_field,
        geometry, directories,
        output_name, identifier, params,
        precision, delimiter }
      {}
      
      void operator()() override
      {
        for (auto const& part : Output::subject.particles())
        {
          auto const& state = part.state_new();
          if (state.info.absorbed)
          {
            Output::output << state.time << Output::delimiter;
            Output::output << state.tag << Output::delimiter;
            Output::output << state.mass << "\n";
          }
        }
      }
    };
    
    Subject const& subject;                              /**< Subject to measure.                   */
    VelocityField const& velocity_field;                 /**< Compute velocity from particle state. */
    Geometry const& geometry;                            /**< Mesh and other domain information.    */
    std::unique_ptr<Criterion<Subject>> end_criterion;   /**< To check if end criterion is met.     */
    std::unique_ptr<NextMeasureTime> next_measure;       /**< Handle next measure time.             */
    std::vector<std::unique_ptr<OutputTime>> output_time;/**< Handle each output type given time.   */
    std::vector<std::unique_ptr<Output>> output;         /**< Handle each output type given nothing.*/
  };
}
        
#endif /* PTOF_OUTPUT_H */
