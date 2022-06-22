//
//  Output.h
//
//  Created by Tomás Aquino on 07/03/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Output_OF_h
#define Output_OF_h

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
#include "general/useful.h"
#include "ParticleTrackingOF/Directories.h"

namespace ptof
{
  // Compute number of absorbed particles
  template <typename Subject>
  auto nr_absorbed(Subject const& subject)
  {
    std::size_t absorbed = 0;
    for (auto const& part : subject.particles())
      if (part.state_new().info.absorbed)
        ++absorbed;
    return absorbed;
  };
  
  // Compute total mass
  template <typename Subject>
  auto mass(Subject const& subject)
  {
    double mass = 0.;
    for (auto const& part : subject.particles())
      if (!part.state_new().info.absorbed)
        mass += part.state_new().mass;
    return mass;
  };
  
  // Compute mean position (weighted by mass)
  template <typename Subject>
  auto position_mean(Subject const& subject)
  {
    decltype(subject.particles(0).state_new().position)
      position_mean = Foam::zero{};
    for (auto const& part : subject.particles())
    {
      auto const& state = part.state_new();
      if (!state.info.absorbed)
        position_mean += state.mass*state.position;
    }
    return position_mean/mass(subject);
  };
  
  // Compute second moment of position (weighted by mass)
  template <typename Subject>
  auto position_second_moment(Subject const& subject)
  {
    decltype(subject.particles(0).state_new().position)
      second_moment = Foam::zero{};
    for (auto const& part : subject.particles())
    {
      auto const& state = part.state_new();
      if (!state.info.absorbed)
      {
        for (std::size_t dd = 0; dd < second_moment.size(); ++dd)
          second_moment[dd] += state.position[dd]*state.position[dd];
        second_moment *= state.mass;
      }
    }
    return second_moment/mass(subject);
  };
  
  // Compute position variance (weighted by mass)
  template <typename Subject>
  auto position_variance(Subject const& subject)
  {
    auto position_mean_sq = position_mean(subject);
    for (std::size_t dd = 0; dd < position_mean_sq.size(); ++dd)
      position_mean_sq[dd] *= position_mean_sq[dd];
    return (position_second_moment(subject)-position_mean_sq)/mass(subject);
  };
  
  // Keep track of names and
  // types of end criteria
  struct EndCriterium
  {
    // Implemented types
    enum class Type
    {
      time,
      time_max,
      mass_below,
      mass_above,
      all_absorbed,
      one_absorbed,
      fraction_absorbed
    };
    
    // Type from name
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    // Name from type
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    // Check if name exists
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    // Map names to types
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
    
    // Map types to names
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
  
  // Keep track of names and
  // types of measure spacing
  struct MeasureSpacing
  {
    // Implemented types
    enum class Type
    {
      linear,
      log,
      step
    };
    
    // Type from name
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    // Name from type
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    // Check if name exists
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    // Map names to types
    inline static const
    std::unordered_map<std::string, Type> name_to_type
    {
      { "linear", Type::linear },
      { "log", Type::log },
      { "step", Type::step }
    };
    
    // Map types to names
    inline static const
    std::unordered_map<Type, std::string> type_to_name
    {
      { Type::linear, "linear" },
      { Type::log, "log" },
      { Type::step, "step" }
    };
  };
  
  // Keep track of names and
  // types of measure spacing units (scaling)
  struct MeasureSpacingUnits
  {
    // Implemented types
    enum class Type
    {
      diffusion,
      advection,
      reaction,
      arbitrary
    };
    
    // Type from name
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    // Name from type
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    // Check if name exists
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    // Map names to types
    inline static const
    std::unordered_map<std::string, Type>
    name_to_type
    {
      { "diffusion", Type::diffusion },
      { "advection", Type::advection },
      { "reaction", Type::reaction },
      { "arbitrary", Type::arbitrary }
    };
    
    // Map types to names
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
  
  // Keep track of measurement types
  struct Measure
  {
    // Implemented types
    enum class Type
    {
      position,
      position_mean,
      position_second_moment,
      position_variance,
      mass,
      absorption_time
    };
    
    // Type from name
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    // Name from type
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    // Check if name exists
    static bool exists(std::string const& name)
    { return name_to_type.count(name); }
    
    // Map names to types
    inline static const
    std::unordered_map<std::string, Type>
    name_to_type
    {
      { "position", Type::position },
      { "position_mean", Type::position_mean },
      { "position_second_moment", Type::position_second_moment },
      { "position_variance", Type::position_variance },
      { "mass", Type::mass },
      { "absorption_time", Type::absorption_time }
    };
    
    // Map types to names
    inline static const
    std::unordered_map<Type, std::string>
    type_to_name
    {
      { Type::position, "position" },
      { Type::position_mean, "position_mean" },
      { Type::position_second_moment, "position_second_moment" },
      { Type::position_variance, "position_variance" },
      { Type::mass, "mass" },
      { Type::absorption_time, "absorption_time" }
    };
  };
  
  // Polymorphic criteria
  template <typename Subject>
  struct Criterium
  {
    virtual bool operator()
    (Subject const& subject, double time) const;
    
    virtual ~Criterium()
    {}
    
  protected:
    Criterium()
    {}
  };
  
  // Check if time is greater than value
  template <typename Subject>
  struct Criterium_time final : Criterium<Subject>
  {
    Criterium_time(double end_value)
    : end_value{ end_value }
    {}
    
    bool operator()
    (Subject const& subject, double time) const override
    {
      return time > end_value;
    }
      
    double end_value;
  };
  
  // Check if mass is less than or equal to value
  template <typename Subject>
  struct Criterium_mass_below final : Criterium<Subject>
  {
    Criterium_mass_below(double end_value)
    : end_value{ end_value }
    {}
    
    bool operator()
    (Subject const& subject, double time) const override
    {
      return mass(subject) <= end_value;
    }
      
    double end_value;
  };
  
  // Check if mass is greater than or equal to value
  template <typename Subject>
  struct Criterium_mass_above final : Criterium<Subject>
  {
    Criterium_mass_above(double end_value)
    : end_value{ end_value }
    {}
    
    bool operator()
    (Subject const& subject, double time) const override
    {
      return mass(subject) >= end_value;
    }
      
    double end_value;
  };
  
  // Check if all particles have been absorbed
  template <typename Subject>
  struct Criterium_all_absorbed final : Criterium<Subject>
  {
    Criterium_all_absorbed()
    {}
    
    bool operator()
    (Subject const& subject, double time) const override
    {
      return nr_absorbed(subject) == subject.size();
    }
  };
  
  // Check if at least one particle has been absorbed
  template <typename Subject>
  struct Criterium_one_absorbed final : Criterium<Subject>
  {
    Criterium_one_absorbed()
    {}
    
    bool operator()(Subject const& subject, double time) const
    {
      return nr_absorbed(subject) > 0;
    }
  };
  
  // Check if at a least a given fraction of particles have been absorbed
  template <typename Subject>
  struct Criterium_fraction_absorbed final : Criterium<Subject>
  {
    Criterium_fraction_absorbed(double end_value)
    : end_value{ end_value }
    {}
    
    bool operator()
    (Subject const& subject, double time) const override
    {
      return nr_absorbed(subject) >= end_value*subject.size();
    }
      
    double end_value;
  };
  
  // Output object to handle implemented output options
  template <typename Subject>
  class Output_Cases
  {
  public:
    // Output parameters
    struct Parameters
    {
      std::size_t run_nr;
      std::string end_criterium;
      double end_value;
      std::string time_units;
      double time_unit_factor;
      std::string measure_spacing;
      double time_min;
      double time_max;
      double time_increment;
      std::vector<std::string> measure_names;
      
      // Get Output parameters from Directories info,
      // parameter set name, transport parameters,
      // reaction parameters, and solver parameters
      template
      <typename TransportParameters,
      typename ReactionParameters,
      typename SolverParameters>
      Parameters
      (Directories const& directories,
       std::string const& name,
       TransportParameters const& params_transport,
       ReactionParameters const& params_reaction,
       SolverParameters const& params_solvers)
      {
        auto input = useful::open_read(directories.dir_parameters
                                       + "/parameters_output_"
                                       + name + ".dat");
        useful::read(input, run_nr);
        read_end_criterium(input);
        read_measure_spacing(input,
                             params_transport, params_reaction);
        read_output_types(input);
        input.close();
      }
      
      // Output generic information about object
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Output parameters\n"
          "--------------------------------------------------\n"
          "- Run number (nonnegative integer to index output)\n"
          "- End criterium to finish dynamics:\n"
          "\ttime: Specified time\n"
          "\ttime_max: Maximum output time\n"
          "\tmass_below: Total mass below value\n"
          "\tmass_above: Total mass above value\n"
          "\tall_absorbed: All particles absorbed\n"
          "\tone_absorbed: One particle absorbed\n"
          "\tfraction_absorbed: Fraction of particles absorbed\n"
          "- End value, if required by end criterium\n"
          "- Time units to rescale measurement times:\n"
          "\tdiffusion: Rescale by diffusion time\n"
          "\tadvection: Rescale by reaction time\n"
          "\treaction: Rescale by reaction time\n"
          "\tuncscaled: Do not rescale\n"
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
      // Read end criterium from input stream
      template <typename IStream>
      void read_end_criterium(IStream& input)
      {
        useful::read(input, end_criterium);
        switch (EndCriterium::type(end_criterium))
        {
          case EndCriterium::Type::time:
          {
            useful::read(input, end_value);
            break;
          }
          case EndCriterium::Type::time_max:
            break;
          case EndCriterium::Type::mass_below:
          {
            useful::read(input, end_value);
            break;
          }
          case EndCriterium::Type::mass_above:
          {
            useful::read(input, end_value);
            break;
          }
          case EndCriterium::Type::all_absorbed:
            break;
          case EndCriterium::Type::one_absorbed:
            break;
          case EndCriterium::Type::fraction_absorbed:
          {
            useful::read(input, end_value);
            break;
          }
          default:
            throw std::runtime_error{
              std::string("End criterium ")
              + end_criterium
              + " not supported" };
        }
      }
      
      // Read measure spacing time units from input stream
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
      
      // Read measure spacing type from input stream
      template <typename IStream,
      typename TransportParameters,
      typename ReactionParameters>
      void read_measure_spacing
      (IStream& input,
       TransportParameters const& params_transport,
       ReactionParameters const& params_reaction)
      {
        read_time_units(input, params_transport, params_reaction);
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
      
      // Read measure spacing type from input stream
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
    
    // Construct given directory information,
    // output parameters,
    // and identifier for output filenames
    // Optionaly set output precision and delimiter character
    Output_Cases
    (Subject const& subject,
     Directories const& directories,
     Parameters params,
     std::string const& identifier,
     int precision = 8, std::string delimiter = "\t")
    : params{ params }
    , subject{ subject }
    {
      set_measure_types(directories,
                        identifier,
                        precision,
                        delimiter);
      set_end_criterium();
      set_next_measure_time();
    }
    
    // Check if end simulation criterium is satisfied
    bool done(double time) const
    {
      return end_criterium->operator()(subject, time);
    }

    // Set time of next measurement
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
            std::string("Measure spacing ")
            + params.measure_spacing
            + " not supported" };
      }
    }
    
    // Get time of next measurement
    double next_measure_time() const
    {
      return next_measure->time();
    }
        
    // Output requested measurements at given time
    // Advance to next measurement
    void operator()(double time)
    {
      for (auto const& measure : output_time)
        measure->operator()(time);
      next_measure->advance();
    }
    
    // Output current information
    void operator()()
    {
      for (auto const& measure : output)
        measure->operator()();
    }
    
    const Parameters params;  // Output parameters
    
    // Output information about current object
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Output\n"
        "--------------------------------------------------\n"
        "- Run number: " << params.run_nr << "\n"
        "- End criterium: " << params.end_criterium << "\n"
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
    // Open output file for a given output type
    std::ofstream open_write
    (Directories const& directories,
     std::string const& output_name,
     std::string const& identifier)
    {
      return
        useful::open_write(directories.dir_output
                           + "/Data"
                           + "_" + output_name
                           + "_" + identifier
                           + "_RUN_" + std::to_string(params.run_nr)
                           + ".dat");
    }
  
    // Set up output streams for requested output types
    void set_measure_types
    (Directories const& directories,
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
            (std::make_unique<
              Output_position>(subject,
                               directories,
                               "position",
                               identifier,
                               params));
            break;
          }
          case Measure::Type::position_mean:
          {
            output_time.emplace_back
            (std::make_unique<
              Output_position_mean>(subject,
                                    directories,
                                    "position_mean",
                                    identifier,
                                    params));
            break;
          }
          case Measure::Type::position_second_moment:
          {
            output_time.emplace_back
            (std::make_unique<
              Output_position_second_moment>(subject,
                                             directories,
                                             "position_second_moment",
                                             identifier,
                                             params));
            break;
          }
          case Measure::Type::position_variance:
          {
            output_time.emplace_back
            (std::make_unique<
              Output_position_variance>(subject,
                                        directories,
                                        "position_variance",
                                        identifier,
                                        params));
            break;
          }
          case Measure::Type::mass:
          {
            output_time.emplace_back
            (std::make_unique<
              Output_mass>(subject,
                           directories,
                           "mass",
                           identifier,
                           params));
            break;
          }
          case Measure::Type::absorption_time:
          {
            output.emplace_back
            (std::make_unique<
              Output_absorption_time>(subject,
                                      directories,
                                      "absorption_time",
                                      identifier,
                                      params));
            break;
          }
          default:
            throw std::runtime_error{
              std::string("Measurement type ")
              + name
              + " not supported" };
        }
      }
    }
    
    // Set end criterium
    void set_end_criterium()
    {
      switch (EndCriterium::type(params.end_criterium))
      {
        case EndCriterium::Type::time:
        {
          end_criterium = std::make_unique<
            Criterium_time<Subject>>(params.end_value);
          break;
        }
        case EndCriterium::Type::time_max:
        {
          end_criterium = std::make_unique<
            Criterium_time<Subject>>(params.time_max);
          break;
        }
        case EndCriterium::Type::mass_below:
        {
          end_criterium = std::make_unique<
            Criterium_mass_below<Subject>>(params.end_value);
          break;
        }
        case EndCriterium::Type::mass_above:
        {
          end_criterium = std::make_unique<
            Criterium_mass_above<Subject>>(params.end_value);
          break;
        }
        case EndCriterium::Type::all_absorbed:
        {
          end_criterium = std::make_unique<
            Criterium_all_absorbed<Subject>>();
          break;
        }
        case EndCriterium::Type::one_absorbed:
        {
          end_criterium = std::make_unique<
            Criterium_one_absorbed<Subject>>();
          break;
        }
        case EndCriterium::Type::fraction_absorbed:
        {
          end_criterium = std::make_unique<
            Criterium_fraction_absorbed<Subject>>(params.end_value);
          break;
        }
        default:
          throw std::runtime_error{
            std::string("End criterium ")
            + params.end_criterium
            + " not supported" };
      }
    }
    
    // Polymorphic object to handle next measurement value
    struct NextMeasureTime
    {
      virtual ~NextMeasureTime()
      {}
      
      // Get next measurement time
      double time() const
      { return next_time; }
      
      // Prepare for next measure
      void advance()
      {
        ++next_measure;
        set_next_time();
      }
      
    protected:
      // Base constructor
      NextMeasureTime
      (Parameters const& params)
      : params{ params }
      , next_time{ params.time_min }
      {}
      
      Parameters const& params; // Output parameters
      std::size_t next_measure{ 0 }; // Index of next measurement
      double next_time; // Time of next measurement
      
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
    
    // Polymorphic output handler for outputting
    // time and some quantity
    struct OutputTime
    {
      virtual ~OutputTime()
      { output.close(); }
      virtual void operator()(double time) = 0;
      
    protected:
      OutputTime
      (Subject const& subject,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : subject{ subject }
      , output{ useful::open_write(directories.dir_output
        + "/Data"
        + "_" + output_name
        + "_" + identifier
        + "_RUN_" + std::to_string(params.run_nr)
        + ".dat") }
      , delimiter{ delimiter }
      {
        output << std::setprecision(precision)
               << std::scientific;
      }
      
      Subject const& subject;
      std::ofstream output;
      std::string delimiter;
    };
        
    // Polymorphic output handler for outputting some quantity
    struct Output
    {
      virtual ~Output()
      { output.close(); }
      virtual void operator()() = 0;
      
    protected:
      Output
      (Subject const& subject,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : subject{ subject }
      , output{ useful::open_write(directories.dir_output
        + "/Data"
        + "_" + output_name
        + "_" + identifier
        + "_RUN_" + std::to_string(params.run_nr)
        + ".dat") }
      , delimiter{ delimiter }
      {
        output << std::setprecision(precision)
               << std::scientific;
      }
      
      Subject const& subject;
      std::ofstream output;
      std::string delimiter;
    };
    
    // Output time, positions, and masses
    struct Output_position final : OutputTime
    {
      Output_position
      (Subject const& subject,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, directories, output_name,
        identifier, params,
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
            useful::print(OutputTime::output, state.position,
                          1, OutputTime::delimiter);
            OutputTime::output << OutputTime::delimiter << state.mass;
          }
        }
        OutputTime::output << "\n";
      }
    };
    
    // Output time and mean position (weighted by mass)
    struct Output_position_mean final : OutputTime
    {
      Output_position_mean
      (Subject const& subject,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, directories, output_name,
        identifier, params,
        precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time;
        useful::print(OutputTime::output,
                      position_mean(OutputTime::subject),
                      1, OutputTime::delimiter);
        OutputTime::output << "\n";
      }
    };
        
    // Output time and second moment of position (weighted by mass)
    struct Output_position_second_moment final : OutputTime
    {
      Output_position_second_moment
      (Subject const& subject,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, directories, output_name,
        identifier, params,
        precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time;
        useful::print(OutputTime::output,
                      position_second_moment(OutputTime::subject),
                      1, OutputTime::delimiter);
        OutputTime::output << "\n";
      }
    };
        
    // Output time and position variance (weighted by mass)
    struct Output_position_variance final : OutputTime
    {
      Output_position_variance
      (Subject const& subject,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, directories, output_name,
        identifier, params,
        precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time;
        useful::print(OutputTime::output,
                      position_variance(OutputTime::subject),
                      1, OutputTime::delimiter);
        OutputTime::output << "\n";
      }
    };
        
    // Output time and total mass
    struct Output_mass final : OutputTime
    {
      Output_mass
      (Subject const& subject,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : OutputTime{ subject, directories, output_name,
        identifier, params,
        precision, delimiter }
      {}
      
      void operator()(double time) override
      {
        OutputTime::output << time << OutputTime::delimiter
                           << mass(OutputTime::subject) << "\n";
      }
    };
        
    // Output absorption times and masses of absorbed particles
    struct Output_absorption_time final : Output
    {
      Output_absorption_time
      (Subject const& subject,
       Directories const& directories,
       std::string const& output_name,
       std::string const& identifier,
       Parameters const& params,
       int precision = 8,
       std::string delimiter = "\t")
      : Output{ subject, directories, output_name,
        identifier, params,
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
            Output::output << state.mass << "\n";
          }
        }
      }
    };
    
    Subject const& subject; // Subject to measure
    std::unique_ptr<Criterium<Subject>> end_criterium;  // To check if end criterium is met
    std::unique_ptr<NextMeasureTime> next_measure;      // Handle next measure time
    std::vector<std::unique_ptr<OutputTime>> output_time; // Handle each output type given time
    std::vector<std::unique_ptr<Output>> output; // Handle each output type given nothing
  };
}
        
#endif /* Output_OF_h */
