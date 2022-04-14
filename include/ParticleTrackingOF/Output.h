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
  template <typename Subject>
  auto nr_absorbed(Subject const& subject)
  {
    std::size_t absorbed = 0;
    for (auto const& part : subject.particles())
      if (part.state_new().info.absorbed)
        ++absorbed;
    return absorbed;
  };
  
  template <typename Subject>
  auto mass(Subject const& subject)
  {
    double mass = 0.;
    for (auto const& part : subject.particles())
      if (!part.state_new().info.absorbed)
        mass += part.state_new().mass;
    return mass;
  };
  
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
  
  template <typename Subject>
  auto position_variance(Subject const& subject)
  {
    auto position_mean_sq = position_mean(subject);
    for (std::size_t dd = 0; dd < position_mean_sq.size(); ++dd)
      position_mean_sq[dd] *= position_mean_sq[dd];
    return (position_second_moment(subject)-position_mean_sq)/mass(subject);
  };
  
  // Output object to handle implemented output options
  struct Output_Cases
  {
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
      std::vector<std::string> to_measure;
      
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
        
        // Get run number identifier
        useful::read(input, run_nr);
        
        // Set end simulation criterium
        useful::read(input, end_criterium);
        if (end_criterium == "time")
        {
          useful::read(input, end_value);
        }
        else if (end_criterium == "time_max")
        {}
        else if (end_criterium == "mass_below")
        {
          useful::read(input, end_value);
        }
        else if (end_criterium == "mass_above")
        {
          useful::read(input, end_value);
        }
        else if (end_criterium == "all_absorbed")
        {}
        else if (end_criterium == "one_absorbed")
        {}
        else if (end_criterium == "fraction_absorbed")
        {
          useful::read(input, end_value);
        }
        else
        throw std::runtime_error{
          std::string("End criterium ")
          + end_criterium
          + " not supported"
        };
        
        // Set output time units
        useful::read(input, time_units);
        time_unit_factor = 1.;
        if (time_units == "diffusion")
        {
          time_unit_factor = params_transport.diffusion_time;
        }
        else if (time_units == "advection")
        {
          time_unit_factor = params_transport.advection_time;
        }
        else if (time_units == "reaction")
        {
          time_unit_factor = params_reaction.reaction_time;
        }
        else if (time_units == "unscaled")
        {}
        else
          throw std::runtime_error{
            std::string("Time units ")
            + time_units
            + " not supported"
          };
        
        // Set spacing type mesurements
        useful::read(input, measure_spacing);
        useful::read(input, time_min);
        time_min *= time_unit_factor;
        if (measure_spacing == "step")
        {
          useful::read(input, time_increment);
          time_increment *= time_unit_factor;
          time_max = std::numeric_limits<double>::infinity();
        }
        else if (measure_spacing == "linear")
        {
          useful::read(input, time_max);
          time_max *= time_unit_factor;
          std::size_t nr_measures;
          useful::read(input, nr_measures);
          time_increment = (time_max - time_min)/(nr_measures - 1);
        }
        else if (measure_spacing == "log")
        {
          useful::read(input, time_max);
          time_max *= time_unit_factor;
          std::size_t nr_measures;
          useful::read(input, nr_measures);
          time_increment = std::pow(time_max/time_min, 1./(nr_measures - 1));
        }
        else
          throw std::runtime_error{
            std::string("Measure spacing ")
            + measure_spacing
            + " not supported"
          };
        
        // Get requested output types
        while (1)
        {
          std::string measurement;
          input >> measurement;
          if (!input.fail())
            to_measure.push_back(measurement);
          else
            break;
        }
        
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
    };
    
    // Construct given directory information,
    // output parameters,
    // and identifier for output filenames
    // Optionaly set output precision and delimiter character
    Output_Cases
    (Directories const& directories,
     Parameters params,
     std::string const& identifier,
     int precision = 8, std::string delimiter = "\t")
    : params{ params }
    , delimiter{ delimiter }
    {
      // Set up output streams for requested output types
      if (std::any_of(params.to_measure.begin(), params.to_measure.end(),
                      [](std::string const& string)
                      { return string == "position"; }))
      {
        output_position = open_write(directories, identifier,
                                     "position");
        output_position << std::setprecision(precision)
                        << std::scientific;
      }
      if (std::any_of(params.to_measure.begin(), params.to_measure.end(),
                      [](std::string const& string)
                      { return string == "position_mean"; }))
      {
        output_position_mean = open_write(directories, identifier,
                                          "position_mean");
        output_position_mean << std::setprecision(precision)
                             << std::scientific;
      }
      if (std::any_of(params.to_measure.begin(), params.to_measure.end(),
                      [](std::string const& string)
                      { return string == "position_second_moment"; }))
      {
        output_position_second_moment = open_write(directories, identifier,
                                                   "position_second_moment");
        output_position_second_moment << std::setprecision(precision)
                                      << std::scientific;
      }
      if (std::any_of(params.to_measure.begin(), params.to_measure.end(),
                      [](std::string const& string)
                      { return string == "position_variance"; }))
      {
        output_position_variance = open_write(directories, identifier,
                                              "position_variance");
        output_position_variance << std::setprecision(precision)
                                 << std::scientific;
      }
      if (std::any_of(params.to_measure.begin(), params.to_measure.end(),
                      [](std::string const& string)
                      { return string == "absorption_time"; }))
      {
        output_absorption_time = open_write(directories, identifier,
                                            "absorption_time");
        output_absorption_time << std::setprecision(precision)
                               << std::scientific;
      }
      if (std::any_of(params.to_measure.begin(), params.to_measure.end(),
                      [](std::string const& string)
                      { return string == "mass"; }))
      {
        output_mass = open_write(directories, identifier,
                                 "mass");
        output_mass << std::setprecision(precision)
                               << std::scientific;
      }
    }
    
    // Destructor: close any open output streams
    ~Output_Cases()
    {
      if (output_position.is_open())
        output_position.close();
      if (output_position_mean.is_open())
        output_position_mean.close();
      if (output_position_second_moment.is_open())
        output_position_second_moment.close();
      if (output_position_variance.is_open())
        output_position_variance.close();
      if (output_absorption_time.is_open())
        output_absorption_time.close();
      if (output_mass.is_open())
        output_mass.close();
    }
    
    // Check if end simulation criterium is satisfied
    template <typename Subject>
    bool done(Subject const& subject, double time) const
    {
      if (params.end_criterium == "time")
      {
        return time > params.end_value;
      }
      else if (params.end_criterium == "time_max")
      {
        return time > params.time_max;
      }
      else if (params.end_criterium == "mass_below")
      {
       return mass(subject) <= params.end_value;
      }
      else if (params.end_criterium == "mass_above")
      {
        return mass(subject) >= params.end_value;
      }
      else if (params.end_criterium == "all_absorbed")
      {
        return nr_absorbed(subject) == subject.size();
      }
      else if (params.end_criterium == "one_absorbed")
      {
        return nr_absorbed(subject) > 0;
      }
      else if (params.end_criterium == "fraction_absorbed")
      {
        return nr_absorbed(subject) >= params.end_value*subject.size();
      }
      else
        throw std::runtime_error{
          std::string("End criterium ")
          + params.end_criterium
          + " not supported"
      };
    }

    // Get time of next measurement
    double next_measure_time() const
    {
      if (params.measure_spacing == "step")
      {
        return params.time_min + next_measure*params.time_increment;
      }
      else if (params.measure_spacing == "linear")
      {
        double next_time = params.time_min*
          std::pow(params.time_increment, next_measure);
        return next_time > params.time_max
          ? std::numeric_limits<double>::infinity()
          : next_time;
      }
      else if (params.measure_spacing == "log")
      {
        double next_time = params.time_min
          + next_measure*params.time_increment;
        return next_time > params.time_max
          ? std::numeric_limits<double>::infinity()
          : next_time;
      }
      else
        throw std::runtime_error{
          std::string("Measure spacing ")
          + params.measure_spacing
          + " not supported"
        };
    }
        
    // Output requested measurements at given time
    // Advance to next measurement
    template <typename Subject>
    void operator()(Subject const& subject, double time)
    {
      ++next_measure;
      
      if (output_position.is_open())
      {
        output_position << time;
        for (auto const& part : subject.particles())
        {
          auto const& state = part.state_new();
          if (!state.info.absorbed)
          {
            useful::print(output_position, state.position, 1, delimiter);
            output_position << delimiter << state.mass;
          }
        }
        output_position << "\n";
      }
      if (output_position_mean.is_open())
      {
        output_position_mean << time;
        useful::print(output_position_mean,
                      position_mean(subject), 1, delimiter);
        output_position_mean << "\n";
      }
      if (output_position_second_moment.is_open())
      {
        output_position_second_moment << time;
        useful::print(output_position_second_moment,
                      position_second_moment(subject), 1, delimiter);
        output_position_second_moment << "\n";
      }
      if (output_position_variance.is_open())
      {
        output_position_variance << time;
        useful::print(output_position_variance,
                      position_variance(subject),
                      1, delimiter);
        output_position_variance << "\n";
      }
      if (output_mass.is_open())
      {
        output_mass << time << delimiter
                    << mass(subject) << "\n";
      }
    }
    
    // Output current information
    template <typename Subject>
    void operator()(Subject const& subject)
    {
      if (output_absorption_time.is_open())
      {
        for (auto const& part : subject.particles())
        {
          auto const& state = part.state_new();
          if (state.info.absorbed)
          {
            output_absorption_time << state.time << delimiter;
            output_absorption_time << state.mass << "\n";
          }
        }
      }
    }
    
    Parameters params;              // Output parameters
    std::string delimiter;          // Delimiter character for output formatting
    std::size_t next_measure{ 0 };  // Number of next measurement
    
    // Output streams for different output types
    std::ofstream output_position;
    std::ofstream output_position_mean;
    std::ofstream output_position_second_moment;
    std::ofstream output_position_variance;
    std::ofstream output_absorption_time;
    std::ofstream output_mass;
        
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
        "- Measurement types:\n";
      useful::print(output, params.to_measure, 1);
      output << "\n";
      output <<
        "--------------------------------------------------\n";
    }
        
    private:
      // Open output file for a given output type
      std::ofstream open_write
      (Directories const& directories,
       std::string const& identifier,
       std::string const& output_name)
      {
        return
          useful::open_write(directories.dir_output
                             + "/Data"
                             + "_" + output_name
                             + "_" + identifier
                             + "_RUN_" + std::to_string(params.run_nr)
                             + ".dat");
      }
  };
}

#endif /* Output_OF_h */
