//
//  Directories.h
//
//  Created by Tomás Aquino on 22/02/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Directories_OF_h
#define Directories_OF_h

#include <string>
#include <Time.H>
#include "general/useful.h"

namespace ptof
{
  // Object to handle simulation case directories
  struct Directories
  {
    std::string dir_case;                 // Case directory
    std::string dir_output;               // Output directory
    std::string dir_parameters;           // Parameters directory
    std::string dir_boundaryconditions;   // Boundary conditions directoty
    
    // Construct given
    // base cases directory (default ../cases if empty string),
    // case name,
    // and output directory (default <case directory>/../cases if empty string)
    Directories
    (std::string const& dir,
     std::string const& case_name,
     std::string const& dir_output)
    : dir_case{
      (useful::empty(dir)
       ? "../cases"
       : dir) + "/" + case_name }
    , dir_output{
      useful::empty(dir_output)
      ? dir_case + "/output"
      : dir_output }
    , dir_parameters{ dir_case + "/parameters" }
    , dir_boundaryconditions{ dir_case + "/boundary_conditions" }
    {}
    
    // Output information about current object
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Directories\n"
        "--------------------------------------------------\n"
        "Case directory                : " + dir_case + "\n"
        "Parameters directory          : " + dir_parameters + "\n"
        "Boundary conditions directory : " + dir_boundaryconditions + "\n"
        "Output directory              : " + dir_output + "\n"
        "--------------------------------------------------\n";
    }
  };

  // Object to handle OpenFOAM case directories
  struct DirectoriesOF
  {
  private:
    std::string dir;              // OpenFOAM base cases directory
    std::string time_name_input;  // Information about time directory
    
  public:
    std::string case_name;        // OpenFOAM case name
    std::string dir_case;         // OpenFOAM case directory
    Foam::Time time;              // OpenFOAM case time
    
    // Construct given Directories information
    DirectoriesOF(Directories const& directories)
    : time{ makeRunTime(directories) }
    {
      if (useful::empty(time_name_input))
        time.setTime(time.times().last(), 0);
      else
        time.setTime(Foam::instant(std::stod(time_name_input)), 0);
    }
    
    // Output generic information about object
    template <typename OStream>
    static void info(OStream& output)
    {
      output <<
        "--------------------------------------------------\n"
        "OpenFOAM directories (pass '' for default in [])\n"
        "--------------------------------------------------\n"
        "OpenFOAM cases directory [$FOAM_RUN]\n"
        "OpenFOAM case name [$FOAM_DIR]\n"
        "OpenFOAM case time [last OpenFOAM case time]\n"
        "--------------------------------------------------\n";
    }

    // Output information about current object
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "OpenFOAM directories\n"
        "--------------------------------------------------\n"
        "Case directory: " + dir_case + "\n"
        "Time directory: " + time.timeName() + "\n"
        "--------------------------------------------------\n";
    }
    
  private:
    // Make OpenFOAM Time instance given Directories information
    Foam::Time makeRunTime(Directories const& directories)
    {
      auto input = useful::open_read(directories.dir_parameters
                                     + "/of.dat");
      useful::read(input, dir);
      useful::read(input, case_name);
      useful::read(input, time_name_input);
      input.close();
      
      // Use OpenFOAM
      if (useful::empty(dir))
        dir = useful::expand_env("${FOAM_RUN}");
      dir_case = dir + "/" + case_name;
      
      return { Foam::Time::controlDictName, dir, case_name };
    }
  };
}


#endif /* Directories_OF_h */
