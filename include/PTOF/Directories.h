/**
   \file PTOF/Directories.h
   \author Tomás Aquino
   \date 22/02/2022
   \brief Handle input and output directories.
*/

#ifndef PTOF_DIRECTORIES_H
#define PTOF_DIRECTORIES_H

#include "General/IO.h"
#include <Time.H>
#include <ostream>
#include <string>

namespace ptof {

/**
   \struct Directories PTOF/Directories.h "PTOF/Directories.h"
   \brief Object to handle simulation case directories.
*/
struct Directories {
  std::string dir_case;               /**< Case directory.                */
  std::string dir_output;             /**<  Output directory.              */
  std::string dir_parameters;         /**< Parameters directory.          */
  std::string dir_boundaryconditions; /**< Boundary conditions directory. */

  /**
     \brief Constructor.
     \param dir Path to current case directory.
     \param case_name Current case directory relative to \p dir.
     \param dir_output Output directory relative to \p dir.
  */
  Directories(std::string const &dir, std::string const &case_name,
              std::string const &dir_output)
      : dir_case{(io::is_empty(dir)
                      ? "../cases"
                      : io::expand_env(io::expand_home_dir(dir))) +
                 "/" + case_name},
        dir_output{io::is_empty(dir_output)
                       ? dir_case + "/output"
                       : io::expand_env(io::expand_home_dir(dir_output))},
        dir_parameters{dir_case + "/parameters"},
        dir_boundaryconditions{dir_case + "/boundary_conditions"} {}

  /** \brief Output information about current object. */
  std::ostream &info_runtime(std::ostream &output) const {
    output << io::line() << "Directories\n"
           << io::line()
           << "Case: " + dir_case +
                  "\n"
                  "Parameters: " +
                  dir_parameters +
                  "\n"
                  "Boundary conditions: " +
                  dir_boundaryconditions +
                  "\n"
                  "Output: " +
                  dir_output + "\n"
           << io::line();
    return output;
  }
};

/**
   \struct DirectoriesOF PTOF/Directories.h "PTOF/Directories.h"
   \brief Object to handle OpenFOAM case directories.
*/
struct DirectoriesOF {
  std::string _time_name; /**< Information about time directory. */

public:
  std::string case_name; /**< OpenFOAM case name.       */
  std::string dir_case;  /**< OpenFOAM case directory.  */
  Foam::Time time;       /**< OpenFOAM case time.       */

  /**
     \brief Constructor.
     \param directories Current case directory information.
  */
  DirectoriesOF(Directories const &directories)
      : time{makeRunTime(directories)} {
    if (io::is_empty(_time_name))
      time.setTime(time.times().last(), 0);
    else
      time.setTime(Foam::instant(std::stod(_time_name)), 0);
  }

  /** \brief Deleted copy constructor. */
  DirectoriesOF(DirectoriesOF const &) = delete;

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  static std::ostream &info(std::ostream &output) {
    output << io::line() << "OpenFOAM directories (pass '' for default in [])\n"
           << io::line()
           << "OpenFOAM cases directory [$FOAM_RUN]\n"
              "OpenFOAM case name\n"
              "OpenFOAM case time [last OpenFOAM case time]\n"
           << io::line();
    return output;
  }

  /**
     \brief Output information about current object.
     \param output Output stream.
  */
  std::ostream &info_runtime(std::ostream &output) const {
    output << io::line() << "OpenFOAM directories\n"
           << io::line()
           << "Case: " + dir_case +
                  "\n"
                  "Time: " +
                  time.timeName() + "\n"
           << io::line();
    return output;
  }

private:
  /**
     \param directories Current case directory information.
     \return OpenFOAM runTime instance given Directories information.
  */
  Foam::Time makeRunTime(Directories const &directories) {
    // openFOAM makes modifications to std::cout when a runTime object is
    // constructed; this restricts them to the local scope
    io::StreamScopeFormat guard{std::cout};

    std::string filename = directories.dir_parameters + "/of.dat";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    auto dir = io::read<std::string>(split_line, param_index,
                                     in_file + "Could not parse directory");
    if (io::is_empty(dir))
      dir = "${FOAM_RUN}";
    io::expand_env_in_place(io::expand_home_dir_in_place(dir));

    split_line = io::split_line(input);
    param_index = 0;
    io::read(split_line, param_index, in_file + "Could not parse case name",
             case_name);
    dir_case = dir + "/" + case_name;

    split_line = io::split_line(input);
    param_index = 0;
    io::read(split_line, param_index, in_file + "Could not parse time name",
             _time_name);

    return {Foam::Time::controlDictName, dir, case_name};
  }
};
} // namespace ptof

#endif /* PTOF_DIRECTORIES_H */
