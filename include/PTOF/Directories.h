/**
   \file PTOF/Directories.h
   \author Tomas Aquino
   \date 22/02/2022
   \brief Handle input and output directories.
*/

#ifndef PTOF_DIRECTORIES_H
#define PTOF_DIRECTORIES_H

#include "General/IO.h"
#include "General/Useful.h"
#include <Time.H>
#include <ostream>
#include <stdexcept>
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
      : dir_case{io::expand_env(io::expand_home_dir(
                     io::read_or_default(dir, std::string{"../cases"}, ""))) +
                 "/" + case_name},
        dir_output{io::expand_env(io::expand_home_dir(
            io::read_or_default(dir_output, dir_case + "/output", "")))},
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
private:
  std::string _start_time_name;

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
    if (time.times().size() == 0) {
      throw std::runtime_error{"No OpenFOAM time folders"};
    } else if (_start_time_name == "last") {
      time.setTime(time.times().back(), 0);
    } else if (_start_time_name == "start") {
      time.setTime(time.startTime(), 0);
    } else if (_start_time_name == "end") {
      time.setTime(time.endTime(), 0);
    } else if (useful::is_numeric(_start_time_name)) {
      time.setTime(Foam::instant(std::stod(_start_time_name)), 0);
    }
  }

  /** \brief Deleted copy constructor. */
  DirectoriesOF(DirectoriesOF const &) = delete;

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  static std::ostream &info(std::ostream &output) {
    output << io::line() << "OpenFOAM directories (pass '' for default in [])\n"
           << io::line() << R"(- OpenFOAM cases directory [${FOAM_RUN}].
- OpenFOAM case name.
- OpenFOAM case time:
  - last
    - Use OpenFOAM case last available time.
  - start
    - Use OpenFOAM case start time.
  - end
    - Use OpenFOAM case end time.
  - Numeric value
    - Use given numeric value as start time name.
)" << io::line();
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

    std::string filename = directories.dir_parameters + "/of.param";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    auto dir =
        io::read_or_default(split_line, param_index, std::string{"${FOAM_RUN}"},
                            in_file + "Could not parse directory");
    io::expand_env_in_place(io::expand_home_dir_in_place(dir));

    split_line = io::split_line(input);
    param_index = 0;
    io::read(split_line, param_index, in_file + "Could not parse case name",
             case_name);
    dir_case = dir + "/" + case_name;

    split_line = io::split_line(input);
    param_index = 0;
    io::read(split_line, param_index, in_file + "Could not parse start time",
             _start_time_name);
    if (_start_time_name != "last" && _start_time_name != "start" &&
        _start_time_name != "end" && !useful::is_numeric(_start_time_name)) {
      throw std::runtime_error{
          in_file +
          "Expected last, start, end, or a numeric value for start time, got " +
          _start_time_name};
    }

    return {Foam::Time::controlDictName, dir, case_name};
  }
};
} // namespace ptof

#endif /* PTOF_DIRECTORIES_H */
