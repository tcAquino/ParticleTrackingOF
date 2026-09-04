/**
 * @file   Directories.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Tue Feb 22 00:00:00 2022
 *
 * @brief Handle input and output directories.
 */

#ifndef PTOF_DIRECTORIES_H
#define PTOF_DIRECTORIES_H

#include "General/IO.h"
#include "General/Useful.h"
#include <Time.H>
#include <filesystem>
#include <ostream>
#include <stdexcept>
#include <string>

namespace ptof {

/** @brief Object to handle simulation case directories. */
struct Directories {
  const std::filesystem::path dir;            /**< Case directory. */
  const std::string case_name;                /**< Case name. */
  const std::filesystem::path dir_output;     /**< Output directory. */
  const std::filesystem::path dir_parameters; /**< Parameters directory. */
  const std::filesystem::path
      dir_boundaryconditions; /**< Boundary conditions directory. */

  /**
   * @brief Constructor.
   *
   * @param dir_case Case directory.
   *
   * @param dir_output Output directory relative to \p dir.
   */
  Directories(std::string const &dir)
      : dir{std::filesystem::canonical(io::expand_home_and_env(dir))},
        case_name{this->dir.filename()}, dir_output{this->dir / "output"},
        dir_parameters{this->dir / "parameters"},
        dir_boundaryconditions{this->dir / "boundary_conditions"} {}

  /** @brief Output information about current object. */
  std::ostream &info_runtime(std::ostream &output) const {
    output << io::line() << "Directories\n"
           << io::line() << "Case: " << dir
           << "\n"
              "Parameters: "
           << dir_parameters
           << "\n"
              "Boundary conditions: "
           << dir_boundaryconditions
           << "\n"
              "Output: "
           << dir_output << "\n"
           << io::line();
    return output;
  }
};

/** @brief Object to handle OpenFOAM case directories. */
struct DirectoriesOF {
private:
  std::string _start_time_name;

public:
  std::filesystem::path dir; /**< OpenFOAM case directory. */
  std::string case_name;     /**< OpenFOAM case name. */
  Foam::Time time;           /**< OpenFOAM case time. */

  /**
   * @brief Constructor.
   *
   * @param directories Current case directory information.
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

  /** @brief Deleted copy constructor. */
  DirectoriesOF(DirectoriesOF const &) = delete;

  /**
   * @brief Output generic information about object.
   *
   * @param output Output stream.
   */
  static std::ostream &info(std::ostream &output) {
    output << io::line() << "OpenFOAM directories\n"
           << io::line() << R"(
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
   * @brief Output information about current object.
   *
   * @param output Output stream.
   */
  std::ostream &info_runtime(std::ostream &output) const {
    output << io::line() << "OpenFOAM directories\n"
           << io::line() << "Case: " << dir
           << "\n"
              "Time: " +
                  time.timeName() + "\n"
           << io::line();
    return output;
  }

private:
  /**
   * @param directories Current case directory information.
   *
   * @return OpenFOAM runTime instance given Directories information.
   */
  Foam::Time makeRunTime(Directories const &directories) {
    // openFOAM makes modifications to std::cout when a runTime object is
    // constructed; this restricts them to the local scope
    io::StreamScopeFormat guard{std::cout};

    std::string filename = directories.dir_parameters / "of.param";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    dir = std::filesystem::canonical(
        io::expand_home_and_env(io::read<std::string>(
            split_line, param_index,
            in_file + "Could not parse base case directory")));
    case_name = dir.filename();

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

    return {Foam::Time::controlDictName, dir.parent_path().string(), case_name};
  }
};
} // namespace ptof

#endif /* PTOF_DIRECTORIES_H */
