/**
 \file PTOF/Directories.h
 \author Tomás Aquino
 \date 22/02/2022
*/

#ifndef PTOF_DIRECTORIES_H
#define PTOF_DIRECTORIES_H

#include "General/Useful.h"
#include <Time.H>
#include <string>

namespace ptof {

/** \struct Directories PTOF/Directories.h "PTOF/Directories.h"
 *  \brief Object to handle simulation case directories.*/
struct Directories {
  std::string dir_case;               /**< Case directory.                */
  std::string dir_output;             /**<  Output directory.              */
  std::string dir_parameters;         /**< Parameters directory.          */
  std::string dir_boundaryconditions; /**< Boundary conditions directory. */

  /** Constructor.
  \param dir Path to current case directory.
  \param case_name Current case directory relative to \p dir.
  \param dir_output Output directory relative to \p dir.
  */
  Directories(std::string const &dir, std::string const &case_name,
              std::string const &dir_output)
      : dir_case{(useful::is_empty(dir)
                      ? "../cases"
                      : useful::expand_env(useful::expand_home_dir(dir))) +
                 "/" + case_name},
        dir_output{
            useful::is_empty(dir_output)
                ? dir_case + "/output"
                : useful::expand_env(useful::expand_home_dir(dir_output))},
        dir_parameters{dir_case + "/parameters"},
        dir_boundaryconditions{dir_case + "/boundary_conditions"} {}

  /** \brief Output information about current object. */
  template <typename OStream> void info_runtime(OStream &output) const {
    output
        << "--------------------------------------------------------------\n"
           "Directories\n"
           "--------------------------------------------------------------\n"
           "Case directory                : " +
               dir_case +
               "\n"
               "Parameters directory          : " +
               dir_parameters +
               "\n"
               "Boundary conditions directory : " +
               dir_boundaryconditions +
               "\n"
               "Output directory              : " +
               dir_output + "\n"
        << "--------------------------------------------------------------\n";
  }
};

/** \struct DirectoriesOF PTOF/Directories.h "PTOF/Directories.h"
 *  \brief Object to handle OpenFOAM case directories.*/
struct DirectoriesOF {
  std::string _time_name_input; /**< Information about time directory. */

public:
  std::string case_name; /**< OpenFOAM case name.       */
  std::string dir_case;  /**< OpenFOAM case directory.  */
  Foam::Time time;       /**< OpenFOAM case time.       */

  /** Constructor
   \param directories Current case directory information.*/
  DirectoriesOF(Directories const &directories)
      : time{makeRunTime(directories)} {
    if (useful::is_empty(_time_name_input))
      time.setTime(time.times().last(), 0);
    else
      time.setTime(Foam::instant(std::stod(_time_name_input)), 0);
  }

  /** \brief Deleted copy constructor. */
  DirectoriesOF(DirectoriesOF const &) = delete;

  /** \brief Construct given Directories information. */
  template <typename OStream> static void info(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "OpenFOAM directories (pass '' for default in [])\n"
           "--------------------------------------------------------------\n"
           "OpenFOAM cases directory [$FOAM_RUN]\n"
           "OpenFOAM case name\n"
           "OpenFOAM case time [last OpenFOAM case time]\n"
           "--------------------------------------------------------------\n";
  }

  /** \brief Output information about current object. */
  template <typename OStream> void info_runtime(OStream &output) const {
    output
        << "--------------------------------------------------------------\n"
           "OpenFOAM directories\n"
           "--------------------------------------------------------------\n"
           "Case directory: " +
               dir_case +
               "\n"
               "Time directory: " +
               time.timeName() + "\n"
        << "--------------------------------------------------------------\n";
  }

  /**
   \param directories Current case directory information.
   \return OpenFOAM runTime instance given Directories information. */
  Foam::Time makeRunTime(Directories const &directories) {
    auto input = useful::open_read(directories.dir_parameters + "/of.dat");
    auto dir =
        useful::read_first_from_line<std::string>(input);
    useful::read_first_from_line(case_name, input);
    useful::read_first_from_line(_time_name_input, input);
    input.close();

    if (useful::is_empty(dir))
      dir = "${FOAM_RUN}";
    useful::expand_env_in_place(useful::expand_home_dir_in_place(dir));

    dir_case = dir + "/" + case_name;

    return {Foam::Time::controlDictName, dir, case_name};
  }
};
} // namespace ptof

#endif /* PTOF_DIRECTORIES_H */
