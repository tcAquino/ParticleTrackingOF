/**
   \file PTOF/Measurer.h
   \author Tomás Aquino
   \date 07/03/2022
   \brief Objects to measure and output quantities.
*/

#ifndef PTOF_MEASURER_H
#define PTOF_MEASURER_H

#include "General/IO.h"
#include "PTOF/Directories.h"
#include "PTOF/Meta.h"
#include <IOobject.H>
#include <algorithm>
#include <array>
#include <fieldTypes.H>
#include <fstream>
#include <iomanip>
#include <string>

namespace ptof {
/**
   \class Measurer PTOF/Measurer.h "PTOF/Measurer.h"
   \brief Abstract polymorphic output handler for outputting some quantity.
*/
template <typename Subject, typename Geometry> struct Measurer {
  virtual ~Measurer() = default;

  /** \brief Make measurement and output. */
  virtual void print() = 0;

protected:
  /**
     \brief Constructor.
     \param subject CTRW object to measure.
     \param geometry Domain geometry info and utilities.
     \param directories Current case directory information.
     \param output_name Name of output type.
     \param identifier String to include in names of output files.
     \param precision Number of digits after decimal point in output, in
     scientific notation.
  */
  Measurer(Subject const &subject, Geometry const &geometry,
           Directories const &directories, std::string const &output_name,
           std::string const &identifier, int precision = 8)
      : _subject{subject}, _geometry{geometry}, _locator{geometry.locator},
        _output{open_write(directories, output_name, identifier)} {
    _output << std::setprecision(precision) << std::scientific;
  }

  Measurer(Subject &&, Geometry const &, Directories const &,
           std::string const &, std::string const &, int = 8) = delete;

  Measurer(Subject const &, Geometry &&, Directories const &,
           std::string const &, std::string const &, int = 8) = delete;

  Measurer(Subject &&, Geometry &&, Directories const &, std::string const &,
           std::string const &, int = 8) = delete;

  Subject const &_subject;   /**< CTRW object to measure. */
  Geometry const &_geometry; /**< Domain geometry info and utilities. */
  typename Geometry::Locator const
      &_locator;         /**< Object to locate positions in mesh. */
  std::ofstream _output; /**< Output stream. */

  /**
     \param directories Current case directory information.
     \param output_name Name of output type.
     \param identifier String to include in names of output files.
     \return Output stream.
  */
  std::ofstream open_write(Directories const &directories,
                           std::string const &output_name,
                           std::string const &identifier) {
    return io::open_write(directories.dir_output + "/Data_" + output_name +
                          (identifier.empty() ? "" : "_" + identifier) +
                          ".dat");
  }
};

/**
   \class Measurer_absorption_time PTOF/Measurer.h "PTOF/Measurer.h"
   \brief  Output absorption times, tags, and masses of absorbed particles.
*/
template <typename Subject, typename Geometry>
struct Measurer_absorption_time final : Measurer<Subject, Geometry> {
  Measurer_absorption_time(Subject const &subject, Geometry const &geometry,
                           Directories const &directories,
                           std::string const &identifier, int precision = 8)
      : Measurer<Subject, Geometry>{subject,           geometry,   directories,
                                    "absorption_time", identifier, precision},
        _patch_names{geometry.mesh().boundaryMesh().names()},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())),
            std::max(_patch_names.empty()
                         ? 0
                         : int(1 + std::max_element(
                                       _patch_names.begin(), _patch_names.end(),
                                       [](Foam::word const &name1,
                                          Foam::word const &name2) {
                                         return name1.length() < name2.length();
                                       })
                                       ->length()),
                     int(1 + std::string{"Patch"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Tag"
            << std::setw(_column_widths[2]) << "Mass";
    if constexpr (meta::has_patch_id_v<
                      typename Subject::Particle::State::Info>) {
      _output << std::setw(_column_widths[3]) << "Patch";
      _patch_names = geometry.mesh().boundaryMesh().names();
    }
    _output << "\n";
  }

  void print() override {
    for (auto const &part : _subject.particles()) {
      auto const &state = part.state_new();
      if (state.info.absorbed) {
        _output << std::setw(_column_widths[0]) << state.time
                << std::setw(_column_widths[1]) << state.tag
                << std::setw(_column_widths[2]) << state.mass;
        if constexpr (meta::has_patch_id_v<
                          typename Subject::Particle::State::Info>)
          _output << std::setw(_column_widths[3])
                  << _patch_names[state.info.patch_id];
        _output << "\n";
      }
    }
  }

private:
  using Measurer<Subject, Geometry>::_output;
  using Measurer<Subject, Geometry>::_subject;
  Foam::List<Foam::word> _patch_names;
  std::array<int, 4> _column_widths;
};
} // namespace ptof

#endif /* PTOF_MEASURER_H */
