/**
   \file PTOF/Measurer.h
   \author Tomas Aquino
   \date 07/03/2022
   \brief Objects to measure and output quantities.
*/

#ifndef PTOF_MEASURER_H
#define PTOF_MEASURER_H

#include "CTRW/Meta.h"
#include "CTRW/StateGetter.h"
#include "General/IO.h"
#include "General/Meta.h"
#include "PTOF/Directories.h"
#include "PTOF/Meta.h"
#include "PTOF/Useful.h"
#include <IOobject.H>
#include <algorithm>
#include <array>
#include <fieldTypes.H>
#include <fstream>
#include <iomanip>
#include <string>
#include <type_traits>

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
   \brief Output absorption time, tag, mass, patch (optional), and
   position (optional, optionally accounting for periodicity) of absorbed
   particles.
*/
template <typename Subject, typename Geometry, bool print_patch = false,
          bool print_position = false, bool periodic_position = false>
struct Measurer_absorption_time final : Measurer<Subject, Geometry> {
  Measurer_absorption_time(Subject const &subject, Geometry const &geometry,
                           Directories const &directories,
                           std::string const &identifier, int precision = 8)
      : Measurer<Subject, Geometry>{subject,
                                    geometry,
                                    directories,
                                    std::string{"absorption_time"} +
                                        (print_patch ? "_patch" : "") +
                                        (print_position ? "_position" : "") +
                                        (periodic_position ? "_periodic" : ""),
                                    identifier,
                                    precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())), 0,
            0} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Tag"
            << std::setw(_column_widths[2]) << "Mass";
    if constexpr (print_patch) {
      _patch_names = geometry.mesh().boundaryMesh().names();
      _column_widths[3] =
          1 + std::max(
                  int(std::max_element(
                          _patch_names.begin(), _patch_names.end(),
                          [](Foam::word const &name1, Foam::word const &name2) {
                            return name1.length() < name2.length();
                          })
                          ->length()),
                  int(std::string{"Patch"}.length()));
      _output << std::setw(_column_widths[3]) << "Patch";
    }
    if constexpr (print_position) {
      _column_widths[4] =
          std::max(9 + precision, int(2 + std::string{"Position_"}.length()));
      for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
        _output << std::setw(_column_widths[4])
                << "Position_" + std::to_string(dd);
      }
    }
    if constexpr (meta::has_reinjections_v<typename State::Info>) {
      _column_widths[5] =
          std::max(12, int(1 + std::string{"Reinjections"}.length()));
    }
    _output << "\n";
  }

  Measurer_absorption_time(Subject const &subject, Geometry const &geometry,
                           Directories const &directories,
                           std::string const &identifier,
                           meta::Selector<bool, true> patch,
                           meta::Selector<bool, true> position,
                           int precision = 8)
      : Measurer_absorption_time<Subject, Geometry>{
            subject, geometry, directories, identifier, precision} {}

  Measurer_absorption_time(Subject const &subject, Geometry const &geometry,
                           Directories const &directories,
                           std::string const &identifier,
                           meta::Selector<bool, true> patch,
                           meta::Selector<bool, false> position,
                           int precision = 8)
      : Measurer_absorption_time<Subject, Geometry>{
            subject, geometry, directories, identifier, precision} {}

  Measurer_absorption_time(Subject const &subject, Geometry const &geometry,
                           Directories const &directories,
                           std::string const &identifier,
                           meta::Selector<bool, false> patch,
                           meta::Selector<bool, true> position,
                           int precision = 8)
      : Measurer_absorption_time<Subject, Geometry>{
            subject, geometry, directories, identifier, precision} {}

  Measurer_absorption_time(Subject const &subject, Geometry const &geometry,
                           Directories const &directories,
                           std::string const &identifier,
                           meta::Selector<bool, false> patch,
                           meta::Selector<bool, false> position,
                           int precision = 8)
      : Measurer_absorption_time<Subject, Geometry>{
            subject, geometry, directories, identifier, precision} {}

  void print() override {
    for (auto const &part : _subject.particles()) {
      auto const &state = part.state_new();
      if (state.info.absorbed) {
        _output << std::setw(_column_widths[0]) << state.time
                << std::setw(_column_widths[1]) << state.tag
                << std::setw(_column_widths[2]) << state.mass;
        print_state_patch(state);
        print_state_position(state);
        print_state_reinjections(state);
        _output << "\n";
      }
    }
  }

private:
  using Measurer<Subject, Geometry>::_output;
  using Measurer<Subject, Geometry>::_subject;
  using State = typename Subject::Particle::State;

  Foam::wordList _patch_names;
  std::array<int, 6> _column_widths;

  void print_state_patch(State const &state) {
    if constexpr (print_patch) {
      Foam::label face;
      if constexpr (meta::has_boundary_face_v<typename State::Info>) {
        face = state.info.boundary_face;
      } else {
        face = this->_geometry.mesh_search().findNearestBoundaryFace(
            make_point(state.position));
      }
      _output << std::setw(_column_widths[3])
              << _patch_names[ptof::patch_id(face, this->_geometry.mesh())];
    }
  }

  void print_state_position(State const &state) {
    if constexpr (print_position) {
      if constexpr (periodic_position) {
        using Boundary = std::decay_t<typename Geometry::BoundaryPeriodic>;
        io::print(_output,
                  ctrw::Get_position_periodic<Boundary const &>{
                      this->_geometry.boundary_periodic}(state),
                  _column_widths[4]);
      } else {
        io::print(_output, state.position, _column_widths[4]);
      }
    }
  }

  void print_state_reinjections(State const &state) {
    if constexpr (meta::has_reinjections_v<typename Subject::State::Info>) {
      _output << std::setw(_column_widths[5]) << state.info.reinjections;
    }
  }
  };
template <typename Subject, typename Geometry>
Measurer_absorption_time(Subject const &, Geometry const &, Directories const &,
                         std::string const &, meta::Selector<bool, true>,
                         meta::Selector<bool, true>, int)
    -> Measurer_absorption_time<Subject, Geometry, true, true>;
template <typename Subject, typename Geometry>
Measurer_absorption_time(Subject const &, Geometry const &, Directories const &,
                         std::string const &, meta::Selector<bool, true>,
                         meta::Selector<bool, true>)
    -> Measurer_absorption_time<Subject, Geometry, true, true>;
template <typename Subject, typename Geometry>
Measurer_absorption_time(Subject const &, Geometry const &, Directories const &,
                         std::string const &, meta::Selector<bool, true>,
                         meta::Selector<bool, false>, int)
    -> Measurer_absorption_time<Subject, Geometry, true, false>;
template <typename Subject, typename Geometry>
Measurer_absorption_time(Subject const &, Geometry const &, Directories const &,
                         std::string const &, meta::Selector<bool, true>,
                         meta::Selector<bool, false>)
    -> Measurer_absorption_time<Subject, Geometry, true, false>;
template <typename Subject, typename Geometry>
Measurer_absorption_time(Subject const &, Geometry const &, Directories const &,
                         std::string const &, meta::Selector<bool, false>,
                         meta::Selector<bool, true>, int)
    -> Measurer_absorption_time<Subject, Geometry, false, true>;
template <typename Subject, typename Geometry>
Measurer_absorption_time(Subject const &, Geometry const &, Directories const &,
                         std::string const &, meta::Selector<bool, false>,
                         meta::Selector<bool, true>)
    -> Measurer_absorption_time<Subject, Geometry, false, true>;
template <typename Subject, typename Geometry>
Measurer_absorption_time(Subject const &, Geometry const &, Directories const &,
                         std::string const &, meta::Selector<bool, false>,
                         meta::Selector<bool, false>, int)
    -> Measurer_absorption_time<Subject, Geometry, false, false>;
template <typename Subject, typename Geometry>
Measurer_absorption_time(Subject const &, Geometry const &, Directories const &,
                         std::string const &, meta::Selector<bool, false>,
                         meta::Selector<bool, false>)
    -> Measurer_absorption_time<Subject, Geometry, false, false>;
} // namespace ptof

#endif /* PTOF_MEASURER_H */
