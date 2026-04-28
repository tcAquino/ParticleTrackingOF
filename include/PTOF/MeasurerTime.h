/**
 * @file   MeasurerTime.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Mon Mar  7 00:00:00 2022
 *
 * @brief  Objects to measure and output quantities given a measurement time.
 */

#ifndef PTOF_MEASURERTIME_H
#define PTOF_MEASURERTIME_H

#include "CTRW/StateGetter.h"
#include "General/IO.h"
#include "General/Meta.h"
#include "PTOF/BoundaryInfo.h"
#include "PTOF/CheckOptions.h"
#include "PTOF/Directories.h"
#include "PTOF/Field.h"
#include "PTOF/Useful.h"
#include <IOobject.H>
#include <algorithm>
#include <array>
#include <fieldTypes.H>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#include <volFields.H>

namespace ptof {
/**
 * @brief Abstract polymorphic output handler for outputting time and some
 *        quantity.
 */
template <typename Subject, typename Geometry> class MeasurerTime {
public:
  /** @brief Destructor. */
  virtual ~MeasurerTime() = default;

  /**
   * @brief Make measurement and output or store, given current time \c time.
   */
  virtual void operator()(double time) = 0;

  /** @brief Output stored information (do nothing if not overriden). */
  virtual void print() {};

  /** @brief Update internal state. */
  virtual void update(double time, Foam::instant const &previous_time_of_change,
                      Foam::instant const &next_time_of_change) {}

protected:
  /**
   * @brief Set up output streams for requested output types.
   *
   * @param subject CTRW object to measure.
   *
   * @param geometry Domain geometry info and utilities.
   *
   * @param directories Current case directory information.
   *
   * @param output_name Name of output type.
   *
   * @param identifier String to include in names of output files.
   *
   * @param precision Number of digits after decimal point in output, in
   *                  scientific notation.
   */
  MeasurerTime(Subject const &subject, Geometry const &geometry,
               Directories const &directories, std::string output_name,
               std::string const &identifier, int precision = 8)
      : _subject{subject}, _geometry{geometry}, _locator{geometry.locator},
        _output_name{output_name},
        _output{open_write(directories, _output_name, identifier)} {
    _output << std::setprecision(precision) << std::scientific;
  }

  MeasurerTime(Subject &&, Geometry const &, Directories const &,
               std::string const &_name, std::string const &, int = 8) = delete;

  MeasurerTime(Subject const &, Geometry &&, Directories const &,
               std::string const &_name, std::string const &, int = 8) = delete;

  MeasurerTime(Subject &&, Geometry &&, Directories const &,
               std::string const &_name, std::string const &, int = 8) = delete;

  Subject const &_subject;   /**< CTRW object to measure. */
  Geometry const &_geometry; /**< Domain geometry info and utilities. */
  typename Geometry::Locator const
      &_locator; /**< Object to locate positions in mesh. */
  std::string _output_name;
  std::ofstream _output; /**< Output stream. */

  /**
   * @brief Open output file for a given output type.
   *
   * @param directories Current case directory information.
   *
   * @param output_name Name of output type.
   *
   * @param identifier String to include in names of output files.
   *
   * @return Output stream.
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
 * @brief Particle tags, masses, cells, and positions, as well as times and
 *        number of particles output per time to a separate file (with \c _times
 *        added to output name).
 *
 * @tparam periodic If \c true, give position accounting for periodicity
 *                  (possibly outside domain).
 */
template <typename Subject, typename Geometry, bool periodic = false>
class MeasurerTime_position final : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_position(Subject const &subject, Geometry const &geometry,
                        Directories const &directories,
                        std::string const &identifier, int precision = 8)
      : MeasurerTime<
            Subject,
            Geometry>{subject,
                      geometry,
                      directories,
                      std::string{"position"} + (periodic ? "_periodic" : ""),
                      identifier,
                      precision},
        _column_widths{
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())),
            std::max(12 + precision, int(1 + std::string{"Cell"}.length())),
            std::max(9 + precision,
                     int(2 + std::string{"Position_"}.length()))},
        _output_times{this->open_write(
            directories, this->_output_name + "_times", identifier)},
        _column_widths_times{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Nr_particles"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Tag"
            << std::setw(_column_widths[1]) << "Mass"
            << std::setw(_column_widths[2]) << "Cell";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
      _output << std::setw(_column_widths[3])
              << "Position_" + std::to_string(dd);
    }
    _output << "\n";

    _output_times << std::setw(_column_widths_times[0]) << "Time"
                  << std::setw(_column_widths_times[1]) << "Nr_particles"
                  << "\n";
  }

  template <bool periodic_v>
  MeasurerTime_position(Subject const &subject, Geometry const &geometry,
                        Directories const &directories,
                        std::string const &identifier,
                        meta::Selector<bool, periodic_v>, int precision = 8)
      : MeasurerTime_position<Subject, Geometry, periodic_v>{
            subject, geometry, directories, identifier, precision} {}

  void operator()(double time) override {
    std::size_t nr_particles = 0;
    for (auto const &part : _subject) {
      auto const &state_old = part.state_old();
      if (adsorbed(state_old) || !brackets_time(part, time)) {
        continue;
      }
      ++nr_particles;

      auto const &state_new = part.state_new();
      auto [position, cell] = position_interpolated_with_cell<true, periodic>(
          part, time, this->_geometry);
      _output << std::setw(_column_widths[0]) << state_old.tag;
      _output << std::setw(_column_widths[1])
              << ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new, state_old);
      _output << std::setw(_column_widths[2]) << cell;
      io::print(_output, position, _column_widths[3]) << "\n";
    }

    _output_times << std::setw(_column_widths_times[0]) << time
                  << std::setw(_column_widths_times[1]) << nr_particles << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 4> _column_widths;
  std::ofstream _output_times;
  std::array<int, 2> _column_widths_times;
};

/**
 * @brief Output tags, masses, cells, and positions within regions specified by
 *        masks, as well as times and number of particles output per time to a
 *        separate file (with \c _times added to output name).
 *
 * @tparam periodic If \c true, give position accounting for periodicity
 *                  (possibly outside domain).
 */
template <typename Subject, typename Geometry, typename Mask,
          bool periodic = false>
class MeasurerTime_position_in_regions final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_position_in_regions(
      Subject const &subject, Geometry const &geometry,
      Directories const &directories, std::string const &identifier,
      std::vector<std::reference_wrapper<const Mask>> masks,
      std::vector<double> thresholds = {}, int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,
                                        geometry,
                                        directories,
                                        std::string{"position_in_regions"} +
                                            (periodic ? "_period"
                                                        "ic"
                                                      : ""),
                                        identifier,
                                        precision},
        _masks{masks}, _thresholds{thresholds},
        _column_widths{
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())),
            std::max(12, int(1 + std::string{"Cell"}.length())),
            std::max(9 + precision, int(2 + std::string{"Position_"}.length())),
            std::max(9, int(4 + std::string{"In_region_"}.length()))},
        _output_times{this->open_write(
            directories, this->_output_name + "_times", identifier)},
        _column_widths_times{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Nr_particles"}.length()))} {
    _thresholds.resize(masks.size(), 0.);
    _output << std::setw(_column_widths[0]) << "Tag";
    _output << std::setw(_column_widths[1]) << "Mass";
    _output << std::setw(_column_widths[2]) << "Cell";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
      _output << std::setw(_column_widths[3])
              << "Position_" + std::to_string(dd);
    }
    if (_masks.size() > 1) {
      for (std::size_t ii = 0; ii < _masks.size(); ++ii) {
        _output << std::setw(_column_widths[4])
                << "In_region_" + std::to_string(ii);
      }
    }
    _output << "\n";

    _output_times << std::setw(_column_widths_times[0]) << "Time"
                  << std::setw(_column_widths_times[1]) << "Nr_particles"
                  << "\n";
  }

  template <bool periodic_v>
  MeasurerTime_position_in_regions(
      Subject const &subject, Geometry const &geometry,
      Directories const &directories, std::string const &identifier,
      std::vector<std::reference_wrapper<const Mask>> masks,
      meta::Selector<bool, periodic_v>, std::vector<double> thresholds = {},
      int precision = 8)
      : MeasurerTime_position_in_regions<Subject, Geometry, Mask, periodic_v>{
            subject, geometry,   directories, identifier,
            masks,   thresholds, precision} {}

  void operator()(double time) override {
    std::size_t nr_particles = 0;
    for (auto const &part : _subject) {
      auto const &state_old = part.state_old();
      if (adsorbed(state_old) || !brackets_time(part, time)) {
        continue;
      }
      auto const &state_new = part.state_new();
      if (outside(state_new.cell) || outside(state_old.cell)) {
        continue;
      }

      std::vector<int> in_region(_masks.size(), 0);
      for (std::size_t ii = 0; ii < _masks.size(); ++ii) {
        if (evaluate(_masks[ii].get(), state_new) >= _thresholds[ii] &&
            evaluate(_masks[ii].get(), state_old) >= _thresholds[ii]) {
          in_region[ii] = 1;
        }
      }
      if (std::any_of(in_region.begin(), in_region.end(),
                      [](int ii) { return ii > 0; })) {
        ++nr_particles;
        auto [position, cell] = position_interpolated_with_cell<true, periodic>(
            part, time, this->_geometry);
        _output << std::setw(_column_widths[0]) << state_old.tag;
        _output << std::setw(_column_widths[1])
                << ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new,
                                                            state_old);
        _output << std::setw(_column_widths[2]) << cell;
        io::print(_output, position, _column_widths[3]);
        if (_masks.size() > 1) {
          io::print(_output, in_region, _column_widths[4]);
        }
        _output << "\n";
      }
    }

    _output_times << std::setw(_column_widths_times[0]) << time
                  << std::setw(_column_widths_times[1]) << nr_particles << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::vector<std::reference_wrapper<const Mask>> _masks;
  std::vector<double> _thresholds;
  std::array<int, 5> _column_widths;
  std::ofstream _output_times;
  std::array<int, 2> _column_widths_times;
};

/**
 * @brief Output time and mean position (weighted by mass).
 *
 * @tparam periodic If \c true, give position accounting for periodicity
 *                  (possibly outside domain).
 */
template <typename Subject, typename Geometry, bool periodic = false>
class MeasurerTime_position_mean final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_position_mean(Subject const &subject, Geometry const &geometry,
                             Directories const &directories,
                             std::string const &identifier, int precision = 8)
      : MeasurerTime<
            Subject,
            Geometry>{subject,
                      geometry,
                      directories,
                      std::string{"position_mean"} + (periodic ? "_periodic" : ""),
                      identifier,
                      precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(2 + std::string{"Position_mean_"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
      _output << std::setw(_column_widths[1])
              << "Position_mean_" + std::to_string(dd);
    }
    _output << "\n";
  }

  template <bool periodic_v>
  MeasurerTime_position_mean(Subject const &subject, Geometry const &geometry,
                             Directories const &directories,
                             std::string const &identifier,
                             meta::Selector<bool, periodic_v>,
                             int precision = 8)
      : MeasurerTime_position_mean<Subject, Geometry, periodic_v>{
            subject, geometry, directories, identifier, precision} {}

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    io::print(_output,
              position_mean<false, periodic>(_subject, time, this->_geometry),
              _column_widths[1]);
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 2> _column_widths;
};

/**
 * @brief Output time and mean of absolute value of position (weighted by mass).
 *
 * @tparam periodic If \c true, give position accounting for periodicity
 *                  (possibly outside domain).
 */
template <typename Subject, typename Geometry, bool periodic = false>
class MeasurerTime_position_abs_mean final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_position_abs_mean(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 int precision = 8)
      : MeasurerTime<
            Subject,
            Geometry>{subject,
                      geometry,
                      directories,
                      std::string{"position_abs_mean"} + (periodic ? "_periodic" : ""),
                      identifier,
                      precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Position_abs_mean"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    _output << std::setw(_column_widths[1]) << "Position_abs_mean";
    _output << "\n";
  }

  template <bool periodic_v>
  MeasurerTime_position_abs_mean(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 meta::Selector<bool, periodic_v>,
                                 int precision = 8)
      : MeasurerTime_position_abs_mean<Subject, Geometry, periodic_v>{
            subject, geometry, directories, identifier, precision} {}

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    io::print(_output,
              position_mean<false, periodic>(
                  _subject, time, this->_geometry,
                  [](auto const &position) { return op::abs(position); }),
              _column_widths[1]);
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 2> _column_widths;
};

/**
 * @brief Output time and second moment of position (weighted by mass).
 *
 * @tparam periodic If \c true, give position accounting for periodicity
 *                  (possibly outside domain).
 */
template <typename Subject, typename Geometry, bool periodic = false>
class MeasurerTime_position_second_moment final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_position_second_moment(Subject const &subject,
                                      Geometry const &geometry,
                                      Directories const &directories,
                                      std::string const &identifier,
                                      int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,
                                        geometry,
                                        directories,
                                        std::string{"position_second_moment"} +
                                            (periodic ? "_periodic" : ""),
                                        identifier,
                                        precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(
                9 + precision,
                int(2 + std::string{"Position_second_moment_"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
      _output << std::setw(_column_widths[1])
              << "Position_second_moment_" + std::to_string(dd);
    }
    _output << "\n";
  }

  template <bool periodic_v>
  MeasurerTime_position_second_moment(Subject const &subject,
                                      Geometry const &geometry,
                                      Directories const &directories,
                                      std::string const &identifier,
                                      meta::Selector<bool, periodic_v>,
                                      int precision = 8)
      : MeasurerTime_position_second_moment<Subject, Geometry, periodic_v>{
            subject, geometry, directories, identifier, precision} {}

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    io::print(_output,
              position_second_moment<false, periodic>(_subject, time,
                                                      this->_geometry),
              _column_widths[1]);
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 2> _column_widths;
};

/**
 * @brief Output time and nth moment of position (weighted by mass).
 *
 * @tparam periodic If \c true, give position accounting for periodicity
 *                  (possibly outside domain).
 */
template <typename Subject, typename Geometry, bool periodic = false>
class MeasurerTime_position_nth_moment final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_position_nth_moment(Subject const &subject, double nn,
                                   Geometry const &geometry,
                                   Directories const &directories,
                                   std::string const &identifier,
                                   int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,       geometry,   directories,
                                        make_name(nn), identifier, precision},
        _nn{nn},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision,
                     int(2 + std::string{"Position_moment_"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
      _output << std::setw(_column_widths[1])
              << "Position_moment_" + std::to_string(dd);
    }
    _output << "\n";
  }

  template <bool periodic_v>
  MeasurerTime_position_nth_moment(Subject const &subject, double nn,
                                   Geometry const &geometry,
                                   Directories const &directories,
                                   std::string const &identifier,
                                   meta::Selector<bool, periodic_v>,
                                   int precision = 8)
      : MeasurerTime_position_nth_moment<Subject, Geometry, periodic_v>{
            subject, geometry, directories, identifier, precision} {}

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    io::print(
        _output,
        position_moment<false, periodic>(_subject, _nn, time, this->_geometry),
        _column_widths[1]);
    _output << "\n";
  }

private:
  static std::string make_name(double nn) {
    std::stringstream stream;
    stream << std::scientific << std::setprecision(2);
    stream << "_" << nn;
    return std::string{"position_moment_"} + stream.str() +
           (periodic ? "_periodic" : "");
  }

  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  double _nn;
  std::array<int, 2> _column_widths;
};

/**
 * @brief Output time and moment of position with specified exponents for each
 *        coordinate (weighted by mass).
 *
 * @tparam periodic If \c true, give position accounting for periodicity
 *                  (possibly outside domain).
 */
template <typename Subject, typename Geometry, bool periodic = false>
class MeasurerTime_position_moment final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_position_moment(Subject const &subject,
                               std::vector<double> exponents,
                               Geometry const &geometry,
                               Directories const &directories,
                               std::string const &identifier, int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,     geometry,
                                        directories, make_name(exponents),
                                        identifier,  precision},
        _exponents{exponents},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision,
                     int(2 + std::string{"Position_moment_"}.length()))} {
    if (_exponents.size() != Geometry::dim) {
      throw std::runtime_error{
          "Number of moment exponents is not the same as spatial dimension"};
    }
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
      _output << std::setw(_column_widths[1])
              << "Position_moment_" + std::to_string(dd);
    }
    _output << "\n";
  }

  template <bool periodic_v>
  MeasurerTime_position_moment(Subject const &subject,
                               std::vector<double> exponents,
                               Geometry const &geometry,
                               Directories const &directories,
                               std::string const &identifier,
                               meta::Selector<bool, periodic_v>,
                               int precision = 8)
      : MeasurerTime_position_moment<Subject, Geometry, periodic_v>{
            subject, geometry, directories, identifier, precision} {}

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    io::print(_output,
              position_moment<false, periodic>(_subject, _exponents, time,
                                               this->_geometry),
              _column_widths[1]);
    _output << "\n";
  }

private:
  static std::string make_name(std::vector<double> const &exponents) {
    std::stringstream stream;
    stream << std::scientific << std::setprecision(2);
    for (auto exp : exponents) {
      stream << "_" << exp;
    }
    return std::string{"position_moment"} + stream.str() +
           (periodic ? "_periodic" : "");
  }

  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::vector<double> _exponents;
  std::array<int, 2> _column_widths;
};

/**
 * @brief Output time and position variance (weighted by mass).
 *
 * @tparam periodic If \c true, give position accounting for periodicity
 *                  (possibly outside domain).
 */
template <typename Subject, typename Geometry, bool periodic = false>
class MeasurerTime_position_variance final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_position_variance(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 int precision = 8)
      : MeasurerTime<
            Subject,
            Geometry>{subject,
                      geometry,
                      directories,
                      std::string{"position_variance"} + (periodic ? "_periodic" : ""),
                      identifier,
                      precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision,
                     int(2 + std::string{"Position_variance_"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
      _output << std::setw(_column_widths[1])
              << "Position_variance_" + std::to_string(dd);
    }
    _output << "\n";
  }

  template <bool periodic_v>
  MeasurerTime_position_variance(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 meta::Selector<bool, periodic_v>,
                                 int precision = 8)
      : MeasurerTime_position_variance<Subject, Geometry, periodic_v>{
            subject,
            geometry,
            directories,
            std::string{"position_variance"} + (periodic ? "_periodic" : ""),
            identifier,
            precision} {}

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    io::print(
        _output,
        position_variance<false, periodic>(_subject, time, this->_geometry),
        _column_widths[1]);
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 2> _column_widths;
};

/** @brief Output time and total mass. */
template <typename Subject, typename Geometry>
class MeasurerTime_mass final : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_mass(Subject const &subject, Geometry const &geometry,
                    Directories const &directories,
                    std::string const &identifier, int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject, geometry,   directories,
                                        "mass",  identifier, precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Mass" << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time
            << std::setw(_column_widths[1]) << mass(_subject, time) << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 2> _column_widths;
};

/** @brief Output time and total absorbed mass. */
template <typename Subject, typename Geometry>
class MeasurerTime_mass_absorbed final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_mass_absorbed(Subject const &subject, Geometry const &geometry,
                             Directories const &directories,
                             std::string const &identifier, int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,     geometry,
                                        directories, "mass_absorbed",
                                        identifier,  precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Mass" << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time
            << std::setw(_column_widths[1]) << mass_absorbed(_subject, time)
            << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 2> _column_widths;
};

/** @brief Output time and total adsorbed mass. */
template <typename Subject, typename Geometry>
class MeasurerTime_mass_adsorbed final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_mass_adsorbed(Subject const &subject, Geometry const &geometry,
                             Directories const &directories,
                             std::string const &identifier, int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,     geometry,
                                        directories, "mass_adsorbed",
                                        identifier,  precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Mass" << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time
            << std::setw(_column_widths[1]) << mass_adsorbed(_subject, time)
            << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 2> _column_widths;
};

/** @brief Output time and total mass in regions specified by maks. */
template <typename Subject, typename Geometry, typename Mask>
class MeasurerTime_mass_in_regions final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_mass_in_regions(
      Subject const &subject, Geometry const &geometry,
      Directories const &directories, std::string const &identifier,
      std::vector<std::reference_wrapper<const Mask>> masks,
      std::vector<double> thresholds = {}, int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,     geometry,
                                        directories, "mass_in_regions",
                                        identifier,  precision},
        _masks{masks}, _thresholds{thresholds},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision,
                     int(4 + std::string{"Mass_region_"}.length()))} {
    _thresholds.resize(masks.size(), 0.);
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t ii = 0; ii < masks.size(); ++ii) {
      _output << std::setw(_column_widths[1])
              << "Mass_region_" + std::to_string(ii);
    }
    _output << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    io::print(_output, mass(_subject, time, _masks, _thresholds),
              _column_widths[1]);
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using MeasurerTime<Subject, Geometry>::_locator;
  std::vector<std::reference_wrapper<const Mask>> _masks;
  std::vector<double> _thresholds;
  std::array<int, 2> _column_widths;
};

/**
 * @brief Output tags, and scalar field values, as well as times and number of
 *        particles output per time to a separate file (with \c _times added to
 *        output name).
 */
template <typename Subject, typename Geometry,
          typename Field = ScalarField_Interpolation<
              Foam::volScalarField, typename Geometry::Locator const &,
              InterpolationTypes::Linear, CheckOptions::Check>>
class MeasurerTime_scalar_field final : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_scalar_field(Subject const &subject, Field &&field,
                            Geometry const &geometry,
                            Directories const &directories,
                            std::string const &identifier,
                            std::string const &field_name, int precision = 8)
      : MeasurerTime<Subject,
                     Geometry>{subject,
                               geometry,
                               directories,
                               std::string{"scalar_field_"} + field_name,
                               identifier,
                               precision},
        _field_name{field_name}, _field{std::forward<Field>(field)},
        _column_widths{
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())),
            std::max(9 + precision, int(1 + field_name.length()))},
        _output_times{this->open_write(
            directories, this->_output_name + "_times", identifier)},
        _column_widths_times{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Nr_particles"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Tag"
            << std::setw(_column_widths[1]) << "Mass"
            << std::setw(_column_widths[2]) << field_name << "\n";

    _output_times << std::setw(_column_widths_times[0]) << "Time"
                  << std::setw(_column_widths_times[1]) << "Nr_particles"
                  << "\n";
  }

  MeasurerTime_scalar_field(Subject const &subject, Geometry const &geometry,
                            Directories const &directories,
                            std::string const &identifier,
                            std::string const &field_name, int precision = 8)
      : MeasurerTime_scalar_field{
            subject,
            ScalarField_Interpolation<
                Foam::volScalarField, typename Geometry::Locator const &,
                InterpolationTypes::Linear, CheckOptions::Check>{
                {Foam::IOobject{field_name, geometry.mesh().time().timeName(),
                                geometry.mesh(), Foam::IOobject::MUST_READ,
                                Foam::IOobject::NO_WRITE},
                 geometry.mesh()},
                geometry.locator},
            geometry,
            directories,
            identifier,
            field_name,
            precision} {}

  void operator()(double time) override {
    std::size_t nr_particles = 0;
    for (auto const &part : _subject) {
      auto const &state_old = part.state_old();
      if (adsorbed(state_old) || !brackets_time(part, time)) {
        continue;
      }
      ++nr_particles;

      auto const &state_new = part.state_new();
      _output << std::setw(_column_widths[0]) << state_old.tag
              << std::setw(_column_widths[1])
              << ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new, state_old)
              << std::setw(_column_widths[2])
              << ctrw::Get_interp{time, ctrw::Get_property{_field}}(state_new,
                                                                    state_old)
              << "\n";
    }

    _output_times << std::setw(_column_widths_times[0]) << time
                  << std::setw(_column_widths_times[1]) << nr_particles << "\n";
  }

  void update(double time = 0.,
              Foam::instant const &previous_time_of_change = {},
              Foam::instant const &next_time_of_change = {}) override {
    if constexpr (!std::is_reference_v<Field>) {
      _field.set(
          {Foam::IOobject{_field_name, _geometry.mesh().time().timeName(),
                          _geometry.mesh(), Foam::IOobject::MUST_READ,
                          Foam::IOobject::NO_WRITE},
           _geometry.mesh()});
    }
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using MeasurerTime<Subject, Geometry>::_geometry;
  std::string _field_name;
  Field _field;
  std::array<int, 3> _column_widths;
  std::ofstream _output_times;
  std::array<int, 2> _column_widths_times;
};
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_scalar_field(Subject const &, Field &&, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &, int)
    -> MeasurerTime_scalar_field<Subject, Geometry, Field>;
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_scalar_field(Subject const &, Field &&, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &)
    -> MeasurerTime_scalar_field<Subject, Geometry, Field>;
template <typename Subject, typename Geometry>
MeasurerTime_scalar_field(Subject const &, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &, int)
    -> MeasurerTime_scalar_field<
        Subject, Geometry,
        ScalarField_Interpolation<
            Foam::volScalarField, typename Geometry::Locator const &,
            InterpolationTypes::Linear, CheckOptions::Check>>;
template <typename Subject, typename Geometry>
MeasurerTime_scalar_field(Subject const &, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &)
    -> MeasurerTime_scalar_field<
        Subject, Geometry,
        ScalarField_Interpolation<
            Foam::volScalarField, typename Geometry::Locator const &,
            InterpolationTypes::Linear, CheckOptions::Check>>;

/**
 * @brief Output tags, and vector field values, as well as times and number of
 *        particles output per time to a separate file (with \c _times added to
 *        output name).
 */
template <typename Subject, typename Geometry,
          typename Field = VectorField_Interpolation<
              Foam::volVectorField, typename Geometry::Locator const &,
              InterpolationTypes::Linear, CheckOptions::Check>>
class MeasurerTime_vector_field final : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_vector_field(Subject const &subject, Field &&field,
                            Geometry const &geometry,
                            Directories const &directories,
                            std::string const &identifier,
                            std::string const &field_name, int precision = 8)
      : MeasurerTime<Subject,
                     Geometry>{subject,
                               geometry,
                               directories,
                               std::string{"vector_field_"} + field_name,
                               identifier,
                               precision},
        _field{std::forward<Field>(field)},
        _column_widths{
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())),
            std::max(9 + precision, int(2 + (field_name + "_").length()))},
        _output_times{this->open_write(
            directories, this->_output_name + "_times", identifier)},
        _column_widths_times{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Nr_particles"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Tag"
            << std::setw(_column_widths[1]) << "Mass";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
      _output << std::setw(_column_widths[2])
              << field_name + "_" + std::to_string(dd);
    }
    _output << "\n";

    _output_times << std::setw(_column_widths_times[0]) << "Time"
                  << std::setw(_column_widths_times[1]) << "Nr_particles"
                  << "\n";
  }

  MeasurerTime_vector_field(Subject const &subject, Geometry const &geometry,
                            Directories const &directories,
                            std::string const &identifier,
                            std::string const &field_name, int precision = 8)
      : MeasurerTime_vector_field{
            subject,
            VectorField_Interpolation<
                Foam::volVectorField, typename Geometry::Locator const &,
                InterpolationTypes::Linear, CheckOptions::Check>{
                {Foam::IOobject{field_name, geometry.mesh().time().timeName(),
                                geometry.mesh(), Foam::IOobject::MUST_READ,
                                Foam::IOobject::NO_WRITE},
                 geometry.mesh()},
                geometry.locator},
            geometry,
            directories,
            identifier,
            field_name,
            precision} {}

  void operator()(double time) override {
    std::size_t nr_particles = 0;
    for (auto const &part : _subject) {
      auto const &state_old = part.state_old();
      if (adsorbed(state_old) || !brackets_time(part, time)) {
        continue;
      }
      ++nr_particles;

      auto const &state_new = part.state_new();
      _output << std::setw(_column_widths[0]) << state_old.tag
              << std::setw(_column_widths[1])
              << ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new, state_old);
      io::print(_output,
                ctrw::Get_interp{time, ctrw::Get_property{_field}}(state_new,
                                                                   state_old),
                _column_widths[2])
          << "\n";
    }

    _output_times << std::setw(_column_widths_times[0]) << time
                  << std::setw(_column_widths_times[1]) << nr_particles << "\n";
  }

  void update(double time = 0.,
              Foam::instant const &previous_time_of_change = {},
              Foam::instant const &next_time_of_change = {}) override {
    if constexpr (!std::is_reference_v<Field>) {
      _field.set(
          {Foam::IOobject{_field_name, _geometry.mesh().time().timeName(),
                          _geometry.mesh(), Foam::IOobject::MUST_READ,
                          Foam::IOobject::NO_WRITE},
           _geometry.mesh()});
    }
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using MeasurerTime<Subject, Geometry>::_geometry;
  std::string _field_name;
  Field _field;
  std::array<int, 3> _column_widths;
  std::ofstream _output_times;
  std::array<int, 2> _column_widths_times;
};
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_vector_field(Subject const &, Field &&, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &, int)
    -> MeasurerTime_vector_field<Subject, Geometry, Field>;
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_vector_field(Subject const &, Field &&, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &)
    -> MeasurerTime_vector_field<Subject, Geometry, Field>;
template <typename Subject, typename Geometry>
MeasurerTime_vector_field(Subject const &, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &, int)
    -> MeasurerTime_vector_field<
        Subject, Geometry,
        VectorField_Interpolation<
            Foam::volVectorField, typename Geometry::Locator const &,
            InterpolationTypes::Linear, CheckOptions::Check>>;
template <typename Subject, typename Geometry>
MeasurerTime_vector_field(Subject const &, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &)
    -> MeasurerTime_vector_field<
        Subject, Geometry,
        VectorField_Interpolation<
            Foam::volVectorField, typename Geometry::Locator const &,
            InterpolationTypes::Linear, CheckOptions::Check>>;

/**
 * @brief Output tags, and tensor field values, as well as times and number of
 *        particles output per time to a separate file (with \c _times added to
 *        output name).
 */
template <typename Subject, typename Geometry,
          typename Field = Foam::volTensorField>
class MeasurerTime_tensor_field final : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_tensor_field(Subject const &subject, Field &&field,
                            Geometry const &geometry,
                            Directories const &directories,
                            std::string const &identifier,
                            std::string const &field_name, int precision = 8)
      : MeasurerTime<Subject,
                     Geometry>{subject,
                               geometry,
                               directories,
                               std::string{"tensor_field_"} + field_name,
                               identifier,
                               precision},
        _field{std::forward<Field>(field)},
        _column_widths{
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())),
            std::max(9 + precision, int(3 + (field_name + "_").length()))},
        _output_times{this->open_write(
            directories, this->_output_name + "_times", identifier)},
        _column_widths_times{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Nr_particles"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Tag"
            << std::setw(_column_widths[1]) << "Mass";
    for (std::size_t dd1 = 0; dd1 < Geometry::dim; ++dd1) {
      for (std::size_t dd2 = 0; dd2 < Geometry::dim; ++dd2) {
        _output << std::setw(_column_widths[2])
                << field_name + "_" + std::to_string(dd1) + std::to_string(dd2);
      }
    }

    _output_times << std::setw(_column_widths_times[0]) << "Time"
                  << std::setw(_column_widths_times[1]) << "Nr_particles"
                  << "\n";
  }

  MeasurerTime_tensor_field(Subject const &subject, Geometry const &geometry,
                            Directories const &directories,
                            std::string const &identifier,
                            std::string const &field_name, int precision = 8)
      : MeasurerTime_tensor_field{
            subject,
            Foam::volTensorField{
                Foam::IOobject{field_name, geometry.mesh().time().timeName(),
                               geometry.mesh(), Foam::IOobject::MUST_READ,
                               Foam::IOobject::NO_WRITE},
                geometry.mesh()},
            geometry,
            directories,
            identifier,
            field_name,
            precision} {}

  void operator()(double time) override {
    std::size_t nr_particles = 0;
    for (auto const &part : _subject) {
      auto const &state_old = part.state_old();
      if (adsorbed(state_old) || !brackets_time(part, time)) {
        continue;
      }
      ++nr_particles;

      auto const &state_new = part.state_new();
      _output << std::setw(_column_widths[0]) << state_old.tag
              << std::setw(_column_widths[1])
              << ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new, state_old);
      if (outside(state_new.cell) || outside(state_old.cell)) {
        for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
          io::print(_output, std::vector<double>(Geometry::dim, 0.),
                    _column_widths[2]);
        }
      } else {
        for (std::size_t dd1 = 0; dd1 < Geometry::dim; ++dd1) {
          for (std::size_t dd2 = 0; dd2 < Geometry::dim; ++dd2) {
            auto field_value =
                ctrw::Get_interp{time, ctrw::Get_cell_arraystyleproperty{
                                           _field}}(state_new, state_old);
            _output << std::setw(_column_widths[2])
                    << field_value.row(dd1)[dd2];
          }
        }
      }
      _output << "\n";
    }

    _output_times << std::setw(_column_widths_times[0]) << time
                  << std::setw(_column_widths_times[1]) << nr_particles << "\n";
  }

  void update(double time = 0.,
              Foam::instant const &previous_time_of_change = {},
              Foam::instant const &next_time_of_change = {}) override {
    if constexpr (!std::is_reference_v<Field>) {
      _field = {Foam::IOobject{_field_name, _geometry.mesh().time().timeName(),
                               _geometry.mesh(), Foam::IOobject::MUST_READ,
                               Foam::IOobject::NO_WRITE},
                _geometry.mesh()};
    }
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using MeasurerTime<Subject, Geometry>::_geometry;
  std::string _field_name;
  Field _field;
  std::array<int, 3> _column_widths;
  std::ofstream _output_times;
  std::array<int, 2> _column_widths_times;
};
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_tensor_field(Subject const &, Field &&, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &, int)
    -> MeasurerTime_tensor_field<Subject, Geometry, Field>;
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_tensor_field(Subject const &, Field &&, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &)
    -> MeasurerTime_tensor_field<Subject, Geometry, Field>;
template <typename Subject, typename Geometry>
MeasurerTime_tensor_field(Subject const &, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &, int)
    -> MeasurerTime_tensor_field<Subject, Geometry, Foam::volTensorField>;
template <typename Subject, typename Geometry>
MeasurerTime_tensor_field(Subject const &, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &)
    -> MeasurerTime_tensor_field<Subject, Geometry, Foam::volTensorField>;

/** @brief Output time and mean of scalar field over particles. */
template <typename Subject, typename Geometry,
          typename Field = ScalarField_Interpolation<
              Foam::volScalarField, typename Geometry::Locator const &,
              InterpolationTypes::Linear, CheckOptions::Check>>
class MeasurerTime_scalar_field_mean final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_scalar_field_mean(Subject const &subject, Field &&field,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 std::string const &field_name,
                                 int precision = 8)
      : MeasurerTime<Subject,
                     Geometry>{subject,
                               geometry,
                               directories,
                               std::string{"scalar_field_mean_"} + field_name,
                               identifier,
                               precision},
        _field_name{field_name}, _field{std::forward<Field>(field)},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision, int(2 + (field_name + "_mean").length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << field_name + "_mean" << "\n";
  }

  MeasurerTime_scalar_field_mean(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 std::string const &field_name,
                                 int precision = 8)
      : MeasurerTime_scalar_field_mean{
            subject,
            ScalarField_Interpolation<
                Foam::volScalarField, typename Geometry::Locator const &,
                InterpolationTypes::Linear, CheckOptions::Check>{
                {Foam::IOobject{field_name, geometry.mesh().time().timeName(),
                                geometry.mesh(), Foam::IOobject::MUST_READ,
                                Foam::IOobject::NO_WRITE},
                 geometry.mesh()},
                geometry.locator},
            geometry,
            directories,
            identifier,
            field_name,
            precision} {}

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time
            << std::setw(_column_widths[1]) << mean(_subject, time, _field)
            << "\n";
  }

  void update(double time = 0.,
              Foam::instant const &previous_time_of_change = {},
              Foam::instant const &next_time_of_change = {}) override {
    if constexpr (!std::is_reference_v<Field>) {
      _field.set(
          {Foam::IOobject{_field_name, _geometry.mesh().time().timeName(),
                          _geometry.mesh(), Foam::IOobject::MUST_READ,
                          Foam::IOobject::NO_WRITE},
           _geometry.mesh()});
    }
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using MeasurerTime<Subject, Geometry>::_geometry;
  std::string _field_name;
  Field _field;
  std::array<int, 2> _column_widths;
};
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_scalar_field_mean(Subject const &, Field &&, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &, int)
    -> MeasurerTime_scalar_field_mean<Subject, Geometry, Field>;
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_scalar_field_mean(Subject const &, Field &&, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &)
    -> MeasurerTime_scalar_field_mean<Subject, Geometry, Field>;
template <typename Subject, typename Geometry>
MeasurerTime_scalar_field_mean(Subject const &, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &, int)
    -> MeasurerTime_scalar_field_mean<
        Subject, Geometry,
        ScalarField_Interpolation<
            Foam::volScalarField, typename Geometry::Locator const &,
            InterpolationTypes::Linear, CheckOptions::Check>>;
template <typename Subject, typename Geometry>
MeasurerTime_scalar_field_mean(Subject const &, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &)
    -> MeasurerTime_scalar_field_mean<
        Subject, Geometry,
        ScalarField_Interpolation<
            Foam::volScalarField, typename Geometry::Locator const &,
            InterpolationTypes::Linear, CheckOptions::Check>>;

/** @brief  Output time and mean of vector field over particles. */
template <typename Subject, typename Geometry,
          typename Field = VectorField_Interpolation<
              Foam::volVectorField, typename Geometry::Locator const &,
              InterpolationTypes::Linear, CheckOptions::Check>>
class MeasurerTime_vector_field_mean final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_vector_field_mean(Subject const &subject, Field &&field,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 std::string const &field_name,
                                 int precision = 8)
      : MeasurerTime<Subject,
                     Geometry>{subject,
                               geometry,
                               directories,
                               std::string{"vector_field_mean_"} + field_name,
                               identifier,
                               precision},
        _field{std::forward<Field>(field)},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision,
                     int(2 + (field_name + "_mean_").length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
      _output << std::setw(_column_widths[1])
              << field_name + "_mean_" + std::to_string(dd) << "\n";
    }
  }

  MeasurerTime_vector_field_mean(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 std::string const &field_name,
                                 int precision = 8)
      : MeasurerTime_vector_field_mean{
            subject,
            VectorField_Interpolation<
                Foam::volVectorField, typename Geometry::Locator const &,
                InterpolationTypes::Linear, CheckOptions::Check>{
                {Foam::IOobject{field_name, geometry.mesh().time().timeName(),
                                geometry.mesh(), Foam::IOobject::MUST_READ,
                                Foam::IOobject::NO_WRITE},
                 geometry.mesh()},
                geometry.locator},
            geometry,
            directories,
            identifier,
            field_name,
            precision} {}

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    io::print(_output, mean(_subject, time, _field), _column_widths[1]);
    _output << "\n";
  }

  void update(double time = 0.,
              Foam::instant const &previous_time_of_change = {},
              Foam::instant const &next_time_of_change = {}) override {
    if constexpr (!std::is_reference_v<Field>) {
      _field.set(
          {Foam::IOobject{_field_name, _geometry.mesh().time().timeName(),
                          _geometry.mesh(), Foam::IOobject::MUST_READ,
                          Foam::IOobject::NO_WRITE},
           _geometry.mesh()});
    }
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using MeasurerTime<Subject, Geometry>::_geometry;
  std::string _field_name;
  Field _field;
  std::array<int, 3> _column_widths;
};
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_vector_field_mean(Subject const &, Field &&, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &, int)
    -> MeasurerTime_vector_field_mean<Subject, Geometry, Field>;
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_vector_field_mean(Subject const &, Field &&, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &)
    -> MeasurerTime_vector_field_mean<Subject, Geometry, Field>;
template <typename Subject, typename Geometry>
MeasurerTime_vector_field_mean(Subject const &, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &, int)
    -> MeasurerTime_vector_field_mean<
        Subject, Geometry,
        VectorField_Interpolation<
            Foam::volVectorField, typename Geometry::Locator const &,
            InterpolationTypes::Linear, CheckOptions::Check>>;
template <typename Subject, typename Geometry>
MeasurerTime_vector_field_mean(Subject const &, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &)
    -> MeasurerTime_vector_field_mean<
        Subject, Geometry,
        VectorField_Interpolation<
            Foam::volVectorField, typename Geometry::Locator const &,
            InterpolationTypes::Linear, CheckOptions::Check>>;

/** @brief Output time and mean of tensor field over particles. */
template <typename Subject, typename Geometry,
          typename Field = Foam::volTensorField>
class MeasurerTime_tensor_field_mean final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_tensor_field_mean(Subject const &subject, Field &&field,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 std::string const &field_name,
                                 int precision = 8)
      : MeasurerTime<Subject,
                     Geometry>{subject,
                               geometry,
                               directories,
                               std::string{"tensor_field_mean_"} + field_name,
                               identifier,
                               precision},
        _field{std::forward<Field>(field)},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision,
                     int(3 + (field_name + "_mean_").length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd1 = 0; dd1 < Geometry::dim; ++dd1) {
      for (std::size_t dd2 = 0; dd2 < Geometry::dim; ++dd2) {
        _output << std::setw(_column_widths[2])
                << field_name + "_mean_" + std::to_string(dd1) +
                       std::to_string(dd2);
      }
    }
    _output << "\n";
  }

  MeasurerTime_tensor_field_mean(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 std::string const &field_name,
                                 int precision = 8)
      : MeasurerTime_tensor_field_mean{
            subject,
            Foam::volTensorField{
                Foam::IOobject{field_name, geometry.mesh().time().timeName(),
                               geometry.mesh(), Foam::IOobject::MUST_READ,
                               Foam::IOobject::NO_WRITE},
                geometry.mesh()},
            geometry,
            directories,
            identifier,
            field_name,
            precision} {}

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    auto field_mean =
        mean(_subject, time, ctrw::Get_cell_arraystyleproperty{_field});
    for (std::size_t dd1 = 0; dd1 < Geometry::dim; ++dd1) {
      for (std::size_t dd2 = 0; dd2 < Geometry::dim; ++dd2) {
        _output << std::setw(_column_widths[1]) << field_mean.row(dd1)[dd2];
      }
    }
    _output << "\n";
  }

  void update(double time = 0.,
              Foam::instant const &previous_time_of_change = {},
              Foam::instant const &next_time_of_change = {}) override {
    if constexpr (!std::is_reference_v<Field>) {
      _field = {Foam::IOobject{_field_name, _geometry.mesh().time().timeName(),
                               _geometry.mesh(), Foam::IOobject::MUST_READ,
                               Foam::IOobject::NO_WRITE},
                _geometry.mesh()};
    }
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using MeasurerTime<Subject, Geometry>::_geometry;
  using MeasurerTime<Subject, Geometry>::_locator;
  using State = typename Subject::Particle::State;
  std::string _field_name;
  Field _field;
  std::array<int, 2> _column_widths;
};
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_tensor_field_mean(Subject const &, Field &&, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &, int)
    -> MeasurerTime_tensor_field_mean<Subject, Geometry, Field>;
template <typename Subject, typename Geometry, typename Field>
MeasurerTime_tensor_field_mean(Subject const &, Field &&, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &)
    -> MeasurerTime_tensor_field_mean<Subject, Geometry, Field>;
template <typename Subject, typename Geometry>
MeasurerTime_tensor_field_mean(Subject const &, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &, int)
    -> MeasurerTime_tensor_field_mean<Subject, Geometry, Foam::volTensorField>;
template <typename Subject, typename Geometry>
MeasurerTime_tensor_field_mean(Subject const &, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &)
    -> MeasurerTime_tensor_field_mean<Subject, Geometry, Foam::volTensorField>;

/** @brief Time and net reacted mass at each boundary face. */
template <typename Subject, typename Geometry>
class MeasurerTime_first_crossing_time final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_first_crossing_time(Subject const &subject,
                                   Geometry const &geometry,
                                   Directories const &directories,
                                   std::string const &identifier,
                                   std::size_t dimension, double position,
                                   int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,
                                        geometry,
                                        directories,
                                        std::string{"first_crossing_time_"} +
                                            std::to_string(dimension) + "_" +
                                            std::to_string(position),
                                        identifier,
                                        precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length()))},
        dimension{dimension}, position{position} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Tag"
            << std::setw(_column_widths[2]) << "Mass" << "\n";
  }

  void operator()(double time) override {
    for (auto const &part : _subject) {
      auto const &state_old = part.state_old();
      if (adsorbed(state_old) || !brackets_time(part, time)) {
        continue;
      }
      auto const &state_new = part.state_new();
      auto pos_new = ctrw::Get_position_component_arg{dimension}(state_new);
      auto pos_old = ctrw::Get_position_component_arg{dimension}(state_old);
      if ((pos_new - position) * (pos_old - position) <= 0.) {
        Time crossing_time = ctrw::Get_time_interp{
            position,
            ctrw::Get_position_component_arg{dimension}}(state_new, state_old);
        Mass crossing_mass = ctrw::Get_interp{crossing_time, ctrw::Get_mass{}}(
            state_new, state_old);
        _crossing_times_masses.insert(std::pair<Tag, std::pair<Time, Mass>>{
            state_old.tag, {crossing_time, crossing_mass}});
      }
    }
  }

  void print() override {
    for (auto const &tag_time_mass : _crossing_times_masses) {
      _output << std::setw(_column_widths[0]) << tag_time_mass.second.first
              << std::setw(_column_widths[1]) << tag_time_mass.first
              << std::setw(_column_widths[2]) << tag_time_mass.second.second
              << "\n";
    }
  }

  auto const &crossing_times_masses() { return crossing_times_masses; }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 3> _column_widths;

  using Time = decltype(Subject::Particle::State::time);
  using Mass = decltype(Subject::Particle::State::mass);
  using Tag = decltype(Subject::Particle::State::tag);
  std::size_t dimension;
  double position;
  std::map<Tag, std::pair<Time, Mass>> _crossing_times_masses;
};

/**
 * @brief Output tags, masses, and positions of adsorbed particles, as well as
 *        times and number of particles output per time to a separate file (with
 *        \c _times added to output name).
 *
 * @tparam periodic If \c true, give position accounting for periodicity
 *                  (possibly outside domain).
 */
template <typename Subject, typename Geometry, bool periodic = false>
class MeasurerTime_position_adsorbed final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_position_adsorbed(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 int precision = 8)
      : MeasurerTime<
            Subject,
            Geometry>{subject,
                      geometry,
                      directories,
                      std::string{"position_adsorbed"} + (periodic ? "_periodic" : ""),
                      identifier,
                      precision},
        _column_widths{
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())),
            std::max(9 + precision,
                     int(2 + std::string{"Position_"}.length()))},
        _output_times{this->open_write(
            directories, this->_output_name + "_times", identifier)},
        _column_widths_times{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Nr_particles"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Tag"
            << std::setw(_column_widths[1]) << "Mass";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
      _output << std::setw(_column_widths[2])
              << "Position_" + std::to_string(dd);
    }

    _output_times << std::setw(_column_widths_times[0]) << "Time"
                  << std::setw(_column_widths_times[1]) << "Nr_particles"
                  << "\n";
  }

  void operator()(double time) override {
    std::size_t nr_particles = 0;
    for (auto const &part : _subject) {
      auto const &state_old = part.state_old();
      if (adsorbed(state_old) || !brackets_time(part, time)) {
        continue;
      }
      ++nr_particles;

      auto const &state_new = part.state_old();
      _output << std::setw(_column_widths[0]) << state_old.tag
              << std::setw(_column_widths[1])
              << ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new, state_old);
      io::print(
          _output,
          ctrw::getter_position(boundary_periodic(this->_geometry))(state_old),
          _column_widths[2])
          << "\n";
    }

    _output_times << std::setw(_column_widths_times[0]) << time
                  << std::setw(_column_widths_times[1]) << nr_particles << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 3> _column_widths;
  std::ofstream _output_times;
  std::array<int, 2> _column_widths_times;
};

/**
 * @brief Output boundary faces and corresponding net reacted masses, as well as
 *        times and number of particles output per time to a separate file (with
 *        \c _times added to output name).
 *
 * @tparam periodic If \c true, give face accounting for periodicity (possibly
 *                  outside domain).
 */
template <typename Subject, typename Geometry, bool periodic = false>
class MeasurerTime_mass_reacted_face final
    : public MeasurerTime<Subject, Geometry> {
public:
  using BoundaryInfo = std::conditional_t<
      periodic,
      BoundaryInfo_Record_mass_reacted_face_periodic<
          typename Subject::Particle::State>,
      BoundaryInfo_Record_mass_reacted_face<typename Subject::Particle::State>>;

  MeasurerTime_mass_reacted_face(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 BoundaryInfo const &boundary_info,
                                 int precision = 8)
      : MeasurerTime<
            Subject,
            Geometry>{subject,
                      geometry,
                      directories,
                      std::string{"mass_reacted_face"} + (periodic ? "_periodic" : ""),
                      identifier,
                      precision},
        _boundary_info{boundary_info},
        _column_widths{
            std::max(12, int(1 + std::string{"Face"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())),
            std::max(12, int(2 + std::string{"Periodicity_"}.length()))},
        _output_times{this->open_write(
            directories, this->_output_name + "_times", identifier)},
        _column_widths_times{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Nr_faces"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Face"
            << std::setw(_column_widths[0]) << "Mass";
    if constexpr (periodic) {
      for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
        _output << std::setw(_column_widths[2])
                << "Periodicity_" + std::to_string(dd);
      }
    }
    _output << "\n";

    _output_times << std::setw(_column_widths_times[0]) << "Time"
                  << std::setw(_column_widths_times[1]) << "Nr_faces" << "\n";
  }

  template <bool periodicity_v>
  MeasurerTime_mass_reacted_face(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 BoundaryInfo const &boundary_info,
                                 meta::Selector<bool, periodicity_v>,
                                 int precision = 8)
      : MeasurerTime_mass_reacted_face<Subject, Geometry, periodicity_v>{
            subject, geometry, directories, identifier, precision} {}

  MeasurerTime_mass_reacted_face(Subject const &, Geometry const &,
                                 Directories const &, std::string const &,
                                 BoundaryInfo &&, int precision = 8) = delete;

  template <bool periodicity_v>
  MeasurerTime_mass_reacted_face(Subject const &, Geometry const &,
                                 Directories const &, std::string const &,
                                 BoundaryInfo &&,
                                 meta::Selector<bool, periodicity_v>,
                                 int precision = 8) = delete;

  void operator()(double time) override {
    if constexpr (periodic) {
      for (auto const &faceperiodicity_mass : _boundary_info.masses()) {
        _output << std::setw(_column_widths[0])
                << faceperiodicity_mass.first.first // Face
                << std::setw(_column_widths[1])
                << faceperiodicity_mass.second; // Face mass
        io::print(_output, faceperiodicity_mass.first.second, _column_widths[2])
            << "\n";
      }
    } else {
      for (auto const &face_mass : _boundary_info.masses()) {
        _output << std::setw(_column_widths[0]) << face_mass.first
                << std::setw(_column_widths[1]) << face_mass.second << "\n";
      }
    }

    _output_times << std::setw(_column_widths_times[0]) << time
                  << std::setw(_column_widths_times[1])
                  << _boundary_info.masses().size() << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  BoundaryInfo const &_boundary_info;
  std::array<int, 3> _column_widths;
  std::ofstream _output_times;
  std::array<int, 2> _column_widths_times;
};

/**
 * @brief Output boundary faces and corresponding adsorbed masses, as well as
 *        times and number of particles output per time to a separate file (with
 *        \c _times added to output name).
 *
 * @tparam periodic If \c true, give face accounting for periodicity (possibly
 *                  outside domain).
 */
template <typename Subject, typename Geometry, bool periodic = false>
class MeasurerTime_mass_adsorbed_face final
    : public MeasurerTime<Subject, Geometry> {
public:
  MeasurerTime_mass_adsorbed_face(Subject const &subject,
                                  Geometry const &geometry,
                                  Directories const &directories,
                                  std::string const &identifier,
                                  int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,
                                        geometry,
                                        directories,
                                        std::string{"mass_adsorbed_face"} +
                                            (periodic ? "_periodi"
                                                        "c"
                                                      : ""),
                                        identifier,
                                        precision},
        _column_widths{
            std::max(12, int(1 + std::string{"Face"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())),
            std::max(12, int(2 + std::string{"Periodicity_"}.length()))},
        _output_times{this->open_write(
            directories, this->_output_name + "_times", identifier)},
        _column_widths_times{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Nr_faces"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Face"
            << std::setw(_column_widths[0]) << "Mass";
    if constexpr (periodic) {
      for (std::size_t dd = 0; dd < Geometry::dim; ++dd) {
        _output << std::setw(_column_widths[2])
                << "Periodicity_" + std::to_string(dd);
      }
    }
    _output << "\n";

    _output_times << std::setw(_column_widths_times[0]) << "Time"
                  << std::setw(_column_widths_times[1]) << "Nr_faces" << "\n";
  }

  template <bool periodicity_v>
  MeasurerTime_mass_adsorbed_face(Subject const &subject,
                                  Geometry const &geometry,
                                  Directories const &directories,
                                  std::string const &identifier,
                                  meta::Selector<bool, periodicity_v>,
                                  int precision = 8)
      : MeasurerTime_mass_adsorbed_face<Subject, Geometry, periodicity_v>{
            subject, geometry, directories, identifier, precision} {}

  void operator()(double time) override {
    auto masses = adsorbed_masses(time);
    if constexpr (periodic) {
      for (auto const &faceperiodicity_mass : masses) {
        _output << std::setw(_column_widths[0])
                << faceperiodicity_mass.first.first // Face
                << std::setw(_column_widths[1])
                << faceperiodicity_mass.second; // Mass
        // Face periodicity
        io::print(_output, faceperiodicity_mass.first.second,
                  _column_widths[2]);
        _output << "\n";
      }
    } else {
      for (auto const &face_mass : masses) {
        _output << std::setw(_column_widths[1]) << face_mass.first
                << std::setw(_column_widths[2]) << face_mass.second << "\n";
      }
    }

    _output_times << std::setw(_column_widths_times[0]) << time
                  << std::setw(_column_widths_times[1]) << masses.size()
                  << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 4> _column_widths;
  std::ofstream _output_times;
  std::array<int, 2> _column_widths_times;

  using State = typename Subject::Particle::State;
  template <typename TT, bool periodic_v> struct PeriodicityOrEmpty {
    using type = meta::Empty;
  };
  template <typename TT> struct PeriodicityOrEmpty<TT, true> {
    using type = typename TT::Periodicity;
  };
  template <typename TT, bool periodic_v>
  using PeriodicityOrEmpty_t =
      typename PeriodicityOrEmpty<TT, periodic_v>::type;

  using Masses = std::conditional_t<
      periodic,
      std::unordered_map<
          std::pair<Foam::label, PeriodicityOrEmpty_t<State, periodic>>, double,
          useful::Hash_pair<Foam::label,
                            PeriodicityOrEmpty_t<State, periodic>>>,
      std::unordered_map<Foam::label, double>>;

  auto insert(Masses &masses, Foam::label face, State const &state,
              double mass_change) const {
    if constexpr (periodic) {
      return masses.insert({{face, state.periodicity}, mass_change});
    } else {
      return masses.insert({face, mass_change});
    }
  }

  auto adsorbed_masses(double time) const {
    Masses masses;

    for (auto const &part : _subject) {
      auto const &state_old = part.state_old();
      if (!adsorbed(state_old) || !brackets_time(part, time)) {
        continue;
      }
      auto const &state_new = part.state_new();
      Foam::label face;
      if constexpr (meta::has_boundary_face_v<
                        typename Subject::Particle::State::Info>) {
        face = state_old.info.boundary_face;
      } else {
        face = this->_geometry.mesh_search().findNearestBoundaryFace(
            make_point(state_old.position));
      }
      double mass_change =
          ctrw::Get_interp{time, ctrw::Get_mass{}}(state_new, state_old);
      auto it_inserted = insert(masses, face, state_old, mass_change);
      if (it_inserted.second) {
        it_inserted.first->second += mass_change;
      }
    }

    return masses;
  }
};
} // namespace ptof

#endif /* PTOF_MEASURERTIME_H */
