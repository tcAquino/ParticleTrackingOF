/**
 \file PTOF/MeasurerTime.h
 \author Tomás Aquino
 \date 07/03/2022
*/

#ifndef PTOF_MEASURERTIME_H
#define PTOF_MEASURERTIME_H

#include "CTRW/StateGetter.h"
#include "General/Useful.h"
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
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <volFieldsFwd.H>

namespace ptof {
/** \class MeasurerTime PTOF/Output.h "PTOF/Output.h"
 \brief Polymorphic output handler for outputting time and some quantity.
 \note To be used only through derived classes. */
template <typename Subject, typename Geometry> struct MeasurerTime {
  /** Destructor. */
  virtual ~MeasurerTime() { _output.close(); }

  /** \brief Make measurement and output or store, given current time \c time.
   */
  virtual void operator()(double time) = 0;

  /** \brief Output stored information (do nothing if not overriden). */
  virtual void print(){};

  /** \brief Update internal state. */
  virtual void update(double time, double time_of_change) {}

protected:
  /** \brief Set up output streams for requested output types.
  \param subject CTRW object to measure.
  \param geometry Domain geometry info and utilities.
  \param directories Current case directory information.
  \param output_name Name of output type.
  \param identifier String to include in names of output files.
  \param precision Number of digits after decimal point in output, in scientific
  notation. */
  MeasurerTime(Subject const &subject, Geometry const &geometry,
               Directories const &directories, std::string const &output_name,
               std::string const &identifier, int precision = 8)
      : _subject{subject}, _geometry{geometry}, _locator{geometry.locator},
        _output{open_write(directories, output_name, identifier)} {
    _output << std::setprecision(precision) << std::scientific;
  }

  Subject const &_subject;   /**< CTRW object to measure. */
  Geometry const &_geometry; /**< Domain geometry info and utilities. */
  typename Geometry::Locator const
      &_locator;         /**< Object to locate positions in mesh. */
  std::ofstream _output; /**< Output stream. */

  /** \brief Open output file for a given output type.
  \param directories Current case directory information.
  \param output_name Name of output type.
  \param identifier String to include in names of output files.
  \return Output stream. */
  std::ofstream open_write(Directories const &directories,
                           std::string const &output_name,
                           std::string const &identifier) {
    return useful::open_write(directories.dir_output + "/Data_" + output_name +
                              (identifier.empty() ? "" : "_" + identifier) +
                              ".dat");
  }
};

/** \class MeasurerTime_position PTOF/Output.h "PTOF/Output.h"
 *  \brief Output time, tags, positions, and masses. */
template <typename Subject, typename Geometry>
struct MeasurerTime_position final : MeasurerTime<Subject, Geometry> {
  MeasurerTime_position(Subject const &subject, Geometry const &geometry,
                        Directories const &directories,
                        std::string const &identifier, int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,    geometry,   directories,
                                        "position", identifier, precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(2 + std::string{"Position_"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Tag";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[2])
              << "Position_" + std::to_string(dd);
    _output << std::setw(_column_widths[3]) << "Mass"
            << std::setw(_column_widths[1]) << "Tag"
            << std::setw(_column_widths[2]) << "..."
            << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    for (auto const &part : _subject.particles()) {
      auto const &state = part.state_new();
      if (!state.info.absorbed && part.state_old().time <= time) {
        _output << std::setw(_column_widths[1]) << state.tag;
        useful::print(_output, state.position, _column_widths[2]);
        _output << std::setw(_column_widths[3]) << state.mass;
      }
    }
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 4> _column_widths;
};

/** \class MeasurerTime_position_in_regions PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time, tags, positions, and masses within regions specified by
 * masks. */
template <typename Subject, typename Geometry, typename Mask>
struct MeasurerTime_position_in_regions final
    : MeasurerTime<Subject, Geometry> {
  MeasurerTime_position_in_regions(
      Subject const &subject, Geometry const &geometry,
      Directories const &directories, std::string const &identifier,
      std::vector<std::reference_wrapper<const Mask>> masks,
      std::vector<double> thresholds = {}, int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,     geometry,
                                        directories, "position_in_regions",
                                        identifier,  precision},
        _masks{masks}, _thresholds{thresholds},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(2 + std::string{"Position_"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())),
            std::max(9, int(4 + std::string{"In_region_"}.length()))} {
    _thresholds.resize(masks.size(), 0.);
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Tag";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[2])
              << "Position_" + std::to_string(dd);
    _output << std::setw(_column_widths[3]) << "Mass";
    for (std::size_t ii = 0; ii < _masks.size(); ++ii)
      _output << std::setw(_column_widths[4])
              << "In_region_" + std::to_string(ii);
    _output << std::setw(_column_widths[1]) << "Tag"
            << std::setw(_column_widths[2]) << "..."
            << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    for (auto const &part : _subject.particles()) {
      auto const &state = part.state_new();
      _output << std::setw(_column_widths[1]) << state.tag;
      if (!state.info.absorbed && part.state_old().time <= time &&
          !outside(state.cell)) {
        std::vector<int> in_region(_masks.size(), 0);
        for (std::size_t ii = 0; ii < _masks.size(); ++ii)
          if (_masks[ii].get()[state.cell] > _thresholds[ii])
            in_region[ii] = 1;
        if (std::any_of(in_region.begin(), in_region.end(),
                        [](int ii) { return ii > 0; })) {
          useful::print(_output, state.position, _column_widths[2]);
          _output << std::setw(_column_widths[3]) << state.mass;
          useful::print(_output, in_region, _column_widths[4]);
        }
      }
    }
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::vector<std::reference_wrapper<const Mask>> _masks;
  std::vector<double> _thresholds;
  std::array<int, 5> _column_widths;
};

/** \class MeasurerTime_position_mean PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time and mean position (weighted by mass). */
template <typename Subject, typename Geometry>
struct MeasurerTime_position_mean final : MeasurerTime<Subject, Geometry> {
  MeasurerTime_position_mean(Subject const &subject, Geometry const &geometry,
                             Directories const &directories,
                             std::string const &identifier, int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,     geometry,
                                        directories, "position_mean",
                                        identifier,  precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(2 + std::string{"Position_mean_"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[1])
              << "Position_mean_" + std::to_string(dd);
    _output << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    useful::print(_output, position_mean(_subject, time), _column_widths[1]);
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 2> _column_widths;
};

/** \class MeasurerTime_position_second_moment PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time and second moment of position (weighted by mass). */
template <typename Subject, typename Geometry>
struct MeasurerTime_position_second_moment final
    : MeasurerTime<Subject, Geometry> {
  MeasurerTime_position_second_moment(Subject const &subject,
                                      Geometry const &geometry,
                                      Directories const &directories,
                                      std::string const &identifier,
                                      int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,     geometry,
                                        directories, "position_second_moment",
                                        identifier,  precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(
                9 + precision,
                int(2 + std::string{"Position_second_moment_"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[1])
              << "Position_second_moment_" + std::to_string(dd);
    _output << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    useful::print(_output, position_second_moment(_subject, time),
                  _column_widths[1]);
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 2> _column_widths;
};

/** \class MeasurerTime_position_variance PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time and position variance (weighted by mass). */
template <typename Subject, typename Geometry>
struct MeasurerTime_position_variance final : MeasurerTime<Subject, Geometry> {
  MeasurerTime_position_variance(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,     geometry,
                                        directories, "position_variance",
                                        identifier,  precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision,
                     int(2 + std::string{"Position_variance_"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[1])
              << "Position_variance_" + std::to_string(dd);
    _output << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    useful::print(_output, position_variance(_subject, time),
                  _column_widths[1]);
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  std::array<int, 2> _column_widths;
};

/** \class MeasurerTime_mass PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time and total mass. */
template <typename Subject, typename Geometry>
struct MeasurerTime_mass final : MeasurerTime<Subject, Geometry> {
  MeasurerTime_mass(Subject const &subject, Geometry const &geometry,
                    Directories const &directories,
                    std::string const &identifier, int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject, geometry,   directories,
                                        "mass",  identifier, precision},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Mass"
            << "\n";
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

/** \class MeasurerTime_mass_in_regions PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time and total mass in regions specified by maks. */
template <typename Subject, typename Geometry, typename Mask>
struct MeasurerTime_mass_in_regions final : MeasurerTime<Subject, Geometry> {
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
    for (std::size_t ii = 0; ii < masks.size(); ++ii)
      _output << std::setw(_column_widths[4])
              << "Mass_region_" + std::to_string(ii);
    _output << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    useful::print(_output, mass(_subject, time, _locator, _masks, _thresholds),
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
\class MeasurerTime_scalar_field PTOF/Output.h "PTOF/Output.h"
\brief Output time, tags, and scalar field values.
*/
template <typename Subject, typename Geometry,
          typename Field = ScalarField_LinearInterpolation_OF<
              Foam::volScalarField, typename Geometry::Locator const &,
              CheckOptions::Check>>
struct MeasurerTime_scalar_field final : MeasurerTime<Subject, Geometry> {
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
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(2 + (field_name + "_").length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Tag"
            << std::setw(_column_widths[2]) << field_name
            << std::setw(_column_widths[1]) << "Tag"
            << std::setw(_column_widths[2]) << "..."
            << "\n";
  }

  MeasurerTime_scalar_field(Subject const &subject, Geometry const &geometry,
                            Directories const &directories,
                            std::string const &identifier,
                            std::string const &field_name, int precision = 8)
      : MeasurerTime_scalar_field{
            subject,
            ScalarField_LinearInterpolation_OF<
                Foam::volScalarField, typename Geometry::Locator const &,
                CheckOptions::Check>{
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
    for (auto const &part : _subject.particles()) {
      auto const &state = part.state_new();
      if (!state.info.absorbed && part.state_old().time <= time)
        _output << std::setw(_column_widths[1]) << state.tag
                << std::setw(_column_widths[2]) << _field(state);
    }
    _output << "\n";
  }

  void update(double time, double time_of_change) override {
    if constexpr (!std::is_reference_v<Field>)
      if (time >= time_of_change)
        _field.set(
            {Foam::IOobject{_field_name, _geometry.mesh().time().timeName(),
                            _geometry.mesh(), Foam::IOobject::MUST_READ,
                            Foam::IOobject::NO_WRITE},
             _geometry.mesh()});
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
        ScalarField_LinearInterpolation_OF<Foam::volScalarField,
                                           typename Geometry::Locator const &,
                                           CheckOptions::Check>>;
template <typename Subject, typename Geometry>
MeasurerTime_scalar_field(Subject const &, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &)
    -> MeasurerTime_scalar_field<
        Subject, Geometry,
        ScalarField_LinearInterpolation_OF<Foam::volScalarField,
                                           typename Geometry::Locator const &,
                                           CheckOptions::Check>>;

/** \class MeasurerTime_vector_field PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time, tags, and vector field values. */
template <typename Subject, typename Geometry,
          typename Field = VectorField_LinearInterpolation_OF<
              Foam::volVectorField, typename Geometry::Locator const &,
              CheckOptions::Check>>
struct MeasurerTime_vector_field final : MeasurerTime<Subject, Geometry> {
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
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(2 + (field_name + "_").length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Tag";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[2])
              << field_name + "_" + std::to_string(dd);
    _output << std::setw(_column_widths[1]) << "Tag"
            << std::setw(_column_widths[2]) << "..."
            << "\n";
  }

  MeasurerTime_vector_field(Subject const &subject, Geometry const &geometry,
                            Directories const &directories,
                            std::string const &identifier,
                            std::string const &field_name, int precision = 8)
      : MeasurerTime_vector_field{
            subject,
            VectorField_LinearInterpolation_OF<
                Foam::volVectorField, typename Geometry::Locator const &,
                CheckOptions::Check>{
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
    for (auto const &part : _subject.particles()) {
      auto const &state = part.state_new();
      if (!state.info.absorbed && part.state_old().time <= time) {
        _output << std::setw(_column_widths[1]) << state.tag;
        useful::print(_output, _field(state), _column_widths[2]);
      }
    }
    _output << "\n";
  }

  void update(double time, double time_of_change) override {
    if constexpr (!std::is_reference_v<Field>)
      if (time >= time_of_change)
        _field.set(
            {Foam::IOobject{_field_name, _geometry.mesh().time().timeName(),
                            _geometry.mesh(), Foam::IOobject::MUST_READ,
                            Foam::IOobject::NO_WRITE},
             _geometry.mesh()});
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
        VectorField_LinearInterpolation_OF<Foam::volVectorField,
                                           typename Geometry::Locator const &,
                                           CheckOptions::Check>>;
template <typename Subject, typename Geometry>
MeasurerTime_vector_field(Subject const &, Geometry const &,
                          Directories const &, std::string const &,
                          std::string const &)
    -> MeasurerTime_vector_field<
        Subject, Geometry,
        VectorField_LinearInterpolation_OF<Foam::volVectorField,
                                           typename Geometry::Locator const &,
                                           CheckOptions::Check>>;

/** \class MeasurerTime_tensor_field PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time, tag, and tensor field values. */
template <typename Subject, typename Geometry,
          typename Field = Foam::volTensorField>
struct MeasurerTime_tensor_field final : MeasurerTime<Subject, Geometry> {
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
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(3 + (field_name + "_").length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Tag";
    for (std::size_t dd1 = 0; dd1 < Geometry::dim; ++dd1)
      for (std::size_t dd2 = 0; dd2 < Geometry::dim; ++dd2)
        _output << std::setw(_column_widths[2])
                << field_name + "_" + std::to_string(dd1) + std::to_string(dd2);
    _output << std::setw(_column_widths[1]) << "Tag"
            << std::setw(_column_widths[2]) << "..."
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
    _output << std::setw(_column_widths[0]) << time;
    for (auto const &part : _subject.particles()) {
      auto const &state = part.state_new();
      if (!state.info.absorbed && part.state_old().time <= time) {
        _output << std::setw(_column_widths[1]) << state.tag;
        if (outside(state.cell))
          for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
            useful::print(_output, std::vector<double>(Geometry::dim, 0.),
                          _column_widths[2]);
        else
          for (std::size_t dd1 = 0; dd1 < Geometry::dim; ++dd1)
            for (std::size_t dd2 = 0; dd2 < Geometry::dim; ++dd2)
              _output << std::setw(_column_widths[2])
                      << _field[state.cell].row(dd1)[dd2];
      }
    }
    _output << "\n";
  }

  void update(double time, double time_of_change) override {
    if constexpr (!std::is_reference_v<Field>)
      if (time >= time_of_change)
        _field = {Foam::IOobject{_field_name,
                                 _geometry.mesh().time().timeName(),
                                 _geometry.mesh(), Foam::IOobject::MUST_READ,
                                 Foam::IOobject::NO_WRITE},
                  _geometry.mesh()};
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

/**
\class MeasurerTime_scalar_field_mean PTOF/Output.h "PTOF/Output.h"
\brief Output time and mean of scalar field over particles.
*/
template <typename Subject, typename Geometry,
          typename Field = ScalarField_LinearInterpolation_OF<
              Foam::volScalarField, typename Geometry::Locator const &,
              CheckOptions::Check>>
struct MeasurerTime_scalar_field_mean final : MeasurerTime<Subject, Geometry> {
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
                               std::string{"scalar_field_"} + field_name,
                               identifier,
                               precision},
        _field_name{field_name}, _field{std::forward<Field>(field)},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision, int(2 + (field_name + "_mean").length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << field_name + "_mean"
            << "\n";
  }

  MeasurerTime_scalar_field_mean(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 std::string const &field_name,
                                 int precision = 8)
      : MeasurerTime_scalar_field_mean{
            subject,
            ScalarField_LinearInterpolation_OF<
                Foam::volScalarField, typename Geometry::Locator const &,
                CheckOptions::Check>{
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

  void update(double time, double time_of_change) override {
    if constexpr (!std::is_reference_v<Field>)
      if (time >= time_of_change)
        _field.set(
            {Foam::IOobject{_field_name, _geometry.mesh().time().timeName(),
                            _geometry.mesh(), Foam::IOobject::MUST_READ,
                            Foam::IOobject::NO_WRITE},
             _geometry.mesh()});
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
        ScalarField_LinearInterpolation_OF<Foam::volScalarField,
                                           typename Geometry::Locator const &,
                                           CheckOptions::Check>>;
template <typename Subject, typename Geometry>
MeasurerTime_scalar_field_mean(Subject const &, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &)
    -> MeasurerTime_scalar_field_mean<
        Subject, Geometry,
        ScalarField_LinearInterpolation_OF<Foam::volScalarField,
                                           typename Geometry::Locator const &,
                                           CheckOptions::Check>>;

/** \class MeasurerTime_vector_field_mean PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time and mean of vector field over particles. */
template <typename Subject, typename Geometry,
          typename Field = VectorField_LinearInterpolation_OF<
              Foam::volVectorField, typename Geometry::Locator const &,
              CheckOptions::Check>>
struct MeasurerTime_vector_field_mean final : MeasurerTime<Subject, Geometry> {
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
                               std::string{"vector_field_"} + field_name,
                               identifier,
                               precision},
        _field{std::forward<Field>(field)},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision,
                     int(2 + (field_name + "_mean_").length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[1])
              << field_name + "_mean_" + std::to_string(dd) << "\n";
  }

  MeasurerTime_vector_field_mean(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 std::string const &field_name,
                                 int precision = 8)
      : MeasurerTime_vector_field_mean{
            subject,
            VectorField_LinearInterpolation_OF<
                Foam::volVectorField, typename Geometry::Locator const &,
                CheckOptions::Check>{
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
    useful::print(_output, mean(_subject, time, _field), _column_widths[1]);
    _output << "\n";
  }

  void update(double time, double time_of_change) override {
    if constexpr (!std::is_reference_v<Field>)
      if (time >= time_of_change)
        _field.set(
            {Foam::IOobject{_field_name, _geometry.mesh().time().timeName(),
                            _geometry.mesh(), Foam::IOobject::MUST_READ,
                            Foam::IOobject::NO_WRITE},
             _geometry.mesh()});
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
        VectorField_LinearInterpolation_OF<Foam::volVectorField,
                                           typename Geometry::Locator const &,
                                           CheckOptions::Check>>;
template <typename Subject, typename Geometry>
MeasurerTime_vector_field_mean(Subject const &, Geometry const &,
                               Directories const &, std::string const &,
                               std::string const &)
    -> MeasurerTime_vector_field_mean<
        Subject, Geometry,
        VectorField_LinearInterpolation_OF<Foam::volVectorField,
                                           typename Geometry::Locator const &,
                                           CheckOptions::Check>>;

/** \class MeasurerTime_tensor_field_mean PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time and mean of tensor field over particles. */
template <typename Subject, typename Geometry,
          typename Field = Foam::volTensorField>
struct MeasurerTime_tensor_field_mean final : MeasurerTime<Subject, Geometry> {
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
                               std::string{"tensor_field_"} + field_name,
                               identifier,
                               precision},
        _field{std::forward<Field>(field)},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision,
                     int(3 + (field_name + "_mean_").length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd1 = 0; dd1 < Geometry::dim; ++dd1)
      for (std::size_t dd2 = 0; dd2 < Geometry::dim; ++dd2)
        _output << std::setw(_column_widths[2])
                << field_name + "_mean_" + std::to_string(dd1) +
                       std::to_string(dd2);
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
    auto field_mean = mean(_subject, time, [this](State const &state) {
      return this->_field[state.cell];
    });
    for (std::size_t dd1 = 0; dd1 < Geometry::dim; ++dd1)
      for (std::size_t dd2 = 0; dd2 < Geometry::dim; ++dd2)
        _output << std::setw(_column_widths[1]) << field_mean.row(dd1)[dd2];
    _output << "\n";
  }

  void update(double time, double time_of_change) override {
    if constexpr (!std::is_reference_v<Field>)
      if (time >= time_of_change)
        _field = {Foam::IOobject{_field_name,
                                 _geometry.mesh().time().timeName(),
                                 _geometry.mesh(), Foam::IOobject::MUST_READ,
                                 Foam::IOobject::NO_WRITE},
                  _geometry.mesh()};
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

/** \class MeasurerTime_position_periodic PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time, tags, positions, and masses, with positions accounting
 * for periodicity. */
template <typename Subject, typename Geometry>
struct MeasurerTime_position_periodic final : MeasurerTime<Subject, Geometry> {
  MeasurerTime_position_periodic(Subject const &subject,
                                 Geometry const &geometry,
                                 Directories const &directories,
                                 std::string const &identifier,
                                 int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,     geometry,
                                        directories, "position_periodic",
                                        identifier,  precision},
        _getter_position{geometry.boundary_periodic},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(3 + std::string{"Position_"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Tag";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[2])
              << "Position_" + std::to_string(dd);
    _output << std::setw(_column_widths[3]) << "Mass"
            << std::setw(_column_widths[1]) << "Tag"
            << std::setw(_column_widths[2]) << "..."
            << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    for (auto const &part : _subject.particles()) {
      auto const &state = part.state_new();
      if (!state.info.absorbed && part.state_old().time <= time) {
        _output << std::setw(_column_widths[1]) << state.tag;
        useful::print(_output, _getter_position(state), _column_widths[2]);
        _output << std::setw(_column_widths[3]) << state.mass;
      }
    }
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using Boundary =
      std::decay_t<decltype(std::declval<Geometry>().boundary_periodic)>;
  ctrw::Get_position_periodic<Boundary const &> _getter_position;
  std::array<int, 4> _column_widths;
};

/** \class MeasurerTime_position_in_regions_periodic PTOF/Output.h
 * "PTOF/Output.h" \brief  Output time, tags, positions, and masses within
 * regions specified by masks. */
template <typename Subject, typename Geometry, typename Mask>
struct MeasurerTime_position_in_regions_periodic final
    : MeasurerTime<Subject, Geometry> {
  MeasurerTime_position_in_regions_periodic(
      Subject const &subject, Geometry const &geometry,
      Directories const &directories, std::string const &identifier,
      std::vector<std::reference_wrapper<const Mask>> masks,
      std::vector<double> thresholds = {}, int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,
                                        geometry,
                                        directories,
                                        "position_in_regions_periodic",
                                        identifier,
                                        precision},
        _getter_position{geometry.boundary_periodic}, _masks{masks},
        _thresholds{thresholds},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(1 + std::string{"Tag"}.length())),
            std::max(9 + precision, int(2 + std::string{"Position_"}.length())),
            std::max(9 + precision, int(1 + std::string{"Mass"}.length())),
            std::max(2, int(4 + std::string{"In_region_"}.length()))} {
    _thresholds.resize(masks.size(), 0.);
    _output << std::setw(_column_widths[0]) << "Time"
            << std::setw(_column_widths[1]) << "Tag";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[2])
              << "Position_" + std::to_string(dd);
    _output << std::setw(_column_widths[3]) << "Mass";
    for (std::size_t ii = 0; ii < _masks.size(); ++ii)
      _output << std::setw(_column_widths[4])
              << "In_region_" + std::to_string(ii);
    _output << std::setw(_column_widths[1]) << "Tag"
            << std::setw(_column_widths[2]) << "..."
            << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    for (auto const &part : _subject.particles()) {
      auto const &state = part.state_new();
      _output << std::setw(_column_widths[1]) << state.tag;
      if (!state.info.absorbed && part.state_old().time <= time &&
          !outside(state.cell)) {
        std::vector<int> in_region(_masks.size(), 0);
        for (std::size_t ii = 0; ii < _masks.size(); ++ii)
          if (_masks[ii].get()[state.cell] > _thresholds[ii])
            in_region[ii] = 1;
        if (std::any_of(in_region.begin(), in_region.end(),
                        [](int ii) { return ii > 0; })) {
          useful::print(_output, _getter_position(state), _column_widths[2]);
          _output << std::setw(_column_widths[3]) << state.mass;
          useful::print(_output, in_region, _column_widths[4]);
        }
      }
    }
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using Boundary =
      std::decay_t<decltype(std::declval<Geometry>().boundary_periodic)>;
  ctrw::Get_position_periodic<Boundary const &> _getter_position;
  std::vector<std::reference_wrapper<const Mask>> _masks;
  std::vector<double> _thresholds;
  std::array<int, 5> _column_widths;
};

/** \class MeasurerTime_position_mean_periodic PTOF/Output.h "PTOF/Output.h"
 *  \brief Output time and mean position (weighted by mass), with positions
 * accounting for periodicity. */
template <typename Subject, typename Geometry>
struct MeasurerTime_position_mean_periodic final
    : MeasurerTime<Subject, Geometry> {
  MeasurerTime_position_mean_periodic(Subject const &subject,
                                      Geometry const &geometry,
                                      Directories const &directories,
                                      std::string const &identifier,
                                      int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,     geometry,
                                        directories, "position_mean_periodic",
                                        identifier,  precision},
        getter_position{geometry.boundary_periodic},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(12, int(2 + std::string{"Position_mean_"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[1])
              << "Position_mean_" + std::to_string(dd);
    _output << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    useful::print(_output, position_mean(_subject, time, getter_position),
                  _column_widths[1]);
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using Boundary =
      std::decay_t<decltype(std::declval<Geometry>().boundary_periodic)>;
  ctrw::Get_position_periodic<Boundary const &> getter_position;
  std::array<int, 2> _column_widths;
};

/** \class MeasurerTime_position_second_moment_periodic PTOF/Output.h
 * "PTOF/Output.h" \brief  Output time and second moment of position (weighted
 * by mass), with position accounting for periodicity. */
template <typename Subject, typename Geometry>
struct MeasurerTime_position_second_moment_periodic final
    : MeasurerTime<Subject, Geometry> {
  MeasurerTime_position_second_moment_periodic(Subject const &subject,
                                               Geometry const &geometry,
                                               Directories const &directories,
                                               std::string const &identifier,
                                               int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,
                                        geometry,
                                        directories,
                                        "position_second_moment_periodic",
                                        identifier,
                                        precision},
        getter_position{geometry.boundary_periodic},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(
                9 + precision,
                int(2 + std::string{"Position_second_moment_"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[1])
              << "Position_second_moment_" + std::to_string(dd);
    _output << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    useful::print(_output,
                  position_second_moment(_subject, time, getter_position),
                  _column_widths[1]);
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using Boundary =
      std::decay_t<decltype(std::declval<Geometry>().boundary_periodic)>;
  ctrw::Get_position_periodic<Boundary const &> getter_position;
  std::array<int, 2> _column_widths;
};

/** \class MeasurerTime_position_variance_periodic PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output time and position variance (weighted by mass), with position
 * accounting for periodicity. */
template <typename Subject, typename Geometry>
struct MeasurerTime_position_variance_periodic final
    : MeasurerTime<Subject, Geometry> {
  MeasurerTime_position_variance_periodic(Subject const &subject,
                                          Geometry const &geometry,
                                          Directories const &directories,
                                          std::string const &identifier,
                                          int precision = 8)
      : MeasurerTime<Subject, Geometry>{subject,
                                        geometry,
                                        directories,
                                        "position_variance_periodic",
                                        identifier,
                                        precision},
        getter_position{geometry.boundary_periodic},
        _column_widths{
            std::max(9 + precision, int(1 + std::string{"Time"}.length())),
            std::max(9 + precision,
                     int(2 + std::string{"Position_variance_"}.length()))} {
    _output << std::setw(_column_widths[0]) << "Time";
    for (std::size_t dd = 0; dd < Geometry::dim; ++dd)
      _output << std::setw(_column_widths[1])
              << "Position_variance_" + std::to_string(dd);
    _output << "\n";
  }

  void operator()(double time) override {
    _output << std::setw(_column_widths[0]) << time;
    useful::print(_output, position_variance(_subject, time, getter_position),
                  _column_widths[1]);
    _output << "\n";
  }

private:
  using MeasurerTime<Subject, Geometry>::_output;
  using MeasurerTime<Subject, Geometry>::_subject;
  using Boundary =
      std::decay_t<decltype(std::declval<Geometry>().boundary_periodic)>;
  ctrw::Get_position_periodic<Boundary const &> getter_position;
  std::array<int, 2> _column_widths;
};

/** \class MeasurerTime_first_crossing_time PTOF/Output.h "PTOF/Output.h"
 *  \brief  Output first crossing times, tags, and masses. */
template <typename Subject, typename Geometry>
struct MeasurerTime_first_crossing_time final
    : MeasurerTime<Subject, Geometry> {
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
            << std::setw(_column_widths[2]) << "Mass"
            << "\n";
  }

  void operator()(double time) override {
    for (auto const &part : _subject.particles()) {
      auto const &state = part.state_new();
      auto const &state_old = part.state_old();
      if (!state.info.absorbed && part.state_old().time <= time) {
        auto pos = ctrw::Get_position_component_arg{dimension}(state);
        auto pos_old = ctrw::Get_position_component_arg{dimension}(state_old);
        if ((pos - position) * (pos_old - position) <= 0.) {
          Time crossing_time = ctrw::Get_time_interp{
              position,
              ctrw::Get_position_component_arg{dimension}}(state, state_old);
          Mass crossing_mass =
              state_old.mass + (state.mass - state_old.mass) /
                                   (state.time - state_old.time) *
                                   (time - state_old.time);
          _crossing_times_masses.insert(std::pair<Tag, std::pair<Time, Mass>>{
              state.tag, {crossing_time, crossing_mass}});
        }
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

} // namespace ptof

#endif /* PTOF_MEASURERTIME_H */
