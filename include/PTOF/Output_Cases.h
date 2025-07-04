/**
   \file PTOF/Output_Cases.h
   \author Tomas Aquino
   \date 07/03/2022
   \brief Handle different types of output measurements.
*/

#ifndef PTOF_OUTPUT_CASES_H
#define PTOF_OUTPUT_CASES_H

#include "CTRW/Meta.h"
#include "General/IO.h"
#include "General/Meta.h"
#include "PTOF/BoundaryInfo.h"
#include "PTOF/Criterion.h"
#include "PTOF/Directories.h"
#include "PTOF/EndCriterionList.h"
#include "PTOF/MeasurementList.h"
#include "PTOF/MeasurementSpacingList.h"
#include "PTOF/Measurer.h"
#include "PTOF/MeasurerTime.h"
#include "PTOF/NextMeasurement.h"
#include "PTOF/OutputParameters_Cases.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <fvcGrad.H>
#include <initializer_list>
#include <limits>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace ptof {
/**
   \class Output_Cases PTOF/Output_Cases.h "PTOF/Output_Cases.h"
   \brief Output object to handle implemented output options.
*/
template <typename Subject, typename Geometry, typename Parameters>
class Output_Cases {
public:
  Parameters parameters; /**< Output parameters .*/

  /**
     \brief Constructor.
     \param subject CTRW object to measure.
     \param velocity_field Velocity field as a function of state.
     \param geometry Domain geometry info and utilities.
     \param boundary Boundary condition enforcer.
     \param directories Current case directory information.
     \param params_output Output parameters.
     \param identifier String to include in names of output files.
     \param masks Scalar fields reference wrappers.
     \param thresholds Thresholds for each mask, such that cells where a mask is
     above or equal to the threshold are considered.
  */
  template <typename Boundary, typename VelocityField = meta::Empty,
            typename Mask = meta::Empty>
  Output_Cases(Subject const &subject, VelocityField &&velocity_field,
               Geometry const &geometry, Boundary &boundary,
               Directories const &directories, Parameters &&params_output,
               std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks = {},
               std::vector<double> thresholds = {})
      : parameters{std::forward<Parameters>(params_output)} {
    using BoundaryInfo =
        BoundaryInfo_IfPresent_boundary_face_contact_point_reinjections<State>;
    boundary.add_boundary_info(std::unique_ptr<BoundaryInfo_Base<State>>(
        std::make_unique<BoundaryInfo>()));
    set_measurement_types(subject, geometry, boundary, directories, identifier,
                          velocity_field, masks, thresholds);
    set_end_criterion(subject);
    set_next_measurement_time();
  }

  /**
     \brief Constructor.
     \details Overload for initializer list arguments.
  */
  template <typename Boundary, typename VelocityField = meta::Empty,
            typename Mask>
  Output_Cases(Subject const &subject, VelocityField &&velocity_field,
               Geometry const &geometry, Boundary &boundary,
               Directories const &directories, Parameters &&params_output,
               std::string const &identifier,
               std::initializer_list<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {})
      : Output_Cases(subject, velocity_field, geometry, boundary, directories,
                     std::forward<Parameters>(params_output), identifier,
                     std::vector<std::reference_wrapper<const Mask>>{masks},
                     std::vector<double>{thresholds}) {}

  /**
     \brief Constructor.
     \details Overload for initializer list arguments.
  */
  template <typename Boundary, typename VelocityField = meta::Empty,
            typename Mask>
  Output_Cases(Subject const &subject, VelocityField &&velocity_field,
               Geometry const &geometry, Boundary &boundary,
               Directories const &directories, Parameters &&params_output,
               std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {})
      : Output_Cases(subject, velocity_field, geometry, boundary, directories,
                     std::forward<Parameters>(params_output), identifier, masks,
                     std::vector<double>{thresholds}) {}

  /**
     \brief Constructor.
     \details Overload for initializer list arguments.
  */
  template <typename Boundary, typename VelocityField = meta::Empty,
            typename Mask>
  Output_Cases(Subject const &subject, VelocityField &&velocity_field,
               Geometry const &geometry, Boundary &boundary,
               Directories const &directories, Parameters &&params_output,
               std::string const &identifier,
               std::initializer_list<std::reference_wrapper<const Mask>> masks,
               std::vector<double> thresholds = {})
      : Output_Cases(subject, velocity_field, geometry, boundary, directories,
                     std::forward<Parameters>(params_output), identifier, masks,
                     std::vector<double>{thresholds}) {}

  /**
     \return \c true  if end simulation criterion is satisfied, \c false
     otherwise.
  */
  bool done(double time) const { return _end_criterion->operator()(time); }

  /** \brief Set time of next measurement. */
  void set_next_measurement_time() {
    switch (MeasurementSpacingList::type(parameters.measurement_spacing)) {
    case MeasurementSpacingList::Type::linear_step: {
      _next_measurement = std::make_unique<
          NextMeasurementTime_linear_step<std::remove_reference_t<Parameters>>>(
          parameters);
      break;
    }
    case MeasurementSpacingList::Type::log_step: {
      _next_measurement = std::make_unique<
          NextMeasurementTime_log_step<std::remove_reference_t<Parameters>>>(
          parameters);
      break;
    }
    case MeasurementSpacingList::Type::linear: {
      _next_measurement = std::make_unique<
          NextMeasurementTime_linear<std::remove_reference_t<Parameters>>>(
          parameters);
      break;
    }
    case MeasurementSpacingList::Type::log: {
      _next_measurement = std::make_unique<
          NextMeasurementTime_log<std::remove_reference_t<Parameters>>>(
          parameters);
      break;
    }
    default: {
      throw std::runtime_error{std::string{"Measurement spacing "} +
                               parameters.measurement_spacing + " : " +
                               "Not supported"};
    }
    }
  }

  /** \return Time of next measurement. */
  double next_measurement_time() const { return _next_measurement->time(); }

  /**
     \brief Output requested measurements at given time and advance to next
     measurement.
  */
  void operator()(double time) {
    for (auto const &output : _output_time) {
      (*output)(time);
    }
    _next_measurement->advance();
  }

  /** \brief Output current information. */
  void operator()() {
    for (auto const &output : _output) {
      output->print();
    }
    for (auto const &output : _output_time) {
      output->print();
    }
  }

  /** \brief Update internal state. */
  void update(double time, double time_of_change) {
    for (auto const &output : _output_time) {
      output->update(time, time_of_change);
    }
  }

  /** \brief Output information about current object. */
  std::ostream &info_runtime(std::ostream &output) const {
    output << io::line() << "Output\n"

           << io::line();
    output << "Measurement time units: " << parameters.time_units << "\n"
           << "Measurement spacing: " << parameters.measurement_spacing << "\n"
           << "Minimum measurement time: " << parameters.time_min << "\n"
           << "Maximum measurement time: ";
    bool time_end_criterion = false;
    for (auto const &end_criterion : parameters.end_criteria) {
      auto criterion_type = EndCriterionList::type(end_criterion);
      if (criterion_type == EndCriterionList::Type::time) {
        time_end_criterion = true;
        break;
      }
    }
    if (!time_end_criterion) {
      output << parameters.time_max << "\n";
    } else {
      double time_end = parameters.end_values.at(EndCriterionList::Type::time);
      if (parameters.time_max < time_end) {
        throw std::runtime_error{
            "Last measurement time is lower than specified end time"};
      }
      output << time_end << "\n";
    }

    if (parameters.end_criteria.size() > 1) {
      output << "Number of end criteria to combine: "
             << parameters.end_criteria.size() << "\n";
      output << "End criterion logical combination: "
             << parameters.end_criterion_logical_combination << "\n";
    }

    for (auto const &end_criterion : parameters.end_criteria) {
      output << "End criterion: " << end_criterion << "\n";
      auto criterion_type = EndCriterionList::type(end_criterion);
      if (parameters.end_values.count(criterion_type)) {
        output << "  - End value: " << parameters.end_values.at(criterion_type)
               << "\n";
      }
    }
    output << "Measurements:\n";
    info_runtime_measurements(output);
    output << io::line();
    return output;
  }

private:
  /** \brief Output information about measurement types. */
  std::ostream &info_runtime_measurements(std::ostream &output) const {
    if (parameters.measurements.empty()) {
      return output << "  None\n";
    }
    for (auto const &measurement : parameters.measurements) {
      measurement->info_runtime(output);
    }
    return output;
  }

  /**
     \brief Set up output streams for requested output types.
     \param subject CTRW object to measure.
     \param directories Current case directory information.
     \param geometry Domain geometry info and utilities.
     \param boundary Boundary condition enforcer.
     \param identifier String to include in names of output files.
     \param velocity_field Velocity field as a function of state.
     \param masks Scalar field reference wrappers.
     \param thresholds Thresholds for each mask, such that cells where a mask
     is above or equal to the threshold are considered.
  */
  template <typename Boundary, typename VelocityField = meta::Empty,
            typename Mask = meta::Empty>
  void set_measurement_types(
      Subject const &subject, Geometry const &geometry, Boundary &boundary,
      Directories const &directories, std::string const &identifier,
      VelocityField &&velocity_field = {},
      std::vector<std::reference_wrapper<const Mask>> masks = {},
      std::vector<double> thresholds = {}) {
    for (auto const &measurement : parameters.measurements) {
      std::string for_measurement_type =
          std::string{"Measurement type "} + measurement->name + " : ";
      switch (MeasurementList::type(measurement->name)) {
      case MeasurementList::Type::position: {
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_position<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::position_in_regions: {
        if constexpr (!std::is_same_v<Mask, meta::Empty>) {
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_in_regions<Subject, Geometry, Mask>>(
                  subject, geometry, directories, identifier, masks, thresholds,
                  measurement->precision));
        } else {
          throw std::runtime_error{for_measurement_type +
                                   "Region masks not provided"};
        }
        break;
      }
      case MeasurementList::Type::position_mean: {
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_position_mean<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::position_second_moment: {
        _output_time.emplace_back(
            std::make_unique<
                MeasurerTime_position_second_moment<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::position_nth_moment: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Order;
        Measurement &measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<
                MeasurerTime_position_nth_moment<Subject, Geometry>>(
                subject, measurement_derived.order, geometry, directories,
                identifier, measurement->precision));
        break;
      }
      case MeasurementList::Type::position_moment: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Orders;
        Measurement &measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_position_moment<Subject, Geometry>>(
                subject, measurement_derived.orders, geometry, directories,
                identifier, measurement->precision));
        break;
      }
      case MeasurementList::Type::position_variance: {
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_position_variance<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::mass: {
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_mass<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::mass_absorbed: {
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_mass_absorbed<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::mass_adsorbed: {
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_mass_adsorbed<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::mass_in_regions: {
        if constexpr (!std::is_same_v<Mask, meta::Empty>) {
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_mass_in_regions<Subject, Geometry, Mask>>(
                  subject, geometry, directories, identifier, masks, thresholds,
                  measurement->precision));
        } else {
          throw std::runtime_error{for_measurement_type +
                                   "Region masks not provided"};
        }
        break;
      }
      case MeasurementList::Type::velocity: {
        if constexpr (!std::is_same_v<useful::remove_cvref_t<VelocityField>,
                                      meta::Empty>) {
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_vector_field<
                  Subject, Geometry, VelocityField const &>>(
                  subject, std::forward<VelocityField>(velocity_field),
                  geometry, directories, identifier, "U",
                  measurement->precision));
        } else {
          throw std::runtime_error{for_measurement_type +
                                   "Velocity field not provided"};
        }
        break;
      }
      case MeasurementList::Type::velocity_mean: {
        if constexpr (!std::is_same_v<useful::remove_cvref_t<VelocityField>,
                                      meta::Empty>) {
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_vector_field_mean<
                  Subject, Geometry, VelocityField const &>>(
                  subject, std::forward<VelocityField>(velocity_field),
                  geometry, directories, identifier, "U",
                  measurement->precision));
        } else {
          throw std::runtime_error{for_measurement_type +
                                   "Velocity field not provided"};
        }
        break;
      }
      case MeasurementList::Type::velocity_gradient: {
        if constexpr (!std::is_same_v<useful::remove_cvref_t<VelocityField>,
                                      meta::Empty>) {
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_tensor_field<Subject, Geometry>>(
                  subject, Foam::fvc::grad(velocity_field.field()), geometry,
                  directories, identifier, "gradU", measurement->precision));
        } else {
          throw std::runtime_error{for_measurement_type +
                                   "Velocity field not provided"};
        }
        break;
      }
      case MeasurementList::Type::velocity_gradient_mean: {
        if constexpr (!std::is_same_v<useful::remove_cvref_t<VelocityField>,
                                      meta::Empty>) {
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_tensor_field_mean<Subject, Geometry>>(
                  subject, Foam::fvc::grad(velocity_field.field()), geometry,
                  directories, identifier, "gradU", measurement->precision));
        } else {
          throw std::runtime_error{for_measurement_type +
                                   "Velocity field not provided"};
        }
        break;
      }
      case MeasurementList::Type::scalar_field: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement &measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_scalar_field<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived.field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::vector_field: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement &measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_vector_field<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived.field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::tensor_field: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement &measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_tensor_field<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived.field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::scalar_field_mean: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement &measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_scalar_field_mean<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived.field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::vector_field_mean: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement &measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_vector_field_mean<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived.field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::tensor_field_mean: {
        using Measurement =
            typename std::remove_reference_t<Parameters>::Measurement_Field;
        Measurement &measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_tensor_field_mean<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived.field_name, measurement->precision));
        break;
      }
      case MeasurementList::Type::position_periodic: {
        if constexpr (meta::has_periodicity_v<State>) {
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_periodic<Subject, Geometry>>(
                  subject, geometry, directories, identifier,
                  measurement->precision));
        } else {
          throw std::runtime_error{
              for_measurement_type +
              "Particle state does not define periodicity"};
        }
        break;
      }
      case MeasurementList::Type::position_in_regions_periodic: {
        if constexpr (meta::has_periodicity_v<State> &&
                      !std::is_same_v<Mask, meta::Empty>) {
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_position_in_regions_periodic<
                  Subject, Geometry, Mask>>(subject, geometry, directories,
                                            identifier, masks, thresholds,
                                            measurement->precision));
        } else {
          if constexpr (!meta::has_periodicity_v<State>) {
            throw std::runtime_error{
                for_measurement_type +
                "Particle state does not define periodicity"};
          }
          throw std::runtime_error{for_measurement_type +
                                   "Region masks not provided"};
        }
        break;
      }
      case MeasurementList::Type::position_mean_periodic: {
        if constexpr (meta::has_periodicity_v<State>) {
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_mean_periodic<Subject, Geometry>>(
                  subject, geometry, directories, identifier,
                  measurement->precision));
        } else {
          throw std::runtime_error{
              for_measurement_type +
              "Particle state does not define periodicity"};
        }
        break;
      }
      case MeasurementList::Type::position_second_moment_periodic: {
        if constexpr (meta::has_periodicity_v<State>) {
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_position_second_moment_periodic<
                  Subject, Geometry>>(subject, geometry, directories,
                                      identifier, measurement->precision));
        } else {
          throw std::runtime_error{
              for_measurement_type +
              "Particle state does not define periodicity"};
        }
        break;
      }
      case MeasurementList::Type::position_variance_periodic: {
        if constexpr (meta::has_periodicity_v<State>) {
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_variance_periodic<Subject, Geometry>>(
                  subject, geometry, directories, identifier,
                  measurement->precision));
        } else {
          throw std::runtime_error{
              for_measurement_type +
              "Particle state does not define periodicity"};
        }
        break;
      }
      case MeasurementList::Type::position_nth_moment_periodic: {
        if constexpr (meta::has_periodicity_v<State>) {
          using Measurement =
              typename std::remove_reference_t<Parameters>::Measurement_Order;
          Measurement &measurement_derived = cast<Measurement>(measurement);
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_nth_moment_periodic<Subject, Geometry>>(
                  subject, measurement_derived.order, geometry, directories,
                  identifier, measurement->precision));
        } else {
          throw std::runtime_error{
              for_measurement_type +
              "Particle state does not define periodicity"};
        }
        break;
      }
      case MeasurementList::Type::position_moment_periodic: {
        if constexpr (meta::has_periodicity_v<State>) {
          using Measurement =
              typename std::remove_reference_t<Parameters>::Measurement_Orders;
          Measurement &measurement_derived = cast<Measurement>(measurement);
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_moment_periodic<Subject, Geometry>>(
                  subject, measurement_derived.orders, geometry, directories,
                  identifier, measurement->precision));
        } else {
          throw std::runtime_error{
              for_measurement_type +
              "Particle state does not define periodicity"};
        }
        break;
      }
      case MeasurementList::Type::first_crossing_time: {
        using Measurement = typename std::remove_reference_t<
            Parameters>::Measurement_Dim_Position;
        Measurement &measurement_derived = cast<Measurement>(measurement);
        _output_time.emplace_back(
            std::make_unique<
                MeasurerTime_first_crossing_time<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement_derived.dim, measurement_derived.position,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::position_adsorbed: {
        _output_time.emplace_back(
            std::make_unique<MeasurerTime_position_adsorbed<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::position_adsorbed_periodic: {
        if constexpr (meta::has_periodicity_v<State>) {
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_position_adsorbed<Subject, Geometry>>(
                  subject, geometry, directories, identifier,
                  measurement->precision));
        } else {
          throw std::runtime_error{
              for_measurement_type +
              "Particle state does not define periodicity"};
        }
        break;
      }
      case MeasurementList::Type::mass_adsorbed_face: {
        _output_time.emplace_back(
            std::make_unique<
                MeasurerTime_mass_adsorbed_face<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::mass_adsorbed_face_periodic: {
        if constexpr (meta::has_periodicity_v<State>) {
          _output_time.emplace_back(
              std::make_unique<
                  MeasurerTime_mass_adsorbed_face_periodic<Subject, Geometry>>(
                  subject, geometry, directories, identifier,
                  measurement->precision));
        } else {
          throw std::runtime_error{
              for_measurement_type +
              "Particle state does not define periodicity"};
        }
        break;
      }
      case MeasurementList::Type::surface_reacted_mass: {
        using BoundaryInfo = BoundaryInfo_Record_surface_reacted_mass<State>;
        boundary.add_boundary_info(std::unique_ptr<BoundaryInfo_Base<State>>(
            std::make_unique<BoundaryInfo>()));
        _output_time.emplace_back(
            std::make_unique<
                MeasurerTime_surface_reacted_mass<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                dynamic_cast<BoundaryInfo const &>(
                    boundary.boundary_info_back()),
                measurement->precision));
        break;
      }
      case MeasurementList::Type::surface_reacted_mass_periodic: {
        if constexpr (meta::has_periodicity_v<State>) {
          using BoundaryInfo =
              BoundaryInfo_Record_surface_reacted_mass_periodic<State>;
          boundary.add_boundary_info(std::unique_ptr<BoundaryInfo_Base<State>>(
              std::make_unique<BoundaryInfo>()));
          _output_time.emplace_back(
              std::make_unique<MeasurerTime_surface_reacted_mass_periodic<
                  Subject, Geometry>>(subject, geometry, directories,
                                      identifier,
                                      dynamic_cast<BoundaryInfo const &>(
                                          boundary.boundary_info_back()),
                                      measurement->precision));
        } else {
          throw std::runtime_error{
              for_measurement_type +
              "Particle state does not define periodicity"};
        }
        break;
      }
      case MeasurementList::Type::absorption_time: {
        _output.emplace_back(
            std::make_unique<Measurer_absorption_time<Subject, Geometry>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::absorption_time_patch: {
        _output.emplace_back(
            std::make_unique<Measurer_absorption_time<Subject, Geometry, true>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::absorption_time_position: {
        _output.emplace_back(
            std::make_unique<
                Measurer_absorption_time<Subject, Geometry, false, true>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::absorption_time_patch_position: {
        _output.emplace_back(
            std::make_unique<
                Measurer_absorption_time<Subject, Geometry, true, true>>(
                subject, geometry, directories, identifier,
                measurement->precision));
        break;
      }
      case MeasurementList::Type::absorption_time_position_periodic: {
        if constexpr (meta::has_periodicity_v<State>) {
          _output.emplace_back(
              std::make_unique<Measurer_absorption_time<Subject, Geometry,
                                                        false, true, true>>(
                  subject, geometry, directories, identifier,
                  measurement->precision));
        } else {
          throw std::runtime_error{
              for_measurement_type +
              "Particle state does not define periodicity"};
        }
        break;
      }
      case MeasurementList::Type::absorption_time_patch_position_periodic: {
        if constexpr (meta::has_periodicity_v<State>) {
          _output.emplace_back(
              std::make_unique<Measurer_absorption_time<Subject, Geometry, true,
                                                        true, true>>(
                  subject, geometry, directories, identifier,
                  measurement->precision));
        } else {
          throw std::runtime_error{
              for_measurement_type +
              "Particle state does not define periodicity"};
        }
        break;
      }
      default: {
        throw std::runtime_error{for_measurement_type + "Not supported"};
      }
      }
    }
  }

  /** \brief Set end criterion. */
  void set_end_criterion(Subject const &subject) {
    // Note: The non-standard control flow structure of the if statements is
    // needed to deal with the quirks of if constexpr

    std::vector<std::unique_ptr<Criterion<Subject>>> end_criteria;
    for (auto const &end_criterion : parameters.end_criteria) {
      auto criterion_type = EndCriterionList::type(end_criterion);
      switch (criterion_type) {
      case EndCriterionList::Type::time: {
        end_criteria.emplace_back(std::make_unique<Criterion_time<Subject>>(
            subject, parameters.end_values.at(criterion_type)));
        break;
      }
      case EndCriterionList::Type::time_max: {
        if (parameters.time_max == std::numeric_limits<double>::infinity()) {
          throw std::runtime_error{
              std::string{"Infinite maximum time with end criterion "} +
              end_criterion};
        }
        end_criteria.emplace_back(std::make_unique<Criterion_time<Subject>>(
            subject, parameters.time_max));
        break;
      }
      case EndCriterionList::Type::mass_below: {
        end_criteria.emplace_back(
            std::make_unique<Criterion_mass_below<Subject>>(
                subject, parameters.end_values.at(criterion_type)));
        break;
      }
      case EndCriterionList::Type::mass_above: {
        end_criteria.emplace_back(
            std::make_unique<Criterion_mass_above<Subject>>(
                subject, parameters.end_values.at(criterion_type)));
        break;
      }
      case EndCriterionList::Type::all_absorbed: {
        end_criteria.emplace_back(
            std::make_unique<Criterion_all_absorbed<Subject>>(subject));
        break;
      }
      case EndCriterionList::Type::one_absorbed: {
        end_criteria.emplace_back(
            std::make_unique<Criterion_one_absorbed<Subject>>(subject));
        break;
      }
      case EndCriterionList::Type::fraction_not_absorbed: {
        end_criteria.emplace_back(
            std::make_unique<Criterion_fraction_not_absorbed<Subject>>(
                subject, parameters.end_values.at(criterion_type)));
        break;
      }
      default: {
        throw std::runtime_error{std::string{"End criterion "} + end_criterion +
                                 " : " + "Not supported"};
      }
      }
    }

    if (end_criteria.empty()) {
      throw std::runtime_error{"No end criterion"};
    }

    if (parameters.end_criterion_logical_combination == "single") {
      _end_criterion = std::move(end_criteria[0]);
    } else if (parameters.end_criterion_logical_combination == "and") {
      _end_criterion = std::make_unique<Criterion_and<Subject>>(
          subject, std::move(end_criteria));
    } else if (parameters.end_criterion_logical_combination == "or")
      _end_criterion = std::make_unique<Criterion_or<Subject>>(
          subject, std::move(end_criteria));
  }

  template <typename Measurement>
  Measurement &cast(
      std::unique_ptr<OutputParameters_Cases::Measurement> const &measurement) {
    return dynamic_cast<Measurement &>(*measurement.get());
  }

  using State = typename Subject::Particle::State;
  std::unique_ptr<Criterion<Subject>>
      _end_criterion; /**< To check if end criterion is met. */
  std::unique_ptr<NextMeasurementTime<std::remove_reference_t<Parameters>>>
      _next_measurement; /**< Handle next measurement time. */
  std::vector<std::unique_ptr<MeasurerTime<Subject, Geometry>>>
      _output_time; /**< Handle each output type given time. */
  std::vector<std::unique_ptr<Measurer<Subject, Geometry>>>
      _output; /**< Handle each output type given nothing. */
};
template <typename Subject, typename Geometry, typename Parameters,
          typename Boundary, typename VelocityField, typename Mask>
Output_Cases(
    Subject const &, VelocityField const &, Geometry const &, Boundary &,
    Directories const &, Parameters &&, std::string const &,
    std::vector<std::reference_wrapper<const Mask>>,
    std::vector<double>) -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename Boundary, typename VelocityField, typename Mask>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Boundary &, Directories const &, Parameters &&,
             std::string const &,
             std::vector<std::reference_wrapper<const Mask>>)
    -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename Boundary, typename VelocityField>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Boundary &, Directories const &, Parameters &&,
             std::string const &)
    -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename Boundary, typename VelocityField, typename Mask>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Boundary &, Directories const &, Parameters &&,
             std::string const &,
             std::initializer_list<std::reference_wrapper<const Mask>>,
             std::initializer_list<double>)
    -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename Boundary, typename VelocityField, typename Mask>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Boundary &, Directories const &, Parameters &&,
             std::string const &,
             std::initializer_list<std::reference_wrapper<const Mask>>)
    -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename Boundary, typename VelocityField, typename Mask>
Output_Cases(Subject const &, VelocityField const &, Geometry const &,
             Boundary &, Directories const &, Parameters &&,
             std::string const &,
             std::vector<std::reference_wrapper<const Mask>>,
             std::initializer_list<double>)
    -> Output_Cases<Subject, Geometry, Parameters>;
template <typename Subject, typename Geometry, typename Parameters,
          typename Boundary, typename VelocityField, typename Mask>
Output_Cases(
    Subject const &, VelocityField const &, Geometry const &, Boundary &,
    Directories const &, Parameters &&, std::string const &,
    std::initializer_list<std::reference_wrapper<const Mask>>,
    std::vector<double>) -> Output_Cases<Subject, Geometry, Parameters>;
} // namespace ptof

#endif /* PTOF_OUTPUT_CASES_H */
