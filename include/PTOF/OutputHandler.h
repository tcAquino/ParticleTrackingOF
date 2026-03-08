/**
 * @file   OutputHandler.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Tue Feb  4 00:00:00 2025
 *
 * @brief Output helpers.
 */

#ifndef PTOF_OUTPUTHANDLER_H
#define PTOF_OUTPUTHANDLER_H

#include "General/Meta.h"
#include "PTOF/Output_Cases.h"

namespace ptof {
/** @brief Class for handling output of different types. */
struct OutputHandler_Generic {
  OutputHandler_Generic() = delete;

  using Parameters = OutputParameters_Cases;

  template <typename Subject, typename VelocityField, typename Geometry,
            typename Boundary, typename Mask = meta::Empty>
  static auto
  makeOutput(Subject const &subject, VelocityField const &velocity_field,
             Geometry const &geometry, Boundary &boundary,
             Directories const &directories, Parameters const &params_output,
             std::string const &identifier,
             std::vector<std::reference_wrapper<const Mask>> masks = {},
             std::vector<double> thresholds = {}) {
    return Output_Cases{subject,    velocity_field, geometry,
                        boundary,   directories,    params_output,
                        identifier, masks,          thresholds};
  }

  template <typename Subject, typename VelocityField, typename Geometry,
            typename Boundary, typename Mask>
  static auto
  makeOutput(Subject const &subject, VelocityField const &velocity_field,
             Geometry const &geometry, Boundary &boundary,
             Directories const &directories, Parameters const &params_output,
             std::string const &identifier,
             std::initializer_list<std::reference_wrapper<const Mask>> masks,
             std::initializer_list<double> thresholds = {}) {
    return Output_Cases{subject,    velocity_field, geometry,
                        boundary,   directories,    params_output,
                        identifier, masks,          thresholds};
  }

  template <typename Subject, typename VelocityField, typename Geometry,
            typename Boundary, typename Mask>
  static auto
  makeOutput(Subject const &subject, VelocityField const &velocity_field,
             Geometry const &geometry, Boundary &boundary,
             Directories const &directories, Parameters const &params_output,
             std::string const &identifier,
             std::vector<std::reference_wrapper<const Mask>> masks,
             std::initializer_list<double> thresholds = {}) {
    return Output_Cases{subject,    velocity_field, geometry,
                        boundary,   directories,    params_output,
                        identifier, masks,          thresholds};
  }

  template <typename Subject, typename VelocityField, typename Geometry,
            typename Boundary, typename Mask>
  static auto
  makeOutput(Subject const &subject, VelocityField const &velocity_field,
             Geometry const &geometry, Boundary &boundary,
             Directories const &directories, Parameters const &params_output,
             std::string const &identifier,
             std::initializer_list<std::reference_wrapper<const Mask>> masks,
             std::vector<double> thresholds = {}) {
    return Output_Cases{subject,    velocity_field, geometry,
                        boundary,   directories,    params_output,
                        identifier, masks,          thresholds};
  }
};
} // namespace ptof

#endif /* PTOF_OUTPUTHANDLER_H */
