/**
   \file PTOF/OutputHandler.h
   \author Tomas Aquino
   \date 04/02/2025
   \brief Output helpers.
*/

#ifndef PTOF_OUTPUTHANDLER_H
#define PTOF_OUTPUTHANDLER_H

#include "General/Meta.h"
#include "PTOF/Output_Cases.h"

namespace ptof {
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
