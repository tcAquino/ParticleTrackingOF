/**
   \file PTOF/TransportHandler.h
   \author Tomas Aquino
   \date 19/01/2025
   \brief Transport parameters and utilities.
*/

#ifndef PTOF_TRANSPORTHANDLER_H
#define PTOF_TRANSPORTHANDLER_H

#include "General/IO.h"
#include "General/Meta.h"
#include "PTOF/Advection.h"
#include "PTOF/Transitions.h"
#include "PTOF/TransportParameters.h"
#include <memory>
#include <type_traits>
#include <utility>

namespace ptof {
template <typename Solvers_t, typename Parameters_t>
struct TransportHandler_LinearInterp {
  TransportHandler_LinearInterp() = delete;

  using Solvers = Solvers_t;
  using Parameters = Parameters_t;

  template <typename VelocityField, typename Geometry, typename Boundary,
            typename ReactionParameters>
  static auto
  makeTransitions(VelocityField const &velocity_field, Geometry const &geometry,
                  Boundary &boundary, Parameters const &params_transport,
                  ReactionParameters const &params_reaction,
                  typename Solvers::Parameters const &params_solvers) {
    return makeTransportTransitions<Solvers>(velocity_field, geometry, boundary,
                                             params_transport, params_reaction,
                                             params_solvers);
  }

  template <typename Geometry, typename VelocityData>
  static auto makeVelocityInterpolator(Geometry const &geometry,
                                       VelocityData &&velocity_data) {
    if constexpr (Solvers::Parameters::advection) {
      return makeLinearVelocityInterpolator(
          geometry, std::forward<VelocityData>(velocity_data));
    } else {
      return meta::Empty{};
    }
  }

  template <typename Geometry, typename VelocityData>
  static auto makeVelocityTimeInterpolator(Geometry const &geometry,
                                           VelocityData &&velocity_data_new,
                                           VelocityData &&velocity_data_old,
                                           double time_new, double time_old) {
    if constexpr (Solvers::Parameters::advection) {
      using Type = decltype(makeVelocityInterpolator(
          geometry, std::forward<VelocityData>(velocity_data_new)));
      return ptof::Field_TimeInterpolation{
          std::unique_ptr<Type>{new Type{makeVelocityInterpolator(
              geometry, std::forward<VelocityData>(velocity_data_new))}},
          std::unique_ptr<Type>{new Type{makeVelocityInterpolator(
              geometry, std::forward<VelocityData>(velocity_data_old))}},
          time_new, time_old};
    } else {
      return meta::Empty{};
    }
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Transport\n" << io::line() << "Advection: ";
    if constexpr (Solvers::Parameters::advection) {
      output << "Yes\n";
    } else {
      output << "No\n";
    }
    output << "Diffusion: ";
    if constexpr (Solvers::Parameters::diffusion) {
      output << "Yes\n";
    } else {
      output << "No\n";
    }
    if constexpr (Solvers::Parameters::advection) {
      output << "Velocity interpolation: Linear\n";
    }
    output << io::line();
    return output;
  }
  };

template <typename Solvers>
using TransportHandler_Generic = TransportHandler_LinearInterp<
    Solvers,
    std::conditional_t<Solvers::Parameters::advection &&
                           Solvers::Parameters::diffusion,
                       TransportParameters_AdvectionDiffusion,
                       std::conditional_t<Solvers::Parameters::advection,
                                          TransportParameters_Advection,
                                          TransportParameters_Diffusion>>>;

template <typename Solvers>
using TransportHandler_Bcc = TransportHandler_LinearInterp<
    Solvers,
    std::conditional_t<Solvers::Parameters::advection &&
                           Solvers::Parameters::diffusion,
                       TransportParameters_AdvectionDiffusion_Bcc,
                       std::conditional_t<Solvers::Parameters::advection,
                                          TransportParameters_Advection_Bcc,
                                          TransportParameters_Diffusion_Bcc>>>;

} // namespace ptof

#endif /* PTOF_TRANSPORTHANDLER_H */
