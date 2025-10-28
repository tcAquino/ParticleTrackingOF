/**
   \file PTOF/Solvers.h
   \author Tomas Aquino
   \date 19/01/2025
   \brief Solver parameters and utilities.
*/

#ifndef PTOF_SOLVERS_H
#define PTOF_SOLVERS_H

#include "General/IO.h"
#include "General/Meta.h"
#include "PTOF/SolverParameters.h"
#include "PTOF/Steppers.h"
#include <ostream>
#include <type_traits>
#include <utility>

namespace ptof {
template <typename Stepper_Advection_t, typename Stepper_Diffusion_t,
          typename Stepper_CTRW_t>
struct Solvers_Generic {
  Solvers_Generic() = delete;

  using Stepper_Advection = Stepper_Advection_t;
  using Stepper_Diffusion = Stepper_Diffusion_t;
  using Stepper_CTRW = Stepper_CTRW_t;
  using Parameters = SolverParameters_Generic<Stepper_Advection,
                                              Stepper_Diffusion, Stepper_CTRW>;
  static_assert(Parameters::advection || Parameters::diffusion,
                "No advection or diffusion must be present");
  static_assert(!(Parameters::diffusion &&
                  !std::is_same_v<Stepper_Diffusion, Steppers::Euler>),
                "Currently only no diffusion or stochastic forward Euler "
                "stepping are supported for diffusion");
  static_assert(
      !(Parameters::advection &&
        !std::is_same_v<Stepper_Advection, Steppers::Euler> &&
        !std::is_same_v<Stepper_Advection, Steppers::RK2> &&
        !std::is_same_v<Stepper_Advection, Steppers::RK4> &&
        !std::is_same_v<Stepper_Advection, Steppers::Heun>),
      "Currently only no advection, forward Euler, RK2, RK4, or Heun stepping "
      "are supported for advection");
  using Steppers_Transport = std::conditional_t<
      std::is_same_v<Stepper_Advection, Steppers::Heun>,
      Steppers_Advection_Heun_Diffusion_Euler,
      std::conditional_t<
          std::is_same_v<Stepper_Advection, Steppers::RK2>,
          Steppers_Advection_RK2_Diffusion_Euler,
          std::conditional_t<std::is_same_v<Stepper_Advection, Steppers::RK4>,
                             Steppers_Advection_RK4_Diffusion_Euler,
                             Steppers_Advection_Euler_Diffusion_Euler>>>;

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Solvers\n" << io::line();
    if constexpr (Parameters::advection) {
      output << "Advection: " << Steppers_Transport::advection_method << "\n";
    }
    if constexpr (Parameters::diffusion) {
      output << "Diffusion: " << Steppers_Transport::diffusion_method << "\n";
    }
    output << "CTRW: " << Stepper_CTRW::method << "\n";
    output << io::line();
    return output;
  }

  template <typename CTRW, typename Transitions, typename Time,
            typename Update = meta::DoNothing>
  static void evolve(CTRW &ctrw, Transitions &&transitions,
                     Parameters const &parameters_solvers, Time new_time,
                     Time old_time, Update &&update = {}) {
    return Stepper_CTRW::evolve(ctrw, transitions, parameters_solvers, new_time,
                                old_time, update);
  }

  template <typename Geometry, typename VelocityField, typename Boundary,
            typename TransportParameters, typename ReactionParameters,
            typename SolverParameters>
  static auto makeTimeGenerator(Geometry const &geometry,
                                VelocityField const &velocity_field,
                                Boundary &boundary,
                                TransportParameters const &params_transport,
                                ReactionParameters const &params_reaction,
                                SolverParameters const &params_solvers) {
    return Steppers_Transport::makeTimeGenerator(
        geometry, velocity_field, boundary, params_transport, params_reaction,
        params_solvers);
  }

  template <typename ParallelOption, typename VelocityField, typename Boundary,
            typename TransportParameters>
  static auto
  makeJumpGenerator(VelocityField &&velocity_field, Boundary &&boundary,
                    TransportParameters const &params_transport,
                    Parameters const &params_solvers, std::size_t dim) {
    if constexpr (Parameters::advection && Parameters::diffusion) {
      return Steppers_Transport::template makeJumpGenerator<ParallelOption>(
          std::forward<VelocityField>(velocity_field),
          std::forward<Boundary>(boundary), params_transport, params_solvers,
          dim);
    } else if constexpr (Parameters::advection) {
      return Steppers_Transport::template makeJumpGenerator_Advection<
          ParallelOption>(std::forward<VelocityField>(velocity_field),
                          std::forward<Boundary>(boundary), params_transport,
                          params_solvers, dim);
    } else if constexpr (Parameters::diffusion) {
      return Steppers_Transport::template makeJumpGenerator_Diffusion<
          ParallelOption>(std::forward<Boundary>(boundary), params_transport,
                          params_solvers, dim);
    }
  }
};
} // namespace ptof

#endif /* PTOF_SOLVERS_H */
