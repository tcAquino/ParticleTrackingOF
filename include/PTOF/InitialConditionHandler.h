/**
 * @file   InitialConditionHandler.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Tue Feb  4 00:00:00 2025
 *
 * @brief Initial condition helpers.
 */

#ifndef PTOF_INITIALCONDITIONHANDLER_H
#define PTOF_INITIALCONDITIONHANDLER_H

#include "General/Meta.h"
#include "PTOF/InitialCondition_Cases.h"

namespace ptof {
/** @brief Generic initial condition handler. */
struct InitialConditionHandler_Generic {
  InitialConditionHandler_Generic() = delete;

  using Parameters = InitialConditionParameters_Cases;

  template <typename Particle, typename Geometry, typename VelocityField,
            typename SolverParameters, typename Mask = meta::Empty>
  static auto makeInitialCondition(Geometry const &geometry,
                                   VelocityField const &velocity_field,
                                   Parameters const &params_initial_condition,
                                   SolverParameters const &params_solvers,
                                   Mask &&mask = {}, double threshold = 0.) {
    return InitialCondition_Cases{meta::Selector_t<Particle>{},
                                  geometry,
                                  velocity_field,
                                  params_solvers.nr_particles,
                                  params_initial_condition,
                                  std::forward<Mask>(mask),
                                  threshold};
  }
};
} // namespace ptof

#endif /* PTOF_INITIALCONDITIONHANDLER_H */
