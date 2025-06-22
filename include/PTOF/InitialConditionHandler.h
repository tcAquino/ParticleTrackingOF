/**
   \file PTOF/InitialConditionHandler.h
   \author Tomas Aquino
   \date 04/02/2025
   \brief Initial condition helpers.
*/

#ifndef PTOF_INITIALCONDITIONHANDLER_H
#define PTOF_INITIALCONDITIONHANDLER_H

#include "PTOF/InitialCondition_Cases.h"
#include "General/Meta.h"

namespace ptof {
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
