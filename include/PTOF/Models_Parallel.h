/**
 \file PTOF/Models_Parallel.h
 \author Tomás Aquino
 \date 12/02/2024
 */

#ifndef PTOF_MODELS_PARALLEL_H
#define PTOF_MODELS_PARALLEL_H

#include "CTRW/CTRW_Parallel.h"
#include "General/Parallel.h"
#include "PTOF/Models.h"

namespace ptof {
/** \namespace ptof::model_advection_diffusion_2d_parallel
 Definitions for 2D advective--diffusive transport, in parallel. */
namespace model_advection_diffusion_2d_parallel {
struct Model {
  inline static const std::string name{"advection_diffusion_2d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry = Geometry<2, meta::ParallelOptions::Parallel>;
using model_advection_diffusion_2d::Info;
using model_advection_diffusion_2d::State;
using CTRW = ctrw::CTRW_Parallel<State>;
using model_advection_diffusion_2d::InitialCondition;
using model_advection_diffusion_2d::Output;
using model_advection_diffusion_2d::Reaction;
using model_advection_diffusion_2d::Solvers;
using model_advection_diffusion_2d::Transport;

template <typename Transport, typename Solvers, typename VelocityField,
          typename Geometry, typename Boundary, typename ReactionParameters,
          typename BulkReaction>
auto makeTransitions(VelocityField const &velocity_field,
                     Geometry const &geometry, Boundary &boundary,
                     typename Transport::Parameters const &params_transport,
                     ReactionParameters const &params_reaction,
                     typename Solvers::Parameters const &params_solvers,
                     BulkReaction const &bulk_reaction) {
  using Transitions = decltype(ctrw::Transitions_CTRW_Transport_Reaction{
      Transport::template makeTransitions<Solvers>(
          velocity_field, geometry, boundary, params_transport, params_reaction,
          params_solvers),
      bulk_reaction});
  std::vector<Transitions> transitions;
  for (std::size_t thread = 0;
       thread < useful::get_num_threads(meta::ParallelOptions::Parallel{});
       ++thread)
    transitions.emplace_back(Transport::template makeTransitions<Solvers>(
                                 velocity_field, geometry, boundary,
                                 params_transport, params_reaction,
                                 params_solvers),
                             bulk_reaction);

  return transitions;
}
} // namespace model_advection_diffusion_2d_parallel

/** \namespace ptof::model_advection_2d_parallel
 Definitions for 2D advective transport, in parallel. */
namespace model_advection_2d_parallel {
struct Model {
  inline static const std::string name{"advection_2d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_2d::Info;
using model_advection_2d::InitialCondition;
using model_advection_2d::Output;
using model_advection_2d::Reaction;
using model_advection_2d::Solvers;
using model_advection_2d::State;
using model_advection_2d::Transport;
using model_advection_diffusion_2d_parallel::CTRW;
using model_advection_diffusion_2d_parallel::Geometry;
using model_advection_diffusion_2d_parallel::makeTransitions;
}; // namespace model_advection_2d_parallel

/** \namespace ptof::model_advection_diffusion_fpt_2d_parallel
 Definitions for first-passage times under 2D advective--diffusive transport, in
 parallel. */
namespace model_advection_diffusion_fpt_2d_parallel {
struct Model {
  inline static const std::string name{"advection_diffusion_fpt_2d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry =
    Geometry<2, meta::ParallelOptions::Parallel, Dynamics::Type::firstpassage>;
using model_advection_diffusion_fpt_2d::Info;
using model_advection_diffusion_fpt_2d::State;
using CTRW = ctrw::CTRW_Parallel<State>;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_advection_diffusion_fpt_2d::InitialCondition;
using model_advection_diffusion_fpt_2d::Output;
using model_advection_diffusion_fpt_2d::Reaction;
using model_advection_diffusion_fpt_2d::Solvers;
using model_advection_diffusion_fpt_2d::Transport;
} // namespace model_advection_diffusion_fpt_2d_parallel

/** \namespace ptof::model_advection_diffusion_surface_decay_2d_parallel
 Definitions for 2D advective--diffusive transport with surface reaction \f$ A_F
 + B_S \to B_S\f$, in parallel. */
namespace model_advection_diffusion_surface_decay_2d_parallel {
struct Model {
  inline static const std::string name{
      "advection_diffusion_surface_decay_2d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d_parallel::CTRW;
using model_advection_diffusion_2d_parallel::Geometry;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_advection_diffusion_surface_decay_2d::Info;
using model_advection_diffusion_surface_decay_2d::InitialCondition;
using model_advection_diffusion_surface_decay_2d::Output;
using model_advection_diffusion_surface_decay_2d::Reaction;
using model_advection_diffusion_surface_decay_2d::Solvers;
using model_advection_diffusion_surface_decay_2d::State;
using model_advection_diffusion_surface_decay_2d::Transport;
} // namespace model_advection_diffusion_surface_decay_2d_parallel

/** \namespace ptof::model_periodic_cartesian_advection_diffusion_2d_parallel
 Definitions for 2D advective--diffusive transport with some periodic boundaries
 aligned with the Cartesian axes, in parallel. */
namespace model_periodic_cartesian_advection_diffusion_2d_parallel {
struct Model {
  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_2d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry =
    Geometry_Periodic_Cartesian<2, meta::ParallelOptions::Parallel>;
using model_periodic_cartesian_advection_diffusion_2d::Info;
using model_periodic_cartesian_advection_diffusion_2d::State;
using CTRW = ctrw::CTRW_Parallel<State>;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_periodic_cartesian_advection_diffusion_2d::InitialCondition;
using model_periodic_cartesian_advection_diffusion_2d::Output;
using model_periodic_cartesian_advection_diffusion_2d::Reaction;
using model_periodic_cartesian_advection_diffusion_2d::Solvers;
using model_periodic_cartesian_advection_diffusion_2d::Transport;
} // namespace model_periodic_cartesian_advection_diffusion_2d_parallel

/** \namespace ptof::model_periodic_cartesian_advection_2d_parallel
 Definitions for 2D advective transport with some periodic boundaries aligned
 with the Cartesian axes. */
namespace model_periodic_cartesian_advection_2d_parallel {
struct Model {
  inline static const std::string name{
      "periodic_cartesian_advection_2d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d_parallel::makeTransitions;
using model_periodic_cartesian_advection_2d::Info;
using model_periodic_cartesian_advection_2d::InitialCondition;
using model_periodic_cartesian_advection_2d::Output;
using model_periodic_cartesian_advection_2d::Reaction;
using model_periodic_cartesian_advection_2d::Solvers;
using model_periodic_cartesian_advection_2d::State;
using model_periodic_cartesian_advection_2d::Transport;
using model_periodic_cartesian_advection_diffusion_2d_parallel::CTRW;
using model_periodic_cartesian_advection_diffusion_2d_parallel::Geometry;
}; // namespace model_periodic_cartesian_advection_2d_parallel

/** \namespace
 ptof::model_periodic_cartesian_advection_diffusion_fpt_2d_parallel Definitions
 for first-passage times under 2D advective--diffusive transport with some
 periodic boundaries aligned with the Cartesian axes, in parallel. */
namespace model_periodic_cartesian_advection_diffusion_fpt_2d_parallel {
struct Model {
  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_fpt_2d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry = Geometry_Periodic_Cartesian<2, meta::ParallelOptions::Parallel,
                                             Dynamics::Type::firstpassage>;
using model_periodic_cartesian_advection_diffusion_fpt_2d::Info;
using model_periodic_cartesian_advection_diffusion_fpt_2d::State;
using CTRW = ctrw::CTRW_Parallel<State>;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_periodic_cartesian_advection_diffusion_fpt_2d::InitialCondition;
using model_periodic_cartesian_advection_diffusion_fpt_2d::Output;
using model_periodic_cartesian_advection_diffusion_fpt_2d::Reaction;
using model_periodic_cartesian_advection_diffusion_fpt_2d::Solvers;
using model_periodic_cartesian_advection_diffusion_fpt_2d::Transport;
} // namespace model_periodic_cartesian_advection_diffusion_fpt_2d_parallel

/** \namespace
 ptof::model_periodic_cartesian_advection_diffusion_surface_decay_2d_parallel
 Definitions for first-passage times under 2D advective--diffusive transport
 with surface reaction \f$ A_F + B_S \to B_S\f$, with some periodic boundaries
 aligned with the Cartesian axes, in parallel. */
namespace model_periodic_cartesian_advection_diffusion_surface_decay_2d_parallel {
struct Model {
  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_surface_decay_2d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d_parallel::makeTransitions;
using model_periodic_cartesian_advection_diffusion_2d_parallel::CTRW;
using model_periodic_cartesian_advection_diffusion_2d_parallel::Geometry;
using model_periodic_cartesian_advection_diffusion_surface_decay_2d::Info;
using model_periodic_cartesian_advection_diffusion_surface_decay_2d::
    InitialCondition;
using model_periodic_cartesian_advection_diffusion_surface_decay_2d::Output;
using model_periodic_cartesian_advection_diffusion_surface_decay_2d::Reaction;
using model_periodic_cartesian_advection_diffusion_surface_decay_2d::Solvers;
using model_periodic_cartesian_advection_diffusion_surface_decay_2d::State;
using model_periodic_cartesian_advection_diffusion_surface_decay_2d::Transport;
} // namespace
  // model_periodic_cartesian_advection_diffusion_surface_decay_2d_parallel

/** \namespace ptof::model_advection_diffusion_3d_parallel
 Definitions for 2D advective--diffusive transport, in parallel. */
namespace model_advection_diffusion_3d_parallel {
struct Model {
  inline static const std::string name{"advection_diffusion_3d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry = Geometry<3, meta::ParallelOptions::Parallel>;
using model_advection_diffusion_3d::Info;
using model_advection_diffusion_3d::State;
using CTRW = ctrw::CTRW_Parallel<State>;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_advection_diffusion_3d::InitialCondition;
using model_advection_diffusion_3d::Output;
using model_advection_diffusion_3d::Reaction;
using model_advection_diffusion_3d::Solvers;
using model_advection_diffusion_3d::Transport;
} // namespace model_advection_diffusion_3d_parallel

/** \namespace ptof::model_advection_3d_parallel
 Definitions for 3D advective--diffusive transport, in parallel. */
namespace model_advection_3d_parallel {
struct Model {
  inline static const std::string name{"advection_3d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_3d::Info;
using model_advection_3d::InitialCondition;
using model_advection_3d::Output;
using model_advection_3d::Reaction;
using model_advection_3d::Solvers;
using model_advection_3d::State;
using model_advection_3d::Transport;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_advection_diffusion_3d_parallel::CTRW;
using model_advection_diffusion_3d_parallel::Geometry;
}; // namespace model_advection_3d_parallel

/** \namespace ptof::model_advection_diffusion_fpt_3d_parallel
 Definitions for first-passage times under 3D advective--diffusive transport, in
 parallel. */
namespace model_advection_diffusion_fpt_3d_parallel {
struct Model {
  inline static const std::string name{"advection_diffusion_fpt_3d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry =
    Geometry<3, meta::ParallelOptions::Parallel, Dynamics::Type::firstpassage>;
using model_advection_diffusion_fpt_3d::Info;
using model_advection_diffusion_fpt_3d::State;
using CTRW = ctrw::CTRW_Parallel<State>;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_advection_diffusion_fpt_3d::InitialCondition;
using model_advection_diffusion_fpt_3d::Output;
using model_advection_diffusion_fpt_3d::Reaction;
using model_advection_diffusion_fpt_3d::Solvers;
using model_advection_diffusion_fpt_3d::Transport;
} // namespace model_advection_diffusion_fpt_3d_parallel

/** \namespace ptof::model_advection_diffusion_surface_decay_3d_parallel
 Definitions for 3D advective--diffusive transport with surface reaction \f$ A_F
 + B_S \to B_S\f$, in parallel. */
namespace model_advection_diffusion_surface_decay_3d_parallel {
struct Model {
  inline static const std::string name{
      "advection_diffusion_surface_decay_3d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d_parallel::makeTransitions;
using model_advection_diffusion_3d_parallel::CTRW;
using model_advection_diffusion_3d_parallel::Geometry;
using model_advection_diffusion_surface_decay_3d::Info;
using model_advection_diffusion_surface_decay_3d::InitialCondition;
using model_advection_diffusion_surface_decay_3d::Output;
using model_advection_diffusion_surface_decay_3d::Reaction;
using model_advection_diffusion_surface_decay_3d::Solvers;
using model_advection_diffusion_surface_decay_3d::State;
using model_advection_diffusion_surface_decay_3d::Transport;
} // namespace model_advection_diffusion_surface_decay_3d_parallel

/** \namespace ptof::model_periodic_cartesian_advection_diffusion_3d_parallel
 Definitions for 3D advective--diffusive transport with some periodic boundaries
 aligned with the Cartesian axes. */
namespace model_periodic_cartesian_advection_diffusion_3d_parallel {
struct Model {
  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_3d_parallels"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry =
    Geometry_Periodic_Cartesian<3, meta::ParallelOptions::Parallel>;
using model_periodic_cartesian_advection_diffusion_3d::Info;
using model_periodic_cartesian_advection_diffusion_3d::State;
using CTRW = ctrw::CTRW_Parallel<State>;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_periodic_cartesian_advection_diffusion_3d::InitialCondition;
using model_periodic_cartesian_advection_diffusion_3d::Output;
using model_periodic_cartesian_advection_diffusion_3d::Reaction;
using model_periodic_cartesian_advection_diffusion_3d::Solvers;
using model_periodic_cartesian_advection_diffusion_3d::Transport;
} // namespace model_periodic_cartesian_advection_diffusion_3d_parallel

/** \namespace ptof::model_periodic_cartesian_advection_3d_parallel
 Definitions for 3D advective transport with some periodic boundaries aligned
 with the Cartesian axes. */
namespace model_periodic_cartesian_advection_3d_parallel {
struct Model {
  inline static const std::string name{
      "periodic_cartesian_advection_3d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d_parallel::makeTransitions;
using model_periodic_cartesian_advection_3d::Info;
using model_periodic_cartesian_advection_3d::InitialCondition;
using model_periodic_cartesian_advection_3d::Output;
using model_periodic_cartesian_advection_3d::Reaction;
using model_periodic_cartesian_advection_3d::Solvers;
using model_periodic_cartesian_advection_3d::State;
using model_periodic_cartesian_advection_3d::Transport;
using model_periodic_cartesian_advection_diffusion_3d_parallel::CTRW;
using model_periodic_cartesian_advection_diffusion_3d_parallel::Geometry;
}; // namespace model_periodic_cartesian_advection_3d_parallel

/** \namespace
 ptof::model_periodic_cartesian_advection_diffusion_fpt_3d_parallel Definitions
 for first-passage times under 3D advective--diffusive transport with some
 periodic boundaries aligned with the Cartesian axes, in parallel. */
namespace model_periodic_cartesian_advection_diffusion_fpt_3d_parallel {
struct Model {
  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_fpt_3d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry = Geometry_Periodic_Cartesian<3, meta::ParallelOptions::Parallel,
                                             Dynamics::Type::firstpassage>;
using model_periodic_cartesian_advection_diffusion_fpt_3d::Info;
using model_periodic_cartesian_advection_diffusion_fpt_3d::State;
using CTRW = ctrw::CTRW_Parallel<State>;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_periodic_cartesian_advection_diffusion_fpt_3d::InitialCondition;
using model_periodic_cartesian_advection_diffusion_fpt_3d::Output;
using model_periodic_cartesian_advection_diffusion_fpt_3d::Reaction;
using model_periodic_cartesian_advection_diffusion_fpt_3d::Solvers;
using model_periodic_cartesian_advection_diffusion_fpt_3d::Transport;
} // namespace model_periodic_cartesian_advection_diffusion_fpt_3d_parallel

/** \namespace
 ptof::model_periodic_cartesian_advection_diffusion_surface_decay_3d_parallel
 Definitions for 3D advective--diffusive transport with surface reaction \f$ A_F
 + B_S \to B_S\f$, with some periodic boundaries aligned with the Cartesian
 axes, in parallel. */
namespace model_periodic_cartesian_advection_diffusion_surface_decay_3d_parallel {
struct Model {
  inline static const std::string name{
      "periodic_cartesian_advection_diffusion_surface_decay_3d_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d_parallel::makeTransitions;
using model_periodic_cartesian_advection_diffusion_3d_parallel::CTRW;
using model_periodic_cartesian_advection_diffusion_3d_parallel::Geometry;
using model_periodic_cartesian_advection_diffusion_surface_decay_3d::Info;
using model_periodic_cartesian_advection_diffusion_surface_decay_3d::
    InitialCondition;
using model_periodic_cartesian_advection_diffusion_surface_decay_3d::Output;
using model_periodic_cartesian_advection_diffusion_surface_decay_3d::Reaction;
using model_periodic_cartesian_advection_diffusion_surface_decay_3d::Solvers;
using model_periodic_cartesian_advection_diffusion_surface_decay_3d::State;
using model_periodic_cartesian_advection_diffusion_surface_decay_3d::Transport;
} // namespace
  // model_periodic_cartesian_advection_diffusion_surface_decay_3d_parallel

/** \namespace ptof::model_bcc_cartesian_advection_diffusion_parallel
 Definitions for 3D advective--diffusive transport in a body centered cubic
 pack, based on the primitive unit cel, in parallel. */
namespace model_bcc_cartesian_advection_diffusion_parallel {
struct Model {
  inline static const std::string name{
      "bcc_cartesian_advection_diffusion_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry = Geometry_Bcc<meta::ParallelOptions::Parallel>;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_bcc_cartesian_advection_diffusion::Info;
using model_bcc_cartesian_advection_diffusion::InitialCondition;
using model_bcc_cartesian_advection_diffusion::Output;
using model_bcc_cartesian_advection_diffusion::Reaction;
using model_bcc_cartesian_advection_diffusion::Solvers;
using model_bcc_cartesian_advection_diffusion::State;
using model_bcc_cartesian_advection_diffusion::Transport;
using model_periodic_cartesian_advection_diffusion_3d_parallel::CTRW;
} // namespace model_bcc_cartesian_advection_diffusion_parallel

/** \namespace ptof::model_bcc_cartesian_advection_diffusion_fpt_parallel
 Definitions for advective--diffusive transport in a body centered cubic
 beadpack, based on the primitive unit cell, in parallel. */
namespace model_bcc_cartesian_advection_diffusion_fpt_parallel {
struct Model {
  inline static const std::string name{
      "bcc_cartesian_advection_diffusion_fpt_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry =
    Geometry_Bcc<meta::ParallelOptions::Parallel, Periodicity::Type::cartesian,
                 Dynamics::Type::firstpassage>;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_bcc_cartesian_advection_diffusion_fpt::Info;
using model_bcc_cartesian_advection_diffusion_fpt::InitialCondition;
using model_bcc_cartesian_advection_diffusion_fpt::Output;
using model_bcc_cartesian_advection_diffusion_fpt::Reaction;
using model_bcc_cartesian_advection_diffusion_fpt::Solvers;
using model_bcc_cartesian_advection_diffusion_fpt::State;
using model_bcc_cartesian_advection_diffusion_fpt::Transport;
using model_periodic_cartesian_advection_diffusion_fpt_3d_parallel::CTRW;
} // namespace model_bcc_cartesian_advection_diffusion_fpt_parallel

/** \namespace ptof::model_bcc_cartesian_advection_parallel
 Definitions for first-passage times under advective--diffusive transport in a
 body centered cubic beadpack, based on the primitive unit cell. */
namespace model_bcc_cartesian_advection_parallel {
struct Model {
  inline static const std::string name{"bcc_cartesian_advectio_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d_parallel::makeTransitions;
using model_bcc_cartesian_advection::Info;
using model_bcc_cartesian_advection::InitialCondition;
using model_bcc_cartesian_advection::Output;
using model_bcc_cartesian_advection::Reaction;
using model_bcc_cartesian_advection::Solvers;
using model_bcc_cartesian_advection::State;
using model_bcc_cartesian_advection::Transport;
using model_bcc_cartesian_advection_diffusion_parallel::CTRW;
using model_bcc_cartesian_advection_diffusion_parallel::Geometry;
} // namespace model_bcc_cartesian_advection_parallel

/** \namespace
 ptof::model_bcc_cartesian_advection_diffusion_surface_decay_parallel
 Definitions for advective transport in a body centered cubic beadpack, based on
 the primitive unit cell. */
namespace model_bcc_cartesian_advection_diffusion_surface_decay_parallel {
struct Model {
  inline static const std::string name{
      "bcc_cartesian_advection_diffusion_surface_decay_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d_parallel::makeTransitions;
using model_bcc_cartesian_advection_diffusion_parallel::CTRW;
using model_bcc_cartesian_advection_diffusion_parallel::Geometry;
using model_bcc_cartesian_advection_diffusion_surface_decay::Info;
using model_bcc_cartesian_advection_diffusion_surface_decay::InitialCondition;
using model_bcc_cartesian_advection_diffusion_surface_decay::Output;
using model_bcc_cartesian_advection_diffusion_surface_decay::Reaction;
using model_bcc_cartesian_advection_diffusion_surface_decay::Solvers;
using model_bcc_cartesian_advection_diffusion_surface_decay::State;
using model_bcc_cartesian_advection_diffusion_surface_decay::Transport;
} // namespace model_bcc_cartesian_advection_diffusion_surface_decay_parallel

/** \namespace ptof::model_bcc_symmetryplanes_advection_diffusion_parallel
 Definitions for advective--diffusive transport with surface reaction \f$ A_F +
 B_S \to B_S\f$ in a body centered cubic beadpack, based on the primitive unit
 cell, in parallel. */
namespace model_bcc_symmetryplanes_advection_diffusion_parallel {
struct Model {
  inline static const std::string name{
      "bcc_symmetryplanes_advection_diffusion_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry = Geometry_Bcc<meta::ParallelOptions::Parallel,
                              Periodicity::Type::symmetryplanes>;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_bcc_cartesian_advection_diffusion_parallel::CTRW;
using model_bcc_symmetryplanes_advection_diffusion::Info;
using model_bcc_symmetryplanes_advection_diffusion::InitialCondition;
using model_bcc_symmetryplanes_advection_diffusion::Output;
using model_bcc_symmetryplanes_advection_diffusion::Reaction;
using model_bcc_symmetryplanes_advection_diffusion::Solvers;
using model_bcc_symmetryplanes_advection_diffusion::State;
using model_bcc_symmetryplanes_advection_diffusion::Transport;
} // namespace model_bcc_symmetryplanes_advection_diffusion_parallel

/** \namespace ptof::model_bcc_symmetryplanes_advection_diffusion_fpt_parallel
 Definitions for first-passage times under advective--diffusive transport in a
 body centered cubic beadpack, based on the minimal periodic unit cel, in
 parallell. */
namespace model_bcc_symmetryplanes_advection_diffusion_fpt_parallel {
struct Model {
  inline static const std::string name{
      "bcc_cartesian_advection_diffusion_fpt_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using Geometry = Geometry_Bcc<meta::ParallelOptions::Parallel,
                              Periodicity::Type::symmetryplanes,
                              Dynamics::Type::firstpassage>;
using model_advection_diffusion_2d_parallel::makeTransitions;
using model_bcc_cartesian_advection_diffusion_parallel::CTRW;
using model_bcc_symmetryplanes_advection_diffusion_fpt::Info;
using model_bcc_symmetryplanes_advection_diffusion_fpt::InitialCondition;
using model_bcc_symmetryplanes_advection_diffusion_fpt::Output;
using model_bcc_symmetryplanes_advection_diffusion_fpt::Reaction;
using model_bcc_symmetryplanes_advection_diffusion_fpt::Solvers;
using model_bcc_symmetryplanes_advection_diffusion_fpt::State;
using model_bcc_symmetryplanes_advection_diffusion_fpt::Transport;
} // namespace model_bcc_symmetryplanes_advection_diffusion_fpt_parallel

/** \namespace ptof::model_bcc_symmetryplanes_advection_parallel
 Definitions for advective transport in a body centered cubic beadpack, based on
 the minimal periodic unit cell, in parallel. */
namespace model_bcc_symmetryplanes_advection_parallel {
struct Model {
  inline static const std::string name{"bcc_symmetryplanes_advection_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d_parallel::makeTransitions;
using model_bcc_cartesian_advection_parallel::CTRW;
using model_bcc_symmetryplanes_advection::Info;
using model_bcc_symmetryplanes_advection::InitialCondition;
using model_bcc_symmetryplanes_advection::Output;
using model_bcc_symmetryplanes_advection::Reaction;
using model_bcc_symmetryplanes_advection::Solvers;
using model_bcc_symmetryplanes_advection::State;
using model_bcc_symmetryplanes_advection::Transport;
using model_bcc_symmetryplanes_advection_diffusion_parallel::Geometry;
} // namespace model_bcc_symmetryplanes_advection_parallel

/** \namespace ptof::model_bcc_symmetryplanes_advection_diffusion_surface_decay
 Definitions for advective--diffusive transport with surface reaction \f$ A_F +
 B_S \to B_S\f$ in a body centered cubic beadpack, based on the minimal periodic
 unit cell, in parallel. */
namespace model_bcc_symmetryplanes_advection_diffusion_surface_decay_parallel {
struct Model {
  inline static const std::string name{
      "bcc_symmetryplanes_advection_diffusion_surface_decay_parallel"};

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------\n"
              "Model\n"
              "--------------------------------------------------\n" +
                  name +
                  "\n"
                  "--------------------------------------------------\n";
  }
};

using model_advection_diffusion_2d_parallel::makeTransitions;
using model_bcc_symmetryplanes_advection_diffusion_parallel::CTRW;
using model_bcc_symmetryplanes_advection_diffusion_parallel::Geometry;
using model_bcc_symmetryplanes_advection_diffusion_surface_decay::Info;
using model_bcc_symmetryplanes_advection_diffusion_surface_decay::
    InitialCondition;
using model_bcc_symmetryplanes_advection_diffusion_surface_decay::Output;
using model_bcc_symmetryplanes_advection_diffusion_surface_decay::Reaction;
using model_bcc_symmetryplanes_advection_diffusion_surface_decay::Solvers;
using model_bcc_symmetryplanes_advection_diffusion_surface_decay::State;
using model_bcc_symmetryplanes_advection_diffusion_surface_decay::Transport;
} // namespace
  // model_bcc_symmetryplanes_advection_diffusion_surface_decay_parallel
} // namespace ptof

#endif /* PTOF_MODELS_PARALLEL_H */
