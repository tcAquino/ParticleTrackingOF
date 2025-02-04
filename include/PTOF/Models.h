/**
   \file PTOF/Models.h
   \author Tomás Aquino
   \date 22/02/2022
   \brief Type definitions to implement models.
*/

#ifndef PTOF_MODELS_H
#define PTOF_MODELS_H

#include "CTRW/CTRW.h"
#include "CTRW/Meta.h"
#include "General/Meta.h"
#include "PTOF/Advection.h"
#include "PTOF/Geometry.h"
#include "PTOF/Info.h"
#include "PTOF/InitialCondition_Cases.h"
#include "PTOF/Output_Cases.h"
#include "PTOF/ReactionHandler.h"
#include "PTOF/Solvers.h"
#include "PTOF/State.h"
#include "PTOF/Steppers.h"
#include "PTOF/Transport.h"
#include "PTOF/Useful.h"
#include <functional>
#include <ostream>
#include <string>
#include <type_traits>
#include <vector>

/** \namespace ptof Objects and methods for ParticleTrackingOF. */
namespace ptof {
/**
   \namespace model_advection_diffusion_2d Definitions for 2D
   advective--diffusive transport.
*/
namespace model_advection_diffusion_2d {
template <typename ParallelOption> struct Definitions {
  Definitions() = delete;

  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_diffusion_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport in 2D\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Geometry = Geometry_Generic<2, ParallelOption>;
  using Info = Info_Absorbed_Patch;
  using State = State_Generic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using Solvers = Solvers_Generic<Stepper::Euler, Stepper::Euler,
                                  CTRWStepper::Asynchronous>;
  using Transport = Transport_Generic<Solvers>;
  using Reaction = ReactionHandler_NoBulk_NoSurface;

  struct InitialCondition {
    InitialCondition() = delete;

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

  struct Output {
    Output() = delete;

    using Parameters = OutputParameters_Cases;

    template <typename Subject, typename VelocityField, typename Geometry,
              typename Mask = meta::Empty>
    static auto
    makeOutput(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters const &params_output, std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks = {},
               std::vector<double> thresholds = {}) {
      return Output_Cases{subject,       velocity_field, geometry, directories,
                          params_output, identifier,     masks,    thresholds};
    }

    template <typename Subject, typename VelocityField, typename Geometry,
              typename Mask>
    static auto
    makeOutput(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters const &params_output, std::string const &identifier,
               std::initializer_list<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {}) {
      return Output_Cases{subject,       velocity_field, geometry, directories,
                          params_output, identifier,     masks,    thresholds};
    }

    template <typename Subject, typename VelocityField, typename Geometry,
              typename Mask>
    static auto
    makeOutput(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters const &params_output, std::string const &identifier,
               std::vector<std::reference_wrapper<const Mask>> masks,
               std::initializer_list<double> thresholds = {}) {
      return Output_Cases{subject,       velocity_field, geometry, directories,
                          params_output, identifier,     masks,    thresholds};
    }

    template <typename Subject, typename VelocityField, typename Geometry,
              typename Mask>
    static auto
    makeOutput(Subject const &subject, VelocityField const &velocity_field,
               Geometry const &geometry, Directories const &directories,
               Parameters const &params_output, std::string const &identifier,
               std::initializer_list<std::reference_wrapper<const Mask>> masks,
               std::vector<double> thresholds = {}) {
      return Output_Cases{subject,       velocity_field, geometry, directories,
                          params_output, identifier,     masks,    thresholds};
    }
  };
};
} // namespace model_advection_diffusion_2d

/**
   \namespace model_advection_2d Definitions for 2D advective transport.
*/
namespace model_advection_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output << "------------------------------------------------------------"
                "--\n"
                "Model\n"
                "------------------------------------------------------------"
                "--\n"
             << "Name: " << name << "\n"
             << "Description: Advective transport in 2D\n"
             << "Parallelization: "
             << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                     ? "Serial"
                     : "Parallel")
             << "\n"
             << "------------------------------------------------------------"
                "--\n";
      return output;
    }
  };

  using CTRW =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::CTRW;
  using Geometry = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Geometry;
  using Info =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Reaction;
  using State =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::State;

  using Solvers =
      Solvers_Generic<Stepper::Euler, meta::Empty, CTRWStepper::Asynchronous>;
  using Transport = Transport_Generic<Solvers>;
};
} // namespace model_advection_2d

/**
   \namespace model_advection_diffusion_fpt_2d Definitions for
   first-passage times under 2D advective--diffusive transport.
*/
namespace model_advection_diffusion_fpt_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_diffusion_fpt_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: First passage times under advective-diffusive\n"
             "             transport in 2D\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Geometry =
      Geometry_Generic<2, ParallelOption, Dynamics::Type::firstpassage>;
  using Info = Info_Absorbed_Reinjections;
  using State = State_Generic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Solvers;
  using Transport = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Transport;
};
} // namespace model_advection_diffusion_fpt_2d

/**
   \namespace model_advection_diffusion_surface_decay_2d Definitions for
   2D advective--diffusive transport with surface reaction \f$ A_F
   + B_S \to B_S\f$.
*/
namespace model_advection_diffusion_surface_decay_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "advection_diffusion_surface_decay_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport with surface\n"
             "             reaction A_F + B_S -> B_S\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::CTRW;
  using Geometry = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Geometry;
  using Info =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Output;
  using Solvers = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Solvers;
  using State =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::State;
  using Transport = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Transport;
  using Reaction = ReactionHandler_NoBulk_SurfaceDecay<false>;
};
} // namespace model_advection_diffusion_surface_decay_2d

/**
   \namespace model_advection_diffusion_surface_order2_2d Definitions for
   2D advective--diffusive transport with surface reaction \f$ A_F
   + B_S \to \emptyset\f$.
*/
namespace model_advection_diffusion_surface_order2_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "advection_diffusion_surface_order2_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport with surface\n"
             "             reaction A_F + B_S -> Nothing\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::CTRW;
  using Geometry = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Geometry;
  using Info =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Output;
  using Solvers = Solvers_Generic<Stepper::Euler, Stepper::Euler,
                                  CTRWStepper::TimeStep>;
  using State =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::State;
  using Transport = Transport_Generic<Solvers>;
  using Reaction = ReactionHandler_NoBulk_SurfaceDecay<true>;
};
} // namespace model_advection_diffusion_surface_order2_2d

/**
   \namespace model_periodic_cartesian_advection_diffusion_2d Definitions
   for 2D advective--diffusive transport with some periodic boundaries aligned
   with the Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport in 2D, with some\n"
             "             periodic boundaries aligned with the Cartesian\n"
             "             axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption>;
  using Info =
      typename model_advection_diffusion_2d::Definitions<ParallelOption>::Info;
  using State =
      State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Solvers;
  using Transport = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Transport;
};
} // namespace model_periodic_cartesian_advection_diffusion_2d

/**
   \namespace model_periodic_cartesian_advection_2d Definitions for 2D
   advective transport with some periodic boundaries aligned with the Cartesian
   axes.
*/
namespace model_periodic_cartesian_advection_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"periodic_cartesian_advection_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advective transport in 2D, with some periodic\n"
             "             boundaries aligned with the Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using InitialCondition = typename model_advection_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output =
      typename model_advection_2d::Definitions<ParallelOption>::Output;
  using Reaction =
      typename model_advection_2d::Definitions<ParallelOption>::Reaction;
  using Solvers =
      typename model_advection_2d::Definitions<ParallelOption>::Solvers;
  using Transport =
      typename model_advection_2d::Definitions<ParallelOption>::Transport;
  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Geometry;
  using mInfo =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Info;
  using State =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::State;
};
} // namespace model_periodic_cartesian_advection_2d

/**
   \namespace model_periodic_cartesian_advection_diffusion_fpt_2d
   Definitions for first-passage times under 2D advective--diffusive transport
   with some periodic boundaries aligned with the Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_fpt_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_fpt_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: First-passage times under advective-diffusive\n"
             "             transport in 2D, with some periodic boundaries\n"
             "             aligned with the Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption,
                                               Dynamics::Type::firstpassage>;
  using Info = typename model_advection_diffusion_fpt_2d::Definitions<
      ParallelOption>::Info;
  using State =
      State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Output;
  using Reaction =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Reaction;
  using Solvers =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Solvers;
  using Transport =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_periodic_cartesian_advection_diffusion_fpt_2d

/**
   \namespace model_periodic_cartesian_advection_diffusion_surface_decay_2d
   Definitions for first-passage times under 2D advective--diffusive transport
   with surface reaction \f$ A_F + B_S \to B_S\f$, with some periodic boundaries
   aligned with the Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_surface_decay_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_decay_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: First-passage times under advection-diffusion in\n"
             "             2D, with surface reaction A_F + B_S -> B_S, with\n"
             "             some periodic boundaries aligned with the\n"
             "             Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Reaction =
      typename model_advection_diffusion_surface_decay_2d::Definitions<
          ParallelOption>::Reaction;
  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Geometry;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Output;
  using Solvers =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Solvers;
  using State =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::State;
  using Transport =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_periodic_cartesian_advection_diffusion_surface_decay_2d

/**
   \namespace model_periodic_cartesian_advection_diffusion_order2_2d
   Definitions for first-passage times under 2D advective--diffusive transport
   with surface reaction \f$ A_F + B_S \to \emptyset\f$, with some periodic
   boundaries aligned with the Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_surface_order2_2d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_order2_2d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: First-passage times under advection-diffusion in\n"
             "             2D, with surface reaction A_F + B_S -> Nothing,\n"
             "             with some periodic boundaries aligned with the\n"
             "             Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Geometry;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::Output;
  using Solvers =
      typename model_advection_diffusion_surface_order2_2d::Definitions<
          ParallelOption>::Solvers;
  using State =
      typename model_periodic_cartesian_advection_diffusion_2d::Definitions<
          ParallelOption>::State;
  using Transport =
      typename model_advection_diffusion_surface_order2_2d::Definitions<
          ParallelOption>::Transport;
  using Reaction =
      typename model_advection_diffusion_surface_order2_2d::Definitions<
          ParallelOption>::Reaction;
};
} // namespace model_periodic_cartesian_advection_diffusion_surface_order2_2d

/**
   \namespace model_advection_diffusion_3d Definitions for 3D
   advective--diffusive transport.
*/
namespace model_advection_diffusion_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_diffusion_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in 3D\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Geometry = Geometry_Generic<3, ParallelOption>;
  using Info = Info_Absorbed_Patch;
  using State = State_Generic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Solvers;
  using Transport = typename model_advection_diffusion_2d::Definitions<
      ParallelOption>::Transport;
};
} // namespace model_advection_diffusion_3d

/**
   \namespace model_advection_3d Definitions for 3D advective transport.
*/
namespace model_advection_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection in 3D\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Solvers =
      typename model_advection_2d::Definitions<ParallelOption>::Solvers;
  using Transport =
      typename model_advection_2d::Definitions<ParallelOption>::Transport;
  using CTRW =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::CTRW;
  using Geometry = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Geometry;
  using Info =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Reaction;
  using State =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::State;
};
} // namespace model_advection_3d

/**
   \namespace model_advection_diffusion_fpt_3d Definitions for
   first-passage times under 3D advective--diffusive transport.
*/
namespace model_advection_diffusion_fpt_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"advection_diffusion_fpt_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: First-passage times under advection-diffusion in\n"
             "             3D\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Geometry =
      Geometry_Generic<3, ParallelOption, Dynamics::Type::firstpassage>;
  using Info = typename model_advection_diffusion_fpt_2d::Definitions<
      ParallelOption>::Info;
  using State = State_Generic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Solvers;
  using Transport = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Transport;
};
} // namespace model_advection_diffusion_fpt_3d

/**
   \namespace model_advection_diffusion_surface_decay_3d Definitions for
   3D advective--diffusive transport with surface reaction \f$ A_F + B_S \to
   B_S\f$.
*/
namespace model_advection_diffusion_surface_decay_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "advection_diffusion_surface_decay_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Description: Advection-diffusion in 3D, with surface reaction\n"
             "             A_F + B_S -> B_S\n"
          << "Name: " << name << "\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::CTRW;
  using Geometry = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Geometry;
  using Info =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Output;
  using Solvers = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Solvers;
  using State =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::State;
  using Transport = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Transport;
  using Reaction =
      typename model_advection_diffusion_surface_decay_2d::Definitions<
          ParallelOption>::Reaction;
};
} // namespace model_advection_diffusion_surface_decay_3d

/**
   \namespace model_advection_diffusion_surface_order2_3d Definitions for
   3D advective--diffusive transport with surface reaction \f$ A_F + B_S \to
   \emptyset\f$.
*/
namespace model_advection_diffusion_surface_order2_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "advection_diffusion_surface_order2_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Description: Advection-diffusion in 3D, with surface reaction\n"
             "             A_F + B_S -> Nothing\n"
          << "Name: " << name << "\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::CTRW;
  using Geometry = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Geometry;
  using Info =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Output;
  using Solvers =
      typename model_advection_diffusion_surface_order2_2d::Definitions<
          ParallelOption>::Solvers;
  using State =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::State;
  using Transport =
      typename model_advection_diffusion_surface_order2_2d::Definitions<
          ParallelOption>::Transport;
  using Reaction =
      typename model_advection_diffusion_surface_order2_2d::Definitions<
          ParallelOption>::Reaction;
};
} // namespace model_advection_diffusion_surface_order2_3d

/**
   \namespace model_periodic_cartesian_advection_diffusion_3d Definitions
   for 3D advective transport with some periodic boundaries aligned with the
   Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in 3D, with some periodic\n"
             "             boundaries aligned with the Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption>;
  using Info =
      typename model_advection_diffusion_3d::Definitions<ParallelOption>::Info;
  using State =
      State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Solvers;
  using Transport = typename model_advection_diffusion_3d::Definitions<
      ParallelOption>::Transport;
};
} // namespace model_periodic_cartesian_advection_diffusion_3d

/**
   \namespace model_periodic_cartesian_advection_3d
   Definitions for 3D advective transport with some periodic boundaries aligned
   with the Cartesian axes.
*/
namespace model_periodic_cartesian_advection_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"periodic_cartesian_advection_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection in 3D, with some periodic boundaries\n"
             "             aligned with the Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Reaction =
      typename model_advection_3d::Definitions<ParallelOption>::Reaction;
  using Solvers =
      typename model_advection_3d::Definitions<ParallelOption>::Solvers;
  using Transport =
      typename model_advection_3d::Definitions<ParallelOption>::Transport;
  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Geometry;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Output;
  using State =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::State;
};
} // namespace model_periodic_cartesian_advection_3d

/**
   \namespace model_periodic_cartesian_advection_diffusion_fpt_3d
   Definitions for first-passage times under 3D advective--diffusive transport
   with some periodic boundaries aligned with the Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_fpt_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_fpt_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: First-passage times under advection-diffusion in\n"
             "             3D, with some periodic boundaries aligned with the\n"
             "             Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption,
                                               Dynamics::Type::firstpassage>;
  using Info = typename model_advection_diffusion_fpt_3d::Definitions<
      ParallelOption>::Info;
  using State =
      State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
  using CTRW = ctrw::CTRW<State, ParallelOption>;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Output;
  using Reaction =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Reaction;
  using Solvers =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Solvers;
  using Transport =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_periodic_cartesian_advection_diffusion_fpt_3d

/**
   \namespace model_periodic_cartesian_advection_diffusion_surface_decay_3d
   Definitions for 3D advective--diffusive transport with surface reaction \f$
   A_F + B_S \to B_S\f$, with some periodic boundaries aligned with the
   Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_surface_decay_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_decay_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in 3D, with surface reaction\n"
             "             A_F + B_S -> B_S, with some periodic boundaries\n"
             "             aligned with the Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Geometry;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Output;
  using Solvers =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Solvers;
  using State =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::State;
  using Transport =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Transport;
  using Reaction =
      typename model_advection_diffusion_surface_decay_3d::Definitions<
          ParallelOption>::Reaction;
};
} // namespace model_periodic_cartesian_advection_diffusion_surface_decay_3d

/**
   \namespace model_periodic_cartesian_advection_diffusion_surface_order2_3d
   Definitions for 3D advective--diffusive transport with surface reaction \f$
   A_F + B_S \to B_S\f$, with some periodic boundaries aligned with the
   Cartesian axes.
*/
namespace model_periodic_cartesian_advection_diffusion_surface_order2_3d {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_order2_3d"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in 3D, with surface reaction\n"
             "             A_F + B_S -> Nothing, with some periodic\n"
             "             boundaries aligned with the Cartesian axes\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Geometry;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Output;
  using Solvers =
      typename model_advection_diffusion_surface_order2_3d::Definitions<
          ParallelOption>::Solvers;
  using State =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::State;
  using Transport =
      typename model_advection_diffusion_surface_order2_3d::Definitions<
          ParallelOption>::Transport;
  using Reaction =
      typename model_advection_diffusion_surface_order2_3d::Definitions<
          ParallelOption>::Reaction;
};
} // namespace model_periodic_cartesian_advection_diffusion_surface_order2_3d

/**
   \namespace model_bcc_cartesian_advection_diffusion Definitions for
   advective--diffusive transport in a body centered cubic beadpack, based on
   the primitive unit cell.
*/
namespace model_bcc_cartesian_advection_diffusion {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"bcc_cartesian_advection_diffusion"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in a body centered cubic\n"
             "             beadpack, based on the primitive unit cell\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Geometry = Geometry_Bcc<ParallelOption>;
  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::CTRW;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Output;
  using Reaction =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Reaction;
  using Solvers =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::Solvers;
  using State =
      typename model_periodic_cartesian_advection_diffusion_3d::Definitions<
          ParallelOption>::State;
  using Transport = Transport_Bcc<Solvers>;
};
} // namespace model_bcc_cartesian_advection_diffusion

/**
   \namespace model_bcc_cartesian_advection_diffusion_fpt Definitions for
   first-passage times under advective--diffusive transport in a body centered
   cubic beadpack, based on the primitive unit cell.
*/
namespace model_bcc_cartesian_advection_diffusion_fpt {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_fpt"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Description: First-passage times under advection-diffusion in\n"
             "             body centered cubic beadpack, based on the\n"
             "             primitive unit cell\n"
          << "Name: " << name << "\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "-------------------------------------------------------------"
             "-\n";
      return output;
    }
  };

  using Geometry = Geometry_Bcc<ParallelOption, Periodicity::Type::cartesian,
                                Dynamics::Type::firstpassage>;
  using InitialCondition =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Output;
  using Reaction =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Reaction;
  using Solvers = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Solvers;
  using Transport =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Transport;
  using CTRW =
      typename model_periodic_cartesian_advection_diffusion_fpt_3d::Definitions<
          ParallelOption>::CTRW;
  using Info =
      typename model_periodic_cartesian_advection_diffusion_fpt_3d::Definitions<
          ParallelOption>::Info;
  using State =
      typename model_periodic_cartesian_advection_diffusion_fpt_3d::Definitions<
          ParallelOption>::State;
};
} // namespace model_bcc_cartesian_advection_diffusion_fpt

/**
   \namespace model_bcc_cartesian_advection Definitions for advective
   transport in a body centered cubic beadpack, based on the primitive unit
   cell.
*/
namespace model_bcc_cartesian_advection {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"bcc_cartesian_advection"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Description: Advection in a body centered cubic beadpack, based\n"
             "             on the primitive unit cell\n"
          << "Name: " << name << "\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::CTRW;
  using Geometry =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Geometry;
  using Info = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Info;
  using InitialCondition =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Output;
  using State = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::State;
  using Reaction = typename model_periodic_cartesian_advection_3d::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_periodic_cartesian_advection_3d::Definitions<
      ParallelOption>::Solvers;
  using Transport = Transport_Bcc<Solvers>;
};
} // namespace model_bcc_cartesian_advection

/**
   \namespace model_bcc_cartesian_advection_diffusion_surface_decay
   Definitions for advective--diffusive transport with surface reaction \f$ A_F
   + B_S \to B_S\f$ in a body centered cubic beadpack, based on the primitive
   unit cell.
*/
namespace model_bcc_cartesian_advection_diffusion_surface_decay {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_surface_decay"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in a body centered cubic\n"
             "             beadpack, based on the primitive unit cell, with\n"
             "             surface reaction A_F + B_S -> B_S\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::CTRW;
  using Geometry =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Geometry;
  using Info = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Info;
  using InitialCondition =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Output;
  using Solvers = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Solvers;
  using State = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::State;
  using Transport =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Transport;
  using Reaction =
      typename model_periodic_cartesian_advection_diffusion_surface_decay_3d::
          Definitions<ParallelOption>::Reaction;
};
} // namespace model_bcc_cartesian_advection_diffusion_surface_decay

/**
   \namespace model_bcc_cartesian_advection_diffusion_surface_decay
   Definitions for advective--diffusive transport with surface reaction \f$ A_F
   + B_S \to \emptyset\f$ in a body centered cubic beadpack, based on the
   primitive unit cell.
*/
namespace model_bcc_cartesian_advection_diffusion_surface_order2 {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_surface_order2"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in a body centered cubic\n"
             "             beadpack, based on the primitive unit cell, with\n"
             "             surface reaction A_F + B_S -> Nothing\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::CTRW;
  using Geometry =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Geometry;
  using Info = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Info;
  using InitialCondition =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Output;
  using Solvers =
      typename model_periodic_cartesian_advection_diffusion_surface_decay_3d::
          Definitions<ParallelOption>::Solvers;
  using State = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::State;
  using Transport =
      typename model_periodic_cartesian_advection_diffusion_surface_decay_3d::
          Definitions<ParallelOption>::Transport;
  using Reaction =
      typename model_periodic_cartesian_advection_diffusion_surface_decay_3d::
          Definitions<ParallelOption>::Reaction;
};
} // namespace model_bcc_cartesian_advection_diffusion_surface_order2

/**
   \namespace model_bcc_symmetryplanes_advection_diffusion Definitions for
   advective--diffusive transport in a body centered cubic beadpack, based on
   the minimal unit cell.
*/
namespace model_bcc_symmetryplanes_advection_diffusion {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in a body centered cubic\n"
             "             beadpack, based on the minimal unit cell\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Geometry =
      Geometry_Bcc<ParallelOption, Periodicity::Type::symmetryplanes>;
  using CTRW = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::CTRW;
  using Info = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Info;
  using InitialCondition =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Output;
  using Reaction =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Reaction;
  using Solvers = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::Solvers;
  using State = typename model_bcc_cartesian_advection_diffusion::Definitions<
      ParallelOption>::State;
  using Transport =
      typename model_bcc_cartesian_advection_diffusion::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_bcc_symmetryplanes_advection_diffusion

/**
   \namespace model_bcc_symmetryplanes_advection_diffusion_fpt Definitions
   for first-passage times under advective--diffusive transport in a body
   centered cubic beadpack, based on the minimal unit cell.
*/
namespace model_bcc_symmetryplanes_advection_diffusion_fpt {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_fpt"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Description: First-passage times under advection-diffusion in\n"
             "             body centered cubic beadpack, based on the minimal\n"
             "             unit cell\n"
          << "Name: " << name << "\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Geometry =
      Geometry_Bcc<ParallelOption, Periodicity::Type::symmetryplanes,
                   Dynamics::Type::firstpassage>;
  using CTRW =
      typename model_bcc_cartesian_advection_diffusion_fpt::Definitions<
          ParallelOption>::CTRW;
  using Info =
      typename model_bcc_cartesian_advection_diffusion_fpt::Definitions<
          ParallelOption>::Info;
  using State =
      typename model_bcc_cartesian_advection_diffusion_fpt::Definitions<
          ParallelOption>::State;
  using InitialCondition =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Output;
  using Reaction =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Reaction;
  using Solvers =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Solvers;
  using Transport =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_bcc_symmetryplanes_advection_diffusion_fpt

/**
   \namespace model_bcc_symmetryplanes_advection Definitions for advective
   transport in a body centered cubic beadpack, based on the minimal unit cell.
*/
namespace model_bcc_symmetryplanes_advection {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{"bcc_symmetryplanes_advection"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW =
      typename model_bcc_cartesian_advection::Definitions<ParallelOption>::CTRW;
  using Info =
      typename model_bcc_cartesian_advection::Definitions<ParallelOption>::Info;
  using InitialCondition = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::InitialCondition;
  using Output = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::Output;
  using Reaction = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::Reaction;
  using Solvers = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::Solvers;
  using State = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::State;
  using Transport = typename model_bcc_cartesian_advection::Definitions<
      ParallelOption>::Transport;
  using Geometry =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Geometry;
};
} // namespace model_bcc_symmetryplanes_advection

/**
   \namespace model_bcc_symmetryplanes_advection_diffusion_surface_decay
   Definitions for advective--diffusive transport with surface reaction \f$ A_F
   + B_S \to B_S\f$ in a body centered cubic beadpack, based on the minimal unit
   cell.
*/
namespace model_bcc_symmetryplanes_advection_diffusion_surface_decay {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion_surface_decay"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in a body centered cubic\n"
             "             beadpack, based on the minimal unit cell, with\n"
             "             surface reaction A_F + B_S -> B_S\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using Reaction =
      typename model_bcc_cartesian_advection_diffusion_surface_decay::
          Definitions<ParallelOption>::Reaction;
  using CTRW =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Geometry;
  using Info =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Output;
  using Solvers =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Solvers;
  using State =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::State;
  using Transport =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Transport;
};
} // namespace model_bcc_symmetryplanes_advection_diffusion_surface_decay

/**
   \namespace model_bcc_symmetryplanes_advection_diffusion_surface_order2
   Definitions for advective--diffusive transport with surface reaction \f$ A_F
   + B_S \to \emptyset\f$ in a body centered cubic beadpack, based on the
   minimal unit cell.
*/
namespace model_bcc_symmetryplanes_advection_diffusion_surface_order2 {
template <typename ParallelOption> struct Definitions {
  struct Model {
    Model() = delete;

    inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion_surface_order2"};

    /**
       \brief Output generic information about object.
       \param output Output stream.
    */
    inline static std::ostream &info(std::ostream &output) {
      output
          << "--------------------------------------------------------------\n"
             "Model\n"
             "--------------------------------------------------------------\n"
          << "Name: " << name << "\n"
          << "Description: Advection-diffusion in a body centered cubic\n"
             "             beadpack, based on the minimal unit cell, with\n"
             "             surface reaction A_F + B_S -> Nothing\n"
          << "Parallelization: "
          << (std::is_same_v<ParallelOption, par::ParallelOptions::Serial>
                  ? "Serial"
                  : "Parallel")
          << "\n"
          << "--------------------------------------------------------------\n";
      return output;
    }
  };

  using CTRW =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::CTRW;
  using Geometry =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Geometry;
  using Info =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Info;
  using InitialCondition =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::InitialCondition;
  using Output =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::Output;
  using Solvers =
      typename model_bcc_cartesian_advection_diffusion_surface_order2::
          Definitions<ParallelOption>::Solvers;
  using State =
      typename model_bcc_symmetryplanes_advection_diffusion::Definitions<
          ParallelOption>::State;
  using Transport =
      typename model_bcc_cartesian_advection_diffusion_surface_decay::
          Definitions<ParallelOption>::Transport;
  using Reaction =
      typename model_bcc_cartesian_advection_diffusion_surface_order2::
          Definitions<ParallelOption>::Reaction;
};
} // namespace model_bcc_symmetryplanes_advection_diffusion_surface_order2
} // namespace ptof

#endif /* PTOF_MODELS_H */
