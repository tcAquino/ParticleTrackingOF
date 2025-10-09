/**
   \file PTOF/Model.h
   \author Tomas Aquino
   \date 22/02/2022
   \brief Type definitions to implement models.
*/

#ifndef PTOF_MODEL_H
#define PTOF_MODEL_H

#include "CTRW/CTRW.h"
#include "General/Meta.h"
#include "PTOF/DynamicsList.h"
#include "PTOF/Geometry.h"
#include "PTOF/Info.h"
#include "PTOF/InitialConditionHandler.h"
#include "PTOF/OutputHandler.h"
#include "PTOF/PeriodicityList.h"
#include "PTOF/ReactionHandler.h"
#include "PTOF/Solvers.h"
#include "PTOF/State.h"
#include "PTOF/Steppers.h"
#include "PTOF/TransportHandler.h"
#include "PTOF/TransportParameters.h"
#include <string>

/** \namespace ptof Objects and methods for ParticleTrackingOF. */
namespace ptof {
struct Model {
  inline static std::string banner() {
    return io::line() + "Model\n" + io::line();
  };

  struct advection_diffusion_2d {
    inline static const std::string name{"advection_diffusion_2d"};
    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport in 2D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_2d {
    inline static const std::string name{"advection_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective transport in 2D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct diffusion_2d {
    inline static const std::string name{"diffusion_2d"};
    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport in 2D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_diffusion_fpt_2d {
    inline static const std::string name{"advection_diffusion_fpt_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: First passage times under advective-diffusive\n"
                "             transport in 2D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Generic<2, ParallelOption, DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_fpt_2d {
    inline static const std::string name{"advection_fpt_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: First passage times under advective transport in\n"
             "             2D\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Generic<2, ParallelOption, DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct diffusion_fpt_2d {
    inline static const std::string name{"diffusion_fpt_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: First passage times under diffusive transport in\n"
             "             2D\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Generic<2, ParallelOption, DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_diffusion_surface_decay_2d {
    inline static const std::string name{
        "advection_diffusion_surface_decay_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with surface\n"
                "             reaction A_F + B_S -> B_S in 2D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct diffusion_surface_decay_2d {
    inline static const std::string name{"diffusion_surface_decay_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> B_S in 2D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_diffusion_surface_adsorption_2d {
    inline static const std::string name{
        "advection_diffusion_surface_adsorption_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with linear\n"
                "             reversible surface adsorption in 2D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<2, ParallelOption>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct diffusion_surface_adsorption_2d {
    inline static const std::string name{"diffusion_surface_adsorption_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: Diffusive transport with linear reversible surface\n"
             "             adsorption in 2D\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<2, ParallelOption>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_diffusion_surface_order2decay_2d {
    inline static const std::string name{
        "advection_diffusion_surface_order2decay_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with surface\n"
                "             reaction A_F + B_S -> Nothing in 2D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct diffusion_surface_order2decay_2d {
    inline static const std::string name{"diffusion_surface_order2decay_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> Nothing in 2D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_diffusion_2d {
    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << io::line() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport in 2D, with some\n"
                "             periodic boundaries aligned with the Cartesian\n"
                "             axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_2d {
    inline static const std::string name{"periodic_cartesian_advection_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << io::line() << "Name: " << name << "\n"
             << "Description: Advective transport in 2D, with some periodic\n"
                "             boundaries aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_diffusion_2d {
    inline static const std::string name{"periodic_cartesian_diffusion_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << io::line() << "Name: " << name << "\n"
             << "Description: Diffusive transport in 2D, with some periodic\n"
                "             boundaries aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_diffusion_fpt_2d {
    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_fpt_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << io::line() << "Name: " << name << "\n"
             << "Description: First-passage times under advective-diffusive\n"
                "             transport in 2D, with some periodic boundaries\n"
                "             aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Periodic_Cartesian<2, ParallelOption,
                                      DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_fpt_2d {
    inline static const std::string name{"periodic_cartesian_advection_fpt_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << io::line() << "Name: " << name << "\n"
          << "Description: First-passage times under advective transport in\n"
             "             2D, with some periodic boundaries aligned with\n"
             "             the Cartesian axes\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Periodic_Cartesian<2, ParallelOption,
                                      DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_diffusion_fpt_2d {
    inline static const std::string name{"periodic_cartesian_diffusion_fpt_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << io::line() << "Name: " << name << "\n"
          << "Description: First-passage times under diffusive transport in\n"
             "             2D, with some periodic boundaries aligned with\n"
             "             the Cartesian axes\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Periodic_Cartesian<2, ParallelOption,
                                      DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_diffusion_surface_decay_2d {
    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_decay_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with surface\n"
                "             reaction A_F + B_S -> B_S in 2D, with some\n"
                "             periodic boundaries aligned with the Cartesian\n"
                "             axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_diffusion_surface_decay_2d {
    inline static const std::string name{
        "periodic_cartesian_diffusion_surface_decay_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> B_S in 2D, with some periodic boundaries\n"
                "             aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_diffusion_surface_adsorption_2d {
    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_adsorption_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with surface\n"
                "             reaction A_F + B_S -> B_S in 2D, with some\n"
                "             periodic boundaries aligned with the Cartesian\n"
                "             axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_diffusion_surface_adsorption_2d {
    inline static const std::string name{
        "periodic_cartesian_diffusion_surface_adsorption_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F+\n"
                "             B_S->B_S in 2D, with some periodic boundaries\n"
                "             aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_diffusion_surface_order2decay_2d {
    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_order2decay_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with surface\n"
                "             reaction A_F + B_S -> Nothing in 2D, with some\n"
                "             periodic boundaries aligned with the Cartesian\n"
                "             axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_diffusion_surface_order2decay_2d {
    inline static const std::string name{
        "periodic_cartesian_diffusion_surface_order2decay_2d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> Nothing in 2D, with some periodic\n"
                "             boundaries aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<2, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_diffusion_3d {
    inline static const std::string name{"advection_diffusion_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport in 3D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_3d {
    inline static const std::string name{"advection_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective transport in 3D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct diffusion_3d {
    inline static const std::string name{"diffusion_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport in 3D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_diffusion_fpt_3d {
    inline static const std::string name{"advection_diffusion_fpt_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: First-passage times under advective-diffusive\n"
                "             transport in 3D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Generic<3, ParallelOption, DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_fpt_3d {
    inline static const std::string name{"advection_fpt_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: First-passage times under advective transport 3D\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Generic<3, ParallelOption, DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct diffusion_fpt_3d {
    inline static const std::string name{"diffusion_fpt_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: First-passage times under diffusive transport in\n"
             "             3D\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Generic<3, ParallelOption, DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_diffusion_surface_decay_3d {
    inline static const std::string name{
        "advection_diffusion_surface_decay_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with surface\n"
                "             reaction A_F + B_S -> B_S in 3D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct diffusion_surface_decay_3d {
    inline static const std::string name{"diffusion_surface_decay_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> B_S in 3D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_diffusion_surface_adsorption_3d {
    inline static const std::string name{
        "advection_diffusion_surface_adsorption_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with surface\n"
                "             reaction A_F + B_S -> B_S in 3D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<3, ParallelOption>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct diffusion_surface_adsorption_3d {
    inline static const std::string name{"diffusion_surface_adsorption_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> B_S in 3D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<3, ParallelOption>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct advection_diffusion_surface_order2decay_3d {
    inline static const std::string name{
        "advection_diffusion_surface_order2decay_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with surface\n"
                "             reaction A_F + B_S -> Nothing in 3D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct diffusion_surface_order2decay_3d {
    inline static const std::string name{"diffusion_surface_order2decay_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> Nothing in 3D\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Generic<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Generic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_diffusion_3d {
    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport in 3D, with some\n"
                "             periodic boundaries aligned with the Cartesian\n"
                "             axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_3d {
    inline static const std::string name{"periodic_cartesian_advection_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective transport in 3D, with some periodic\n"
                "             boundaries aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_diffusion_3d {
    inline static const std::string name{"periodic_cartesian_diffusion_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport in 3D, with some periodic\n"
                "             boundaries aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_diffusion_fpt_3d {
    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_fpt_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: First-passage times under advective-diffusive\n"
                "             transport in 3D, with some periodic boundaries\n"
                "             aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Periodic_Cartesian<3, ParallelOption,
                                      DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_fpt_3d {
    inline static const std::string name{"periodic_cartesian_advection_fpt_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: First-passage times under advective transport in\n"
             "             3D, with some periodic boundaries aligned with the\n"
             "             Cartesian axes\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Periodic_Cartesian<3, ParallelOption,
                                      DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_diffusion_fpt_3d {
    inline static const std::string name{"periodic_cartesian_diffusion_fpt_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: First-passage times under diffusive transport in\n"
             "             3D, with some periodic boundaries aligned with the\n"
             "             Cartesian axes\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Periodic_Cartesian<3, ParallelOption,
                                      DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_diffusion_surface_decay_3d {
    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_decay_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with surface\n"
                "             reaction A_F + B_S -> B_S in 3D, with some\n"
                "             periodic boundaries aligned with the Cartesian\n"
                "             axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_diffusion_surface_decay_3d {
    inline static const std::string name{
        "periodic_cartesian_diffusion_surface_decay_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> B_S in 3D, with some periodic boundaries\n"
                "             aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_diffusion_surface_adsorption_3d {
    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_adsorption_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with surface\n"
                "             reaction A_F + B_S -> B_S in 3D, with some\n"
                "             periodic boundaries aligned with the Cartesian\n"
                "             axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_diffusion_surface_adsorption_3d {
    inline static const std::string name{
        "periodic_cartesian_diffusion_surface_adsorption_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> B_S in 3D, with some periodic boundaries\n"
                "             aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_advection_diffusion_surface_order2decay_3d {
    inline static const std::string name{
        "periodic_cartesian_advection_diffusion_surface_order2decay_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective-diffusive transport with surface\n"
                "             reaction A_F + B_S -> Nothing in 3D, with some\n"
                "             periodic boundaries aligned with the Cartesian\n"
                "             axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct periodic_cartesian_diffusion_surface_order2decay_3d {
    inline static const std::string name{
        "periodic_cartesian_diffusion_surface_order2decay_3d"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> Nothing in 3D, with some periodic\n"
                "             boundaries aligned with the Cartesian axes\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Periodic_Cartesian<3, ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Generic<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_advection_diffusion {
    inline static const std::string name{"bcc_cartesian_advection_diffusion"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport in a body centered\n"
             "             cubic beadpack, based on the primitive unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Bcc<ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_advection {
    inline static const std::string name{"bcc_cartesian_advection"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective transport in a body centered cubic\n"
                "             beadpack, based on the primitive unit cell\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Bcc<ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_diffusion {
    inline static const std::string name{"bcc_cartesian_diffusion"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport in a body centered cubic\n"
                "             beadpack, based on the primitive unit cell\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Bcc<ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_advection_diffusion_fpt {
    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_fpt"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: First-passage times under advective-diffusive\n"
                "             transport in a body centered cubic beadpack,\n"
                "             based on the primitive unit cell\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::cartesian,
                       DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_advection_fpt {
    inline static const std::string name{"bcc_cartesian_advection_fpt"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: First-passage times under advective transport in a\n"
             "             body centered cubic beadpack, based on the\n"
             "             primitive unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::cartesian,
                       DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_diffusion_fpt {
    inline static const std::string name{"bcc_cartesian_diffusion_fpt"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: First-passage times under diffusive transport in\n"
             "             body centered cubic beadpack, based on the\n"
             "             primitive unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::cartesian,
                       DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_advection_diffusion_surface_decay {
    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_surface_decay"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: Advection-diffusive transport with surface\n"
             "             reaction A_F + B_S -> B_S in a body centered cubic\n"
             "             beadpack, based on the primitive unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Bcc<ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_diffusion_surface_decay {
    inline static const std::string name{
        "bcc_cartesian_diffusion_surface_decay"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> B_S in a body centered cubic beadpack,\n"
                "             based on the primitive unit cell\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Bcc<ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_advection_diffusion_surface_adsorption {
    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_surface_adsorption"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: Advection-diffusive transport with surface\n"
             "             reaction A_F + B_S -> B_S in a body centered cubic\n"
             "             beadpack, based on the primitive unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Bcc<ParallelOption>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_diffusion_surface_adsorption {
    inline static const std::string name{
        "bcc_cartesian_diffusion_surface_adsorption"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> B_S in a body centered cubic beadpack,\n"
                "             based on the primitive unit cell\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Bcc<ParallelOption>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_advection_diffusion_surface_order2decay {
    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_surface_order2decay"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport with surface\n"
             "             reaction A_F + B_S -> Nothing in a body centered\n"
             "             cubic beadpack, based on the primitive unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Bcc<ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_cartesian_diffusion_surface_order2decay {
    inline static const std::string name{
        "bcc_cartesian_diffusion_surface_order2decay"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: Diffusive transport with surface reaction A_F +\n"
             "             B_S -> Nothing in a body centered cubic beadpack,\n"
             "             based on the primitive unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry = Geometry_Bcc<ParallelOption>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_advection_diffusion {
    inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport in a body centered\n"
             "             cubic beadpack, based on the minimal unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_advection {
    inline static const std::string name{"bcc_symmetryplanes_advection"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Advective transport in a body centered cubic\n"
                "             beadpack, based on the minimal unit cell\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_diffusion {
    inline static const std::string name{"bcc_symmetryplanes_diffusion"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: Diffusive transport in a body centered cubic\n"
                "             beadpack, based on the minimal unit cell\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_advection_diffusion_fpt {
    inline static const std::string name{
        "bcc_cartesian_advection_diffusion_fpt"};

    inline static std::ostream &info(std::ostream &output) {
      output << banner() << "Name: " << name << "\n"
             << "Description: First-passage times under advective-diffusive\n"
                "             transport in a body centered cubic beadpack,\n"
                "             based on the minimal unit cell\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes,
                       DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_advection_fpt {
    inline static const std::string name{"bcc_cartesian_advection_fpt"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: First-passage times under advective transport in a\n"
             "             body centered cubic beadpack, based on the minimal\n"
             "             unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes,
                       DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::RK4, meta::Empty,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_diffusion_fpt {
    inline static const std::string name{"bcc_cartesian_diffusion_fpt"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: First-passage times under advective transport in a\n"
             "             body centered cubic beadpack, based on the minimal\n"
             "             unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes,
                       DynamicsList::Type::firstpassage>;
      using Info = Info_absorbed_reinjections;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler = ReactionHandler_NoBulk_NoSurface;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_advection_diffusion_surface_decay {
    inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion_surface_decay"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << io::line() << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport with surface\n"
             "             reaction A_F + B_S -> B_S in a body centered cubic\n"
             "             beadpack, based on the minimal unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_diffusion_surface_decay {
    inline static const std::string name{
        "bcc_symmetryplanes_diffusion_surface_decay"};

    inline static std::ostream &info(std::ostream &output) {
      output << io::line() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> B_S in a body centered cubic beadpack,\n"
                "             based on the minimal unit cell\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, false, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_advection_diffusion_surface_adsorption {
    inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion_surface_adsorption"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << io::line() << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport with surface\n"
             "             reaction A_F + B_S -> B_S in a body centered cubic\n"
             "             beadpack, based on the minimal unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_diffusion_surface_adsorption {
    inline static const std::string name{
        "bcc_symmetryplanes_diffusion_surface_adsorption"};

    inline static std::ostream &info(std::ostream &output) {
      output << io::line() << "Name: " << name << "\n"
             << "Description: Diffusive transport with surface reaction A_F +\n"
                "             B_S -> B_S in a body centered cubic beadpack,\n"
                "             based on the minimal unit cell\n"
             << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes>;
      using Info = Info_absorbed_adsorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<meta::Empty, Steppers::Euler,
                                      CTRWSteppers::Asynchronous>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceAdsorption<Geometry, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_advection_diffusion_surface_order2decay {
    inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion_surface_order2decay"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: Advective-diffusive transport with surface\n"
             "             reaction A_F + B_S -> Nothing in a body centered\n"
             "             cubic beadpack, based on the minimal unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };

  struct bcc_symmetryplanes_diffusion_surface_order2decay {
    inline static const std::string name{
        "bcc_symmetryplanes_diffusion_surface_order2decay"};

    inline static std::ostream &info(std::ostream &output) {
      output
          << banner() << "Name: " << name << "\n"
          << "Description: Diffusive transport with surface reaction A_F +\n"
             "             B_S -> Nothing in a body centered cubic beadpack,\n"
             "             based on the minimal unit cell\n"
          << io::line();
      return output;
    }

    template <typename ParallelOption> struct Definitions {
      using Geometry =
          Geometry_Bcc<ParallelOption, PeriodicityList::Type::symmetryplanes>;
      using Info = Info_absorbed;
      using State =
          State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
      using CTRW = ctrw::CTRW<State, ParallelOption>;
      using Solvers = Solvers_Generic<Steppers::Euler, Steppers::Euler,
                                      CTRWSteppers::ParticleTime>;
      using TransportHandler =
          TransportHandler_LinearInterp<TransportParameters_Bcc<
              Solvers::Parameters::advection, Solvers::Parameters::diffusion>>;
      using ReactionHandler =
          ReactionHandler_NoBulk_SurfaceDecay<Geometry, true, ParallelOption>;
      using InitialConditionHandler = InitialConditionHandler_Generic;
      using OutputHandler = OutputHandler_Generic;
    };
  };
};
} // namespace ptof

#endif /* PTOF_MODEL::H */
