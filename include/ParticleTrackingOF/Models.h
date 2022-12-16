//
//  Models.h
//
//  Created by Tomás Aquino on 22/02/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Models_OF_h
#define Models_OF_h

#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <string>
#include <type_traits>
#include <meshSearch.H>
#include "general/useful.h"
#include "ParticleTrackingOF/Advection.h"
#include "ParticleTrackingOF/Geometry.h"
#include "ParticleTrackingOF/InitialConditions.h"
#include "ParticleTrackingOF/Info.h"
#include "ParticleTrackingOF/Locator.h"
#include "ParticleTrackingOF/Output.h"
#include "ParticleTrackingOF/Reaction.h"
#include "ParticleTrackingOF/State.h"
#include "ParticleTrackingOF/Steppers.h"
#include "ParticleTrackingOF/Transitions.h"
#include "Stochastic/CTRW/CTRW.h"

namespace ptof
{
  namespace model_advection_diffusion_2d
  {
    struct Model
    {
      inline static const std::string name{ "advection_diffusion_2d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry<2,
      BoundaryConditionSet::Type::transport>;
    using Info = ptof::Info_Absorbed;
    using State = StateDim<
      Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW<State>;
    
    struct Solvers
    {
      using Steppers = Steppers_Advection_RK4_Diffusion_Euler;
      
      struct Parameters
      {
        double time_step;
        
        template
        <typename TransportParameters,
        typename ReactionParameters>
        Parameters
        (Directories const& directories,
         std::string const& name,
         TransportParameters const& params_transport,
         ReactionParameters const& params_reaction
         )
        {
          double time_step_accuracy_adv;
          double time_step_accuracy_diff;
          double time_step_accuracy_react;
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_solvers_"
                                         + name + ".dat");
          useful::read(input, time_step_accuracy_adv);
          useful::read(input, time_step_accuracy_diff);
          useful::read(input, time_step_accuracy_react);
          input.close();
          
          time_step = std::min({
            time_step_accuracy_adv*params_transport.advection_time,
            time_step_accuracy_diff*params_transport.diffusion_time,
            time_step_accuracy_react*params_reaction.reaction_time });
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Solver parameters\n"
            "--------------------------------------------------\n"
            "- Maximum timestep in units of advection time\n"
            "- Maximum timestep in units of diffusion time\n"
            "- Maximum timestep in units of reaction time\n"
            "--------------------------------------------------\n";
        }
      };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Solvers\n"
          "--------------------------------------------------\n"
          "Advection: RK4\n"
          "Diffusion: Stochastic Euler\n"
          "--------------------------------------------------\n";
      }
    };
    
    struct Transport
    {
      struct Parameters
      {
        double lengthscale;
        double peclet;
        double diff_coeff;
        
        double diffusion_time;
        double advection_time;
        
        Parameters
        (Directories const& directories, std::string const& name)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_transport_"
                                         + name + ".dat");
          useful::read(input, lengthscale);
          useful::read(input, peclet);
          useful::read(input, diff_coeff);
          input.close();
          
          diffusion_time = lengthscale*lengthscale/(2.*diff_coeff);
          advection_time = 2.*diffusion_time/peclet;
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Transport parameters\n"
            "--------------------------------------------------\n"
            "- Reference lengthscale\n"
            "- Peclet number\n"
            "- Diffusion coefficient\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Boundary>
      static auto makeTransitions
      (VelocityField&& velocity_field,
       Geometry const& geometry,
       Boundary const& boundary,
       Parameters const& params_transport,
       Solvers::Parameters const& params_solvers)
      {
        return makeTransportTransitions<
          Geometry, Solvers::Steppers>(velocity_field,
                                       geometry,
                                       boundary,
                                       params_transport,
                                       params_solvers);
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry, double average_velocity_magnitude)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data_rescaled(geometry.mesh,
                                           average_velocity_magnitude));
      };
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data(geometry.mesh));
      };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Transport\n"
          "--------------------------------------------------\n"
          "Processes: Advection-diffusion\n"
          "Interpolation: linear\n"
          "--------------------------------------------------\n";
      }
    };
    
    struct Reaction
    {
      struct Parameters
      {
        double damkohler{ 0. };
        double reaction_rate{ 0. };
        double reaction_time{ 1./reaction_rate };
        
        template <typename TransportParameters>
        Parameters
        (Directories const& dir,
         std::string const& name,
         TransportParameters const& params_transport)
        {}
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Reaction parameters\n"
            "--------------------------------------------------\n"
            "None\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename TransportParameters,
      typename SolverParameters>
      static auto makeReaction
      (Geometry const& geometry,
       Parameters const& params,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        return useful::DoNothing{};
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Reaction\n"
          "--------------------------------------------------\n"
          "None\n"
          "--------------------------------------------------\n";
      }
    };
    
    struct InitialCondition
    {
      struct Parameters
      {
        InitialConditions::Type type;
        double initial_mass;
        std::size_t nr_particles;
        double distance_wall;
        
        double time_step_accuracy_adv;
        double time_step_accuracy_diff;
        double time_step_accuracy_react;
        
        double time_min{ 0. };
        double time_max{ 1. };
        double time_step{ 1. };
        
        template
        <typename TransportParameters,
        typename ReactionParameters,
        typename SolverParameters>
        Parameters
        (Directories const& directories,
         std::string const& name,
         TransportParameters const& params_transport,
         ReactionParameters const& params_reaction,
         SolverParameters const& params_solvers)
        : distance_wall{
          10.*std::sqrt(2.*params_transport.diff_coeff*
                        params_solvers.time_step) }
        {
          std::string ic_name;
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_initial_condition_"
                                         + name + ".dat");
          useful::read(input, ic_name);
          verify_initial_condition(ic_name,
                                   InitialConditions{});
          type = InitialConditions::type(ic_name);
          useful::read(input, initial_mass);
          useful::read(input, nr_particles);
          input >> time_min;
          if (input.fail())
            return;
          useful::read(input, time_max);
          useful::read(input, time_step_accuracy_adv);
          useful::read(input, time_step_accuracy_diff);
          useful::read(input, time_step_accuracy_react);
          
          time_step =
            std::min({
              time_step_accuracy_adv*params_transport.advection_time,
              time_step_accuracy_diff*params_transport.diffusion_time,
              time_step_accuracy_react*params_reaction.reaction_time });
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Initial condition parameters (pass '' for default in [])\n"
            "--------------------------------------------------\n"
            "- Initial condition type:\n"
            "\tuniform: Homogeneous throughout the domain\n"
            "\tflux_weighted: Flux-weighted throughout the domain\n"
            "\tuniform_inlet: Homogeneous at the inlet\n"
            "\tflux_weighted_inlet: Flux-weighted at the inlet\n"
            "\tuniform_solid: Homogeneous at the solid surface\n"
            "\tuniform_near_solid: Homogeneous at a fixed distance to the solid interface\n"
            "\tuniform_inlet_continuous: Continuous injection homogeneous at the inlet\n"
            "\tflux_weighted_inlet_continuous: Continuous injection flux-weighted at the inlet\n"
            "- Total transported mass in each injection discretization\n"
            "- Number of Lagrangian particles in each injection discretization\n"
            "- Initial injection time [0.]\n"
            "- Final injection time [1.]\n"
            "- Maximum timestep in units of advection time [unit time step]\n"
            "- Maximum timestep in units of diffusion time [unit time step]\n"
            "- Maximum timestep in units of reaction time [unit time step]\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters params,
       double time = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            time,
            params.initial_mass/params.nr_particles },
          params };
      }
    };
    
    using VelocityField = VectorField_LinearInterpolation_OF
    <Foam::volVectorField, Locator_Cell const&, 1, 1>;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
  }
  
  namespace model_advection_2d
  {
    struct Model
    {
      inline static const std::string name{ "advection_2d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };
    
    using model_advection_diffusion_2d::Geometry;
    using model_advection_diffusion_2d::Info;
    using model_advection_diffusion_2d::State;
    using model_advection_diffusion_2d::CTRW;
    using model_advection_diffusion_2d::Reaction;
    using model_advection_diffusion_2d::InitialCondition;
    using model_advection_diffusion_2d::VelocityField;
    using model_advection_diffusion_2d::Output;
    
    struct Solvers
    {
      using Steppers = Steppers_Advection_RK4_Diffusion_Euler;
      
      struct Parameters
      {
        double time_step = 0.;
        
        std::size_t nr_particles;
        double step_length;
        
        template
        <typename TransportParameters,
        typename ReactionParameters>
        Parameters
        (Directories const& directories,
         std::string const& name,
         TransportParameters const& params_transport,
         ReactionParameters const& params_reaction)
        {
          double step_length_nondim;
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_solvers_"
                                         + name + ".dat");
          useful::read(input, step_length_nondim);
          input.close();
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Solver parameters\n"
            "--------------------------------------------------\n"
            "- Approximate step length in units of characteristic lengthscale\n"
            "--------------------------------------------------\n";
        }
      };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Solvers\n"
          "--------------------------------------------------\n"
          "Advection: RK4\n"
          "--------------------------------------------------\n";
      }
    };
    
    struct Transport
    {
      struct Parameters
      {
        double peclet = std::numeric_limits<double>::infinity();
        double diff_coeff = 0.;
        double diffusion_time = std::numeric_limits<double>::infinity();
        double advection_time = 0.;
        
        double lengthscale;
        
        Parameters(Directories const& directories, std::string const& name)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_transport_"
                                         + name + ".dat");
          useful::read(input, lengthscale);
          input.close();
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Transport parameters\n"
            "--------------------------------------------------\n"
            "- Reference lengthscale\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Boundary>
      static auto makeTransitions
      (VelocityField&& velocity_field,
       Geometry const& geometry,
       Boundary const& boundary,
       Parameters const& params_transport,
       Solvers::Parameters const& params_solvers)
      {
        return makeTransportTransitions_Advection<
          Geometry, Solvers::Steppers>(velocity_field,
                                       geometry,
                                       boundary,
                                       params_transport,
                                       params_solvers);
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry, double average_velocity_magnitude)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data_rescaled(geometry.mesh,
                                           average_velocity_magnitude));
      };
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data(geometry.mesh));
      };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Transport\n"
          "--------------------------------------------------\n"
          "Process: Advection\n"
          "Interpolation: linear\n"
          "--------------------------------------------------\n";
      }
    };
  };
  
  namespace model_advection_diffusion_fpt_2d
  {
    struct Model
    {
      inline static const std::string name{
        "advection_diffusion_fpt_2d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry<2,
      BoundaryConditionSet::Type::firstpassage>;
    using Info = ptof::Info_Absorbed_Reinjections;
    using State = StateDim<Geometry::dim,
      Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW<State>;
    using model_advection_diffusion_2d::Solvers;
    using model_advection_diffusion_2d::Transport;
    using model_advection_diffusion_2d::Reaction;
    using model_advection_diffusion_2d::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    struct InitialCondition final
    : model_advection_diffusion_2d::InitialCondition
    {
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters params,
       double time = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            time,
            params.initial_mass/params.nr_particles },
          params };
      }
    };
  }
  
  namespace model_advection_diffusion_decay_catalytic_2d
  {
    struct Model
    {
      inline static const std::string name{
        "advection_diffusion_decay_catalytic_2d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_advection_diffusion_2d::Geometry;
    using model_advection_diffusion_2d::Info;
    using model_advection_diffusion_2d::State;
    using model_advection_diffusion_2d::CTRW;
    using model_advection_diffusion_2d::Solvers;
    using model_advection_diffusion_2d::Transport;
    using model_advection_diffusion_2d::InitialCondition;
    using model_advection_diffusion_2d::Output;
    
    struct Reaction
    {
      struct Parameters
      {
        double damkohler;
        std::string initial_distribution;
        double surface_concentration;
        double reaction_rate;
        double reaction_time;
      
        template <typename TransportParameters>
        Parameters
        (Directories const& directories,
         std::string const& name,
         TransportParameters const& params_transport)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_reaction_"
                                         + name + ".dat");
          useful::read(input, damkohler);
          useful::read(input, initial_distribution);
          useful::read(input, surface_concentration);
          reaction_rate = params_transport.lengthscale*damkohler/
            (surface_concentration*params_transport.diffusion_time);
          reaction_time = 1./reaction_rate;
          input.close();
        }
      
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Reaction parameters:\n"
            "--------------------------------------------------\n"
            "- Damkohler number\n"
            "- Solid reactant initial distribution type\n"
            "\tuniform: Homogeneous throughout the domain\n"
            "- Initial solid reactant surface concentration\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename TransportParameters,
      typename SolverParameters>
      static auto makeReaction
      (Geometry const& geometry,
       Parameters const& params,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        if (params.initial_distribution == "uniform")
          return Reaction_AFluidPlusASolidtoASolid{
            params.reaction_rate,
            10.*std::sqrt(2.*params_transport.diff_coeff
                          *params_solvers.time_step),
            uniform_solid_reactant_patches(params.surface_concentration,
                                               { "wallFluidSolid" },
                                               geometry.mesh ),
            geometry.mesh_search
          };
        throw std::runtime_error{
          "Initial reactant distribution "
          + params.initial_distribution
          + " not supported"
        };
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        Reaction_AFluidPlusASolidtoASolid<Foam::meshSearch>::
          info(output);
      }
    };
  }
  
  namespace model_advection_diffusion_decay_2d
  {
    struct Model
    {
      inline static const std::string name{
        "advection_diffusion_decay_2d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_advection_diffusion_2d::Geometry;
    using model_advection_diffusion_2d::Info;
    using model_advection_diffusion_2d::State;
    using model_advection_diffusion_2d::CTRW;
    using model_advection_diffusion_2d::Solvers;
    using model_advection_diffusion_2d::Transport;
    using model_advection_diffusion_2d::InitialCondition;
    using model_advection_diffusion_2d::VelocityField;
    using model_advection_diffusion_2d::Output;
    
    struct Reaction
    {
      using Parameters =
        model_advection_diffusion_decay_catalytic_2d::
        Reaction::Parameters;
      
      template
      <typename Geometry,
      typename TransportParameters,
      typename SolverParameters>
      static auto makeReaction
      (Geometry const& geometry,
       Parameters const& params,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        return Reaction_AFluidPlusASolidtoNothing{
          params.reaction_rate,
          10.*std::sqrt(2.*params_transport.diff_coeff*
                        params_solvers.time_step),
          uniform_solid_reactant_patches(params.surface_concentration,
                                             { "wallFluidSolid" },
                                             geometry.mesh ),
          geometry.mesh_search
        };
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        Reaction_AFluidPlusASolidtoNothing<Foam::meshSearch>::
          info(output);
      }
    };
  }
  
  namespace model_periodic_cartesian_advection_diffusion_2d
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_2d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry_Periodic_Cartesian<2,
      BoundaryConditionSet::Type::transport>;
    using model_advection_diffusion_2d::Info;
    using State = StateDim_Periodic<
      Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW<State>;
    using model_advection_diffusion_2d::Solvers;
    using model_advection_diffusion_2d::Transport;
    using model_advection_diffusion_2d::Reaction;
    using model_advection_diffusion_2d::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    struct InitialCondition final
    : model_advection_diffusion_2d::InitialCondition
    {
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters params,
       double time = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker_Periodic{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            geometry.boundary_periodic,
            time,
            params.initial_mass/params.nr_particles },
          params };
      }
    };
  }
  
  namespace model_periodic_cartesian_advection_2d
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_2d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };
    
    using model_periodic_cartesian_advection_diffusion_2d::Geometry;
    using model_periodic_cartesian_advection_diffusion_2d::Info;
    using model_periodic_cartesian_advection_diffusion_2d::State;
    using model_periodic_cartesian_advection_diffusion_2d::CTRW;
    using model_advection_2d::Solvers;
    using model_advection_2d::Transport;
    using model_advection_2d::Reaction;
    using model_periodic_cartesian_advection_diffusion_2d::InitialCondition;
    using model_advection_2d::VelocityField;
    using model_periodic_cartesian_advection_diffusion_2d::Output;
  };
  
  namespace model_periodic_cartesian_advection_diffusion_fpt_2d
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_fpt_2d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry_Periodic_Cartesian<2,
      BoundaryConditionSet::Type::firstpassage>;
    using model_advection_diffusion_fpt_2d::Info;
    using model_periodic_cartesian_advection_diffusion_2d::State;
    using model_periodic_cartesian_advection_diffusion_2d::CTRW;
    using model_periodic_cartesian_advection_diffusion_2d::Solvers;
    using model_periodic_cartesian_advection_diffusion_2d::Transport;
    using model_periodic_cartesian_advection_diffusion_2d::Reaction;
    using model_periodic_cartesian_advection_diffusion_2d::InitialCondition;
    using model_periodic_cartesian_advection_diffusion_2d::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
  }
  
  namespace model_periodic_cartesian_advection_diffusion_decay_catalytic_2d
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_decay_catalytic_2d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_periodic_cartesian_advection_diffusion_2d::Geometry;
    using model_periodic_cartesian_advection_diffusion_2d::Info;
    using model_periodic_cartesian_advection_diffusion_2d::State;
    using model_periodic_cartesian_advection_diffusion_2d::CTRW;
    using model_periodic_cartesian_advection_diffusion_2d::Solvers;
    using model_periodic_cartesian_advection_diffusion_2d::Transport;
    using model_periodic_cartesian_advection_diffusion_2d::InitialCondition;
    using model_periodic_cartesian_advection_diffusion_2d::Output;
    using model_advection_diffusion_decay_catalytic_2d::Reaction;
  }
  
  namespace model_periodic_cartesian_advection_diffusion_decay_2d
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_decay_2d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_periodic_cartesian_advection_diffusion_decay_catalytic_2d::Geometry;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_2d::Info;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_2d::State;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_2d::CTRW;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_2d::Solvers;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_2d::Transport;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_2d::InitialCondition;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_2d::Output;
    using model_advection_diffusion_decay_2d::Reaction;
  }
  
  namespace model_advection_diffusion_3d
  {
    struct Model
    {
      inline static const std::string name{ "advection_diffusion_3d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry<3,
      BoundaryConditionSet::Type::transport>;
    using Info = ptof::Info_Absorbed;
    using State = StateDim<
      Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW<State>;
    
    struct Solvers
    {
      using Steppers = Steppers_Advection_RK4_Diffusion_Euler;
      
      struct Parameters
      {
        double time_step;
        
        template
        <typename TransportParameters,
        typename ReactionParameters>
        Parameters
        (Directories const& directories,
         std::string const& name,
         TransportParameters const& params_transport,
         ReactionParameters const& params_reaction
         )
        {
          double time_step_accuracy_adv;
          double time_step_accuracy_diff;
          double time_step_accuracy_react;
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_solvers_"
                                         + name + ".dat");
          useful::read(input, time_step_accuracy_adv);
          useful::read(input, time_step_accuracy_diff);
          useful::read(input, time_step_accuracy_react);
          input.close();
          
          time_step = std::min({
            time_step_accuracy_adv*params_transport.advection_time,
            time_step_accuracy_diff*params_transport.diffusion_time,
            time_step_accuracy_react*params_reaction.reaction_time });
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Solver parameters\n"
            "--------------------------------------------------\n"
            "- Maximum timestep in units of advection time\n"
            "- Maximum timestep in units of diffusion time\n"
            "- Maximum timestep in units of reaction time\n"
            "--------------------------------------------------\n";
        }
      };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Solvers\n"
          "--------------------------------------------------\n"
          "Advection: RK4\n"
          "Diffusion: Stochastic Euler\n"
          "--------------------------------------------------\n";
      }
    };
    
    struct Transport
    {
      struct Parameters
      {
        double lengthscale;
        double peclet;
        double diff_coeff;
        
        double diffusion_time;
        double advection_time;
        
        Parameters
        (Directories const& directories, std::string const& name)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_transport_"
                                         + name + ".dat");
          useful::read(input, lengthscale);
          useful::read(input, peclet);
          useful::read(input, diff_coeff);
          input.close();
          
          diffusion_time = lengthscale*lengthscale/(2.*diff_coeff);
          advection_time = 2.*diffusion_time/peclet;
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Transport parameters\n"
            "--------------------------------------------------\n"
            "- Reference lengthscale\n"
            "- Peclet number\n"
            "- Diffusion coefficient\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Boundary>
      static auto makeTransitions
      (VelocityField&& velocity_field,
       Geometry const& geometry,
       Boundary const& boundary,
       Parameters const& params_transport,
       Solvers::Parameters const& params_solvers)
      {
        return makeTransportTransitions<
          Geometry, Solvers::Steppers>(velocity_field,
                                       geometry,
                                       boundary,
                                       params_transport,
                                       params_solvers);
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry, double average_velocity_magnitude)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data_rescaled(geometry.mesh,
                                           average_velocity_magnitude));
      };
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data(geometry.mesh));
      };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Transport\n"
          "--------------------------------------------------\n"
          "Processes: Advection-diffusion\n"
          "Interpolation: linear\n"
          "--------------------------------------------------\n";
      }
    };
    
    struct Reaction
    {
      struct Parameters
      {
        double damkohler{ 0. };
        double reaction_rate{ 0. };
        double reaction_time{ 1./reaction_rate };
        
        template <typename TransportParameters>
        Parameters
        (Directories const& dir,
         std::string const& name,
         TransportParameters const& params_transport)
        {}
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Reaction parameters\n"
            "--------------------------------------------------\n"
            "None\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename TransportParameters,
      typename SolverParameters>
      static auto makeReaction
      (Geometry const& geometry,
       Parameters const& params,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        return useful::DoNothing{};
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Reaction\n"
          "--------------------------------------------------\n"
          "None\n"
          "--------------------------------------------------\n";
      }
    };
    
    struct InitialCondition
    {
      struct Parameters
      {
        InitialConditions::Type type;
        double initial_mass;
        std::size_t nr_particles;
        double distance_wall;
        
        double time_step_accuracy_adv;
        double time_step_accuracy_diff;
        double time_step_accuracy_react;
        
        double time_min{ 0. };
        double time_max{ 1. };
        double time_step{ 1. };
        
        template
        <typename TransportParameters,
        typename ReactionParameters,
        typename SolverParameters>
        Parameters
        (Directories const& directories,
         std::string const& name,
         TransportParameters const& params_transport,
         ReactionParameters const& params_reaction,
         SolverParameters const& params_solvers)
        : distance_wall{
          10.*std::sqrt(2.*params_transport.diff_coeff*
                        params_solvers.time_step) }
        {
          std::string ic_name;
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_initial_condition_"
                                         + name + ".dat");
          useful::read(input, ic_name);
          verify_initial_condition(ic_name,
                                   InitialConditions{});
          type = InitialConditions::type(ic_name);
          useful::read(input, initial_mass);
          useful::read(input, nr_particles);
          input >> time_min;
          if (input.fail())
            return;
          useful::read(input, time_max);
          useful::read(input, time_step_accuracy_adv);
          useful::read(input, time_step_accuracy_diff);
          useful::read(input, time_step_accuracy_react);
          
          time_step =
            std::min({
              time_step_accuracy_adv*params_transport.advection_time,
              time_step_accuracy_diff*params_transport.diffusion_time,
              time_step_accuracy_react*params_reaction.reaction_time });
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Initial condition parameters (pass '' for default in [])\n"
            "--------------------------------------------------\n"
            "- Initial condition type:\n"
            "\tuniform: Homogeneous throughout the domain\n"
            "\tflux_weighted: Flux-weighted throughout the domain\n"
            "\tuniform_inlet: Homogeneous at the inlet\n"
            "\tflux_weighted_inlet: Flux-weighted at the inlet\n"
            "\tuniform_solid: Homogeneous at the solid surface\n"
            "\tuniform_near_solid: Homogeneous at a fixed distance to the solid interface\n"
            "\tuniform_inlet_continuous: Continuous injection homogeneous at the inlet\n"
            "\tflux_weighted_inlet_continuous: Continuous injection flux-weighted at the inlet\n"
            "- Total transported mass in each injection discretization\n"
            "- Number of Lagrangian particles in each injection discretization\n"
            "- Initial injection time [0.]\n"
            "- Final injection time [1.]\n"
            "- Maximum timestep in units of advection time [unit time step]\n"
            "- Maximum timestep in units of diffusion time [unit time step]\n"
            "- Maximum timestep in units of reaction time [unit time step]\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters params,
       double time = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            time,
            params.initial_mass/params.nr_particles },
          params };
      }
    };
    
    using VelocityField = VectorField_LinearInterpolation_OF
    <Foam::volVectorField, Locator_Cell const&, 1, 1>;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
  }

  namespace model_advection_3d
  {
    struct Model
    {
      inline static const std::string name{ "advection_3d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };
    
    using model_advection_diffusion_3d::Geometry;
    using model_advection_diffusion_3d::Info;
    using model_advection_diffusion_3d::State;
    using model_advection_diffusion_3d::CTRW;
    using model_advection_diffusion_3d::Reaction;
    using model_advection_diffusion_3d::InitialCondition;
    using model_advection_diffusion_3d::VelocityField;
    using model_advection_diffusion_3d::Output;
    
    struct Solvers
    {
      using Steppers = Steppers_Advection_RK4_Diffusion_Euler;
      
      struct Parameters
      {
        double time_step = 0.;
        
        std::size_t nr_particles;
        double step_length;
        
        template
        <typename TransportParameters,
        typename ReactionParameters>
        Parameters
        (Directories const& directories,
         std::string const& name,
         TransportParameters const& params_transport,
         ReactionParameters const& params_reaction)
        {
          double step_length_nondim;
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_solvers_"
                                         + name + ".dat");
          useful::read(input, step_length_nondim);
          input.close();
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Solver parameters\n"
            "--------------------------------------------------\n"
            "- Approximate step length in units of characteristic lengthscale\n"
            "--------------------------------------------------\n";
        }
      };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Solvers\n"
          "--------------------------------------------------\n"
          "Advection: RK4\n"
          "--------------------------------------------------\n";
      }
    };
    
    struct Transport
    {
      struct Parameters
      {
        double peclet = std::numeric_limits<double>::infinity();
        double diff_coeff = 0.;
        double diffusion_time = std::numeric_limits<double>::infinity();
        double advection_time = 0.;
        
        double lengthscale;
        
        Parameters(Directories const& directories, std::string const& name)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_transport_"
                                         + name + ".dat");
          useful::read(input, lengthscale);
          input.close();
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Transport parameters\n"
            "--------------------------------------------------\n"
            "- Reference lengthscale\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Boundary>
      static auto makeTransitions
      (VelocityField&& velocity_field,
       Geometry const& geometry,
       Boundary const& boundary,
       Parameters const& params_transport,
       Solvers::Parameters const& params_solvers)
      {
        return makeTransportTransitions_Advection<
          Geometry, Solvers::Steppers>(velocity_field,
                                       geometry,
                                       boundary,
                                       params_transport,
                                       params_solvers);
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry, double average_velocity_magnitude)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data_rescaled(geometry.mesh,
                                           average_velocity_magnitude));
      };
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data(geometry.mesh));
      };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Transport\n"
          "--------------------------------------------------\n"
          "Process: Advection\n"
          "Interpolation: linear\n"
          "--------------------------------------------------\n";
      }
    };
  };

  namespace model_advection_diffusion_fpt_3d
  {
    struct Model
    {
      inline static const std::string name{
        "advection_diffusion_fpt_3d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry<3,
      BoundaryConditionSet::Type::firstpassage>;
    using Info = ptof::Info_Absorbed_Reinjections;
    using State = StateDim<Geometry::dim,
      Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW<State>;
    using model_advection_diffusion_3d::Solvers;
    using model_advection_diffusion_3d::Transport;
    using model_advection_diffusion_3d::Reaction;
    using model_advection_diffusion_3d::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    struct InitialCondition final
    : model_advection_diffusion_3d::InitialCondition
    {
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters params,
       double time = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            time,
            params.initial_mass/params.nr_particles },
          params };
      }
    };
  }

  namespace model_advection_diffusion_decay_catalytic_3d
  {
    struct Model
    {
      inline static const std::string name{
        "advection_diffusion_decay_catalytic_3d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_advection_diffusion_3d::Geometry;
    using model_advection_diffusion_3d::Info;
    using model_advection_diffusion_3d::State;
    using model_advection_diffusion_3d::CTRW;
    using model_advection_diffusion_3d::Solvers;
    using model_advection_diffusion_3d::Transport;
    using model_advection_diffusion_3d::InitialCondition;
    using model_advection_diffusion_3d::Output;
    
    struct Reaction
    {
      struct Parameters
      {
        double damkohler;
        std::string initial_distribution;
        double surface_concentration;
        double reaction_rate;
        double reaction_time;
      
        template <typename TransportParameters>
        Parameters
        (Directories const& directories,
         std::string const& name,
         TransportParameters const& params_transport)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_reaction_"
                                         + name + ".dat");
          useful::read(input, damkohler);
          useful::read(input, initial_distribution);
          useful::read(input, surface_concentration);
          reaction_rate = params_transport.lengthscale*damkohler/
            (surface_concentration*params_transport.diffusion_time);
          reaction_time = 1./reaction_rate;
          input.close();
        }
      
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Reaction parameters:\n"
            "--------------------------------------------------\n"
            "- Damkohler number\n"
            "- Solid reactant initial distribution type\n"
            "\tuniform: Homogeneous throughout the domain\n"
            "- Initial solid reactant surface concentration\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename TransportParameters,
      typename SolverParameters>
      static auto makeReaction
      (Geometry const& geometry,
       Parameters const& params,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        if (params.initial_distribution == "uniform")
          return Reaction_AFluidPlusASolidtoASolid{
            params.reaction_rate,
            10.*std::sqrt(2.*params_transport.diff_coeff
                          *params_solvers.time_step),
            uniform_solid_reactant_patches(params.surface_concentration,
                                               { "wallFluidSolid" },
                                               geometry.mesh ),
            geometry.mesh_search
          };
        throw std::runtime_error{
          "Initial reactant distribution "
          + params.initial_distribution
          + " not supported"
        };
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        Reaction_AFluidPlusASolidtoASolid<Foam::meshSearch>::
          info(output);
      }
    };
  }

  namespace model_advection_diffusion_decay_3d
  {
    struct Model
    {
      inline static const std::string name{
        "advection_diffusion_decay_3d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_advection_diffusion_3d::Geometry;
    using model_advection_diffusion_3d::Info;
    using model_advection_diffusion_3d::State;
    using model_advection_diffusion_3d::CTRW;
    using model_advection_diffusion_3d::Solvers;
    using model_advection_diffusion_3d::Transport;
    using model_advection_diffusion_3d::InitialCondition;
    using model_advection_diffusion_3d::VelocityField;
    using model_advection_diffusion_3d::Output;
    
    struct Reaction
    {
      using Parameters =
        model_advection_diffusion_decay_catalytic_3d::
        Reaction::Parameters;
      
      template
      <typename Geometry,
      typename TransportParameters,
      typename SolverParameters>
      static auto makeReaction
      (Geometry const& geometry,
       Parameters const& params,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        return Reaction_AFluidPlusASolidtoNothing{
          params.reaction_rate,
          10.*std::sqrt(2.*params_transport.diff_coeff*
                        params_solvers.time_step),
          uniform_solid_reactant_patches(params.surface_concentration,
                                             { "wallFluidSolid" },
                                             geometry.mesh ),
          geometry.mesh_search
        };
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        Reaction_AFluidPlusASolidtoNothing<Foam::meshSearch>::
          info(output);
      }
    };
  }

  namespace model_periodic_cartesian_advection_diffusion_3d
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_3d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry_Periodic_Cartesian<3,
      BoundaryConditionSet::Type::transport>;
    using model_advection_diffusion_3d::Info;
    using State = StateDim_Periodic<
      Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW<State>;
    using model_advection_diffusion_3d::Solvers;
    using model_advection_diffusion_3d::Transport;
    using model_advection_diffusion_3d::Reaction;
    using model_advection_diffusion_3d::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    struct InitialCondition final
    : model_advection_diffusion_3d::InitialCondition
    {
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters params,
       double time = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker_Periodic{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            geometry.boundary_periodic,
            time,
            params.initial_mass/params.nr_particles },
          params };
      }
    };
  }

  namespace model_periodic_cartesian_advection_3d
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_3d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };
    
    using model_periodic_cartesian_advection_diffusion_3d::Geometry;
    using model_periodic_cartesian_advection_diffusion_3d::Info;
    using model_periodic_cartesian_advection_diffusion_3d::State;
    using model_periodic_cartesian_advection_diffusion_3d::CTRW;
    using model_advection_3d::Solvers;
    using model_advection_3d::Transport;
    using model_advection_3d::Reaction;
    using model_periodic_cartesian_advection_diffusion_3d::InitialCondition;
    using model_advection_3d::VelocityField;
    using model_periodic_cartesian_advection_diffusion_3d::Output;
  };

  namespace model_periodic_cartesian_advection_diffusion_fpt_3d
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_fpt_3d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry_Periodic_Cartesian<3,
      BoundaryConditionSet::Type::firstpassage>;
    using model_advection_diffusion_fpt_3d::Info;
    using model_periodic_cartesian_advection_diffusion_3d::State;
    using model_periodic_cartesian_advection_diffusion_3d::CTRW;
    using model_periodic_cartesian_advection_diffusion_3d::Solvers;
    using model_periodic_cartesian_advection_diffusion_3d::Transport;
    using model_periodic_cartesian_advection_diffusion_3d::Reaction;
    using model_periodic_cartesian_advection_diffusion_3d::InitialCondition;
    using model_periodic_cartesian_advection_diffusion_3d::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
  }

  namespace model_periodic_cartesian_advection_diffusion_decay_catalytic_3d
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_decay_catalytic_3d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_periodic_cartesian_advection_diffusion_3d::Geometry;
    using model_periodic_cartesian_advection_diffusion_3d::Info;
    using model_periodic_cartesian_advection_diffusion_3d::State;
    using model_periodic_cartesian_advection_diffusion_3d::CTRW;
    using model_periodic_cartesian_advection_diffusion_3d::Solvers;
    using model_periodic_cartesian_advection_diffusion_3d::Transport;
    using model_periodic_cartesian_advection_diffusion_3d::InitialCondition;
    using model_periodic_cartesian_advection_diffusion_3d::Output;
    using model_advection_diffusion_decay_catalytic_3d::Reaction;
  }

  namespace model_periodic_cartesian_advection_diffusion_decay_3d
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_decay_3d" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_periodic_cartesian_advection_diffusion_decay_catalytic_3d::Geometry;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_3d::Info;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_3d::State;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_3d::CTRW;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_3d::Solvers;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_3d::Transport;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_3d::InitialCondition;
    using model_periodic_cartesian_advection_diffusion_decay_catalytic_3d::Output;
    using model_advection_diffusion_decay_3d::Reaction;
  }
  
  namespace model_bcc_cartesian_advection_diffusion
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_cartesian_advection_diffusion" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry_Bcc<3,
      BoundaryConditionSet::Type::transport,
      BoundaryConditionSet::Type::cartesian>;
    using model_advection_diffusion_2d::Info;
    using State = StateDim_Periodic<Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW<State>;
    using model_advection_diffusion_2d::Solvers;
    using model_advection_diffusion_2d::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    using model_advection_diffusion_2d::Reaction;
    
    struct InitialCondition final
    : model_advection_diffusion_2d::InitialCondition
    {
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters params,
       double time = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker_Periodic{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            geometry.boundary_periodic,
            time,
            params.initial_mass/params.nr_particles },
          params };
      }
    };
    
    struct Transport
    {
      struct Parameters
      {
        double radius;
        double lengthscale;
        double peclet;
        double diff_coeff;
        
        double cell_side;
        std::vector<std::pair<double, double>> primitive_cell_boundaries;
        double diffusion_time;
        double advection_time;
        
        Parameters(Directories const& directories, std::string const& name)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_transport_"
                                         + name + ".dat");
          useful::read(input, radius);
          cell_side = 4./std::sqrt(3.)*radius;
          std::string lengthscale_definition;
          useful::read(input, lengthscale_definition);
          if (lengthscale_definition == "radius")
            lengthscale = radius;
          else if (lengthscale_definition == "diameter")
            lengthscale = 2.*radius;
          else if (lengthscale_definition == "cell_side")
            lengthscale = cell_side;
          else if (lengthscale_definition == "custom")
            useful::read(input, lengthscale);
          else
            throw std::runtime_error{
              std::string("Lengthscale definition")
              + lengthscale_definition
              + "not supported" };
          std::string primitive_cell_location;
          useful::read(input, primitive_cell_location);
          if (lengthscale_definition == "centered")
            primitive_cell_boundaries
              = std::vector<std::pair<double, double>>
                (3, { -cell_side/2., cell_side/2. });
          else if (lengthscale_definition == "corner")
            primitive_cell_boundaries
              = std::vector<std::pair<double, double>>
                (3, { 0., cell_side });
          else
            throw std::runtime_error{
              std::string("Lengthscale definition")
              + lengthscale_definition
              + "not supported" };
          useful::read(input, peclet);
          useful::read(input, diff_coeff);
          input.close();
          
          diffusion_time = lengthscale*lengthscale/(2.*diff_coeff);
          advection_time = 2.*diffusion_time/peclet;
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Transport parameters\n"
            "--------------------------------------------------\n"
            "- Bead radius\n"
            "- Reference length scale definition\n"
            "\tradius: bead radius\n"
            "\tdiameter: bead diameter\n"
            "\tcell_side: primitive cubic cell side\n"
            "\tcustom: custom value\n"
            "- Reference length scale value, if defined as custom\n"
            "- Location of primitive cell:\n"
            "\tcentered: Centered at the origin\n"
            "\tcorner: Left bottom corner at the origin\n"
            "- Peclet number\n"
            "- Diffusion coefficient\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Boundary>
      static auto makeTransitions
      (VelocityField&& velocity_field,
       Geometry const& geometry,
       Boundary const& boundary,
       Parameters const& params_transport,
       Solvers::Parameters const& params_solvers)
      {
        return makeTransportTransitions<
          Geometry, Solvers::Steppers>(velocity_field,
                                       geometry,
                                       boundary,
                                       params_transport,
                                       params_solvers);
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry, double average_velocity_magnitude)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data_rescaled(geometry.mesh,
                                           average_velocity_magnitude));
      };
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data(geometry.mesh));
      };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Transport\n"
          "--------------------------------------------------\n"
          "Process: Advection-diffusion\n"
          "Interpolation: linear\n"
          "--------------------------------------------------\n";
      }
    };
  }
  
  namespace model_bcc_cartesian_advection_diffusion_fpt
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_cartesian_advection_diffusion_fpt" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry_Bcc<3,
      BoundaryConditionSet::Type::firstpassage,
      BoundaryConditionSet::Type::cartesian>;
    using model_advection_diffusion_fpt_2d::Info;
    using model_bcc_cartesian_advection_diffusion::State;
    using model_bcc_cartesian_advection_diffusion::CTRW;
    using model_bcc_cartesian_advection_diffusion::Solvers;
    using model_bcc_cartesian_advection_diffusion::VelocityField;
    using model_bcc_cartesian_advection_diffusion::Reaction;
    using model_bcc_cartesian_advection_diffusion::Transport;
    using model_bcc_cartesian_advection_diffusion::InitialCondition;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
  }
  
  namespace model_bcc_cartesian_advection
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_cartesian_advection" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_bcc_cartesian_advection_diffusion::Geometry;
    using model_bcc_cartesian_advection_diffusion::Info;
    using model_bcc_cartesian_advection_diffusion::State;
    using model_bcc_cartesian_advection_diffusion::CTRW;
    using model_advection_2d::Solvers;
    using model_bcc_cartesian_advection_diffusion::InitialCondition;
    using model_advection_2d::VelocityField;
    using model_bcc_cartesian_advection_diffusion::Output;
    using model_advection_2d::Reaction;
    
    struct Transport
    {
      struct Parameters
      {
        double peclet = std::numeric_limits<double>::infinity();
        double diff_coeff = 0.;
        double diffusion_time = std::numeric_limits<double>::infinity();
        double advection_time = 0.;
        
        double radius;
        double lengthscale;
        
        double cell_side;
        std::vector<std::pair<double, double>> primitive_cell_boundaries;
        
        Parameters(Directories const& directories, std::string const& name)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_transport_"
                                         + name + ".dat");
          useful::read(input, radius);
          cell_side = 4./std::sqrt(3.)*radius;
          std::string lengthscale_definition;
          useful::read(input, lengthscale_definition);
          if (lengthscale_definition == "radius")
            lengthscale = radius;
          else if (lengthscale_definition == "diameter")
            lengthscale = 2.*radius;
          else if (lengthscale_definition == "cell_side")
            lengthscale = cell_side;
          else if (lengthscale_definition == "custom")
            useful::read(input, lengthscale);
          else
            throw std::runtime_error{
              std::string("Lengthscale definition")
              + lengthscale_definition
              + "not supported" };
          std::string primitive_cell_location;
          useful::read(input, primitive_cell_location);
          if (lengthscale_definition == "centered")
            primitive_cell_boundaries
              = std::vector<std::pair<double, double>>
                (3, { -cell_side/2., cell_side/2. });
          else if (lengthscale_definition == "corner")
            primitive_cell_boundaries
              = std::vector<std::pair<double, double>>
                (3, { 0., cell_side });
          else
            throw std::runtime_error{
              std::string("Lengthscale definition")
              + lengthscale_definition
              + "not supported" };
          input.close();
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
            "--------------------------------------------------\n"
            "Transport parameters\n"
            "--------------------------------------------------\n"
            "- Bead radius\n"
            "- Reference length scale definition\n"
            "\tradius: bead radius\n"
            "\tdiameter: bead diameter\n"
            "\tcell_side: primitive cubic cell side\n"
            "\tcustom: custom value\n"
            "- Reference length scale value, if defined as custom\n"
            "- Location of primitive cell:\n"
            "\tcentered: Centered at the origin\n"
            "\tcorner: Left bottom corner at the origin\n"
            "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Boundary>
      static auto makeTransitions
      (VelocityField&& velocity_field,
       Geometry const& geometry,
       Boundary const& boundary,
       Parameters const& params_transport,
       Solvers::Parameters const& params_solvers)
      {
        return makeTransportTransitions_Advection<
          Geometry, Solvers::Steppers>(velocity_field,
                                       geometry,
                                       boundary,
                                       params_transport,
                                       params_solvers);
      }
     
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry, double average_velocity_magnitude)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data_rescaled(geometry.mesh,
                                           average_velocity_magnitude));
      };
     
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry)
      {
        return makeLinearInterpolator(
          geometry,
          ptof::get_velocity_data(geometry.mesh));
      };
     
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Transport\n"
          "--------------------------------------------------\n"
          "Process: Advection"
          "Interpolation: linear\n"
          "--------------------------------------------------\n";
      }
    };
  }
  
  namespace model_bcc_cartesian_advection_diffusion_decay_catalytic
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_cartesian_advection_diffusion_decay_catalytic" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_bcc_cartesian_advection_diffusion::Geometry;
    using model_bcc_cartesian_advection_diffusion::Info;
    using model_bcc_cartesian_advection_diffusion::State;
    using model_bcc_cartesian_advection_diffusion::CTRW;
    using model_bcc_cartesian_advection_diffusion::Solvers;
    using model_bcc_cartesian_advection_diffusion::Transport;
    using model_bcc_cartesian_advection_diffusion::InitialCondition;
    using model_bcc_cartesian_advection_diffusion::VelocityField;
    using model_bcc_cartesian_advection_diffusion::Output;
    using model_advection_diffusion_decay_catalytic_2d::Reaction;
  }
  
  namespace model_bcc_cartesian_advection_diffusion_decay
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_cartesian_advection_diffusion_decay" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_bcc_cartesian_advection_diffusion::Geometry;
    using model_bcc_cartesian_advection_diffusion::Info;
    using model_bcc_cartesian_advection_diffusion::State;
    using model_bcc_cartesian_advection_diffusion::CTRW;
    using model_bcc_cartesian_advection_diffusion::Solvers;
    using model_bcc_cartesian_advection_diffusion::Transport;
    using model_bcc_cartesian_advection_diffusion::InitialCondition;
    using model_bcc_cartesian_advection_diffusion::VelocityField;
    using model_bcc_cartesian_advection_diffusion::Output;
    using model_advection_diffusion_decay_2d::Reaction;
  }
  
  namespace model_bcc_symmetryplanes_advection_diffusion
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry_Bcc<3,
      BoundaryConditionSet::Type::transport,
      BoundaryConditionSet::Type::symmetryplanes>;
    using model_bcc_cartesian_advection_diffusion::Info;
    using model_bcc_cartesian_advection_diffusion::State;
    using model_bcc_cartesian_advection_diffusion::CTRW;
    using model_bcc_cartesian_advection_diffusion::Solvers;
    using model_bcc_cartesian_advection_diffusion::Transport;
    using model_bcc_cartesian_advection_diffusion::InitialCondition;
    using model_bcc_cartesian_advection_diffusion::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    using model_bcc_cartesian_advection_diffusion::Reaction;
  }
  
  namespace model_bcc_symmetryplanes_advection_diffusion_fpt
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_cartesian_advection_diffusion_fpt" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry_Bcc<3,
      BoundaryConditionSet::Type::firstpassage,
      BoundaryConditionSet::Type::symmetryplanes>;
    using model_advection_diffusion_fpt_2d::Info;
    using model_bcc_symmetryplanes_advection_diffusion::State;
    using model_bcc_symmetryplanes_advection_diffusion::CTRW;
    using model_bcc_symmetryplanes_advection_diffusion::Solvers;
    using model_bcc_symmetryplanes_advection_diffusion::VelocityField;
    using model_bcc_symmetryplanes_advection_diffusion::Reaction;
    using model_bcc_symmetryplanes_advection_diffusion::Transport;
    using model_bcc_symmetryplanes_advection_diffusion::InitialCondition;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
  }
  
  namespace model_bcc_symmetryplanes_advection
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_symmetryplanes_advection" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_bcc_symmetryplanes_advection_diffusion::Geometry;
    using model_bcc_cartesian_advection::Info;
    using model_bcc_cartesian_advection::State;
    using model_bcc_cartesian_advection::CTRW;
    using model_bcc_cartesian_advection::Solvers;
    using model_bcc_cartesian_advection::Transport;
    using model_bcc_cartesian_advection::InitialCondition;
    using model_bcc_cartesian_advection::VelocityField;
    using model_bcc_symmetryplanes_advection_diffusion::Output;
    using model_bcc_cartesian_advection::Reaction;
  }
  
  namespace model_bcc_symmetryplanes_advection_diffusion_decay_catalytic
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion_decay_catalytic" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_bcc_symmetryplanes_advection_diffusion::Geometry;
    using model_bcc_symmetryplanes_advection_diffusion::Info;
    using model_bcc_symmetryplanes_advection_diffusion::State;
    using model_bcc_symmetryplanes_advection_diffusion::CTRW;
    using model_bcc_symmetryplanes_advection_diffusion::Solvers;
    using model_bcc_symmetryplanes_advection_diffusion::Transport;
    using model_bcc_symmetryplanes_advection_diffusion::InitialCondition;
    using model_bcc_symmetryplanes_advection_diffusion::VelocityField;
    using model_bcc_symmetryplanes_advection_diffusion::Output;
    using model_advection_diffusion_decay_catalytic_2d::Reaction;
  }
  
  namespace model_bcc_symmetryplanes_advection_diffusion_decay
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion_decay" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          + name + "\n"
          "--------------------------------------------------\n";
      }
    };

    using model_bcc_symmetryplanes_advection_diffusion::Geometry;
    using model_bcc_symmetryplanes_advection_diffusion::Info;
    using model_bcc_symmetryplanes_advection_diffusion::State;
    using model_bcc_symmetryplanes_advection_diffusion::CTRW;
    using model_bcc_symmetryplanes_advection_diffusion::Solvers;
    using model_bcc_symmetryplanes_advection_diffusion::Transport;
    using model_bcc_symmetryplanes_advection_diffusion::InitialCondition;
    using model_bcc_symmetryplanes_advection_diffusion::VelocityField;
    using model_bcc_symmetryplanes_advection_diffusion::Output;
    using model_advection_diffusion_decay_2d::Reaction;
  }
}

#endif /* Models_OF_h */
