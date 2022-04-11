//
//  Models.h
//  ParticleTracking_OpenFOAM
//
//  Created by Tomás Aquino on 22/02/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Models_PT_OF_h
#define Models_PT_OF_h

#include <cmath>
#include <cstddef>
#include <exception>
#include <string>
#include <type_traits>
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
  namespace model_advection_diffusion
  {
    struct Model
    {
      inline static const std::string name{ "advection_diffusion" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          "Advective-diffusive particle tracking\n"
          "with RK4 stepping for advection\n"
          "and stochastic Euler stepping for diffusion\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry<
      2,
      ImplementedBoundaryConditionSets::Type::transport>;
    using Info = ptof::Info_Absorbed;
    using State = StateDim<Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW<State>;
    
    struct Solvers
    {
      using Steppers = Steppers_Advection_RK4_Diffusion_Euler;
      
      struct Parameters
      {
        std::size_t nr_particles;
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
        
        Parameters(Directories const& directories, std::string const& name)
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
    };
    
    struct InitialCondition
    {
      struct Parameters
      {
        ImplementedInitialConditions::Type type;
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
                                   ImplementedInitialConditions{});
          type = ImplementedInitialConditions::type(ic_name);
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
            geometry.locator, time,
            params.initial_mass/params.nr_particles },
          params };
      }
    };
    
    using Output = ptof::Output_Cases;
  }
  
  namespace model_advection_diffusion_fpt
  {
    struct Model
    {
      inline static const std::string name{ "advection_diffusion_fpt" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          "Advective-diffusive particle tracking\n"
          "for determination of first passage times to solid phase\n"
          "with RK4 stepping for advection\n"
          "and stochastic Euler stepping for diffusion\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = Geometry<
      2,
      ImplementedBoundaryConditionSets::Type::firstpassage>;
    using Info = ptof::Info_Absorbed_Reinjections;
    using State = StateDim<Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW<State>;
    using model_advection_diffusion::Solvers;
    using model_advection_diffusion::Transport;
    using model_advection_diffusion::Reaction;
    
    struct InitialCondition : model_advection_diffusion::InitialCondition
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
            geometry.locator, time,
            params.initial_mass/params.nr_particles },
          params };
      }
    };
    
    using Output = ptof::Output_Cases;
  }
  
  namespace model_advection_diffusion_decay_catalytic
  {
    struct Model
    {
      inline static const std::string name{ "advection_diffusion_decay_catalytic" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          "Advective-diffusive particle tracking\n"
          "with A_Fluid + A_Solid -> A_Solid surface reaction,\n"
          "with RK4 stepping for advection\n"
          "and stochastic Euler stepping for diffusion\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = model_advection_diffusion::Geometry;
    using Info = ptof::Info_Absorbed_Reinjections;
    using model_advection_diffusion::State;
    using model_advection_diffusion::CTRW;
    using model_advection_diffusion::Solvers;
    using model_advection_diffusion::Transport;
    using model_advection_diffusion::InitialCondition;
    using Output = ptof::Output_Cases;
    
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
            10.*std::sqrt(2.*params_transport.diff_coeff*params_solvers.time_step),
            uniform_solid_reactant_patches(params.surface_concentration,
                                               { "wallFluidSolid" },
                                               geometry.mesh ),
            geometry.mesh_search
          };
        throw std::runtime_error{
          std::string("Initial reactant distribution ")
          + params.initial_distribution
          + " not supported"
        };
      }
    };
  }
  
  namespace model_advection_diffusion_decay
  {
    struct Model
    {
      inline static const std::string name{ "advection_diffusion_decay" };
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
          "--------------------------------------------------\n"
          "Model\n"
          "--------------------------------------------------\n"
          "Advective-diffusive particle tracking\n"
          "with A_Fluid + A_Solid -> Nothing surface reaction,\n"
          "with RK4 stepping for advection\n"
          "and stochastic Euler stepping for diffusion\n"
          "--------------------------------------------------\n";
      }
    };

    using Geometry = model_advection_diffusion::Geometry;
    using Info = model_advection_diffusion::Info;
    using model_advection_diffusion::State;
    using model_advection_diffusion::CTRW;
    using model_advection_diffusion::Solvers;
    using model_advection_diffusion::Transport;
    using model_advection_diffusion::InitialCondition;
    using Output = ptof::Output_Cases;
    
    struct Reaction
    {
      using Parameters
        = model_advection_diffusion_decay_catalytic::Reaction::Parameters;
      
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
    };
  }
}


#endif /* Models_PT_OF_h */
