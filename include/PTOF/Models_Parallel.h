/**
  \file PTOF/Models_Parallel.h
  \author Tomás Aquino
  \date 12/02/2024
 */

#ifndef PTOF_MODELS_PARALLEL_H
#define PTOF_MODELS_PARALLEL_H

#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <string>
#include <type_traits>
#include <boost/algorithm/string.hpp>
#include <meshSearch.H>
#include "CTRW/CTRW_Parallel.h"
#include "General/Useful.h"
#include "PTOF/Advection.h"
#include "PTOF/Geometry_Parallel.h"
#include "PTOF/InitialConditions.h"
#include "PTOF/Info.h"
#include "PTOF/Locator.h"
#include "PTOF/Output.h"
#include "PTOF/Reaction_Parallel.h"
#include "PTOF/State.h"
#include "PTOF/ParticleMaker.h"
#include "PTOF/Steppers.h"
#include "PTOF/Transitions.h"

namespace ptof
{
  /** \namespace ptof::model_advection_diffusion_2d_parallel
   Definitions for 2D advective--diffusive transport, in parallel. */
  namespace model_advection_diffusion_2d_parallel
  {
    struct Model
    {
      inline static const std::string name{ "advection_diffusion_2d_parallel" };
      
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
    
    using Geometry = Geometry_Parallel<2,
      BoundaryConditionSet::Type::transport>;
    using Info = ptof::Info_Absorbed;
    using State = State<Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW_Parallel<State>;
    
    struct Solvers
    {
      using Steppers = Steppers_Advection_Euler_Diffusion_Euler;
      
      struct Parameters
      {
        double local_time_step_adv;
        double local_time_step_diff;
        double local_time_step_react;
        double global_time_step_adv;
        double global_time_step_diff;
        double global_time_step_react;
        
        template
        <typename TransportParameters,
        typename ReactionParameters>
        Parameters
        (Directories const& directories,
         std::string const& name,
         TransportParameters const& params_transport,
         ReactionParameters const& params_reaction)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_solvers_"
                                         + name + ".dat");
          useful::read(input, local_time_step_adv);
          useful::read(input, local_time_step_diff);
          useful::read(input, local_time_step_react);
          useful::read(input, global_time_step_adv);
          useful::read(input, global_time_step_diff);
          useful::read(input, global_time_step_react);
          input.close();
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
          "--------------------------------------------------\n"
          "Solver parameters\n"
          "--------------------------------------------------\n"
          "- Local time step in terms of cell-based advection time\n"
          "- Local time step in terms of cell-based diffusion time\n"
          "- Local time step in terms of cell-based surface reaction time\n"
          "- Global time step in terms of characteristic advection time\n"
          "- Global time step in terms of characteristic diffusion time\n"
          "- Global time step in terms of surface reaction time\n"
          "  (Note: minimum between processes and maximum between local and global is used)\n"
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
        "Advection: Euler\n"
        "Diffusion: Stochastic Euler\n"
        "--------------------------------------------------\n";
      }
    };
    
    struct Transport
    {
      struct Parameters
      {
        double lengthscale;
        double diff_coeff;
        double diffusion_time;
        double advection_time;
        std::string rescale_velocity_field;
        double peclet;
        double mean_velocity;
        double velocity_rescaling_factor{ 1. };
        
        Parameters
        (Directories const& directories, std::string const& name)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_transport_"
                                         + name + ".dat");
          useful::read(input, lengthscale);
          useful::read(input, diff_coeff);
          diffusion_time = lengthscale*lengthscale/(2.*diff_coeff);
          useful::read(input, rescale_velocity_field);
          if (rescale_velocity_field == "rescale_to_peclet")
          {
            useful::read(input, peclet);
            advection_time = 2.*diffusion_time/peclet;
            mean_velocity = lengthscale/advection_time;
          }
          else if (rescale_velocity_field == "rescale_to_mean")
          {
            useful::read(input, mean_velocity);
            peclet = lengthscale*mean_velocity/diff_coeff;
            advection_time = lengthscale/mean_velocity;
          }
          else if (rescale_velocity_field == "rescale_to_advection_time")
          {
            useful::read(input, advection_time);
            peclet = 2.*diffusion_time/advection_time;
            mean_velocity = lengthscale/advection_time;
          }
          else if (rescale_velocity_field == "no_rescale")
          {}
          else
            throw std::runtime_error{ "Flow field rescaling option "
              + rescale_velocity_field
              + " not supported" };
          input.close();
        }
        
        template <typename VelocityField, typename Mesh>
        void rescale
        (VelocityField& velocity_field, Mesh const& mesh)
        {
          double current_mean = ptof::magnitude_of_average(velocity_field.field(), mesh);
          if (rescale_velocity_field == "rescale_to_peclet")
          {
            velocity_rescaling_factor = mean_velocity/current_mean;
            velocity_field.rescale(velocity_rescaling_factor);
          }
          else if (rescale_velocity_field == "rescale_to_mean")
          {
            velocity_rescaling_factor = mean_velocity/current_mean;
            velocity_field.rescale(velocity_rescaling_factor);
          }
          else if (rescale_velocity_field == "rescale_to_advection_time")
          {
            velocity_rescaling_factor = mean_velocity/current_mean;
            velocity_field.rescale(velocity_rescaling_factor);
          }
          else if (rescale_velocity_field == "no_rescale")
          {
            peclet = lengthscale*current_mean/diff_coeff;
            advection_time = lengthscale/current_mean;
            mean_velocity = lengthscale/advection_time;
          }
          else
            throw std::runtime_error{ "Flow field rescaling option "
              + rescale_velocity_field
              + " not supported" };
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
          "--------------------------------------------------\n"
          "Transport parameters\n"
          "--------------------------------------------------\n"
          "- Reference lengthscale\n"
          "- Diffusion coefficient\n"
          "- Whether and how to rescale velocity field:\n"
          "  (Note: Rescaling of the velocity field is not handled when transport parameters are set;\n"
          "         rescale method must be called)\n"
          "\tno_rescale: Do not rescale\n"
          "\t            (Note: rescale method should still be called to compute the Peclet number,\n"
          "\t                   advection time, and mean velocity)\n"
          "\trescale_to_peclet: Rescale according to imposed peclet number\n"
          "\trescale_to_mean: Rescale according to imposed mean flow velocity\n"
          "\trescale_to_advection_time: Rescale according to imposed advection time\n"
          "- Peclet number (pass only if rescaling with rescale_to_peclet)\n"
          "- Mean flow velocity (pass only if rescaling with rescale_to_mean)\n"
          "- Advection time (pass only if rescaling with rescale_to_advection_time)\n"
          "--`------------------------------------------------\n";
        }
      };
      
      template
      <typename VelocityField,
      typename Geometry,
      typename Boundary,
      typename ReactionParameters>
      static auto makeTransitions
      (VelocityField const& velocity_field,
       Geometry const& geometry,
       Boundary& boundary,
       Parameters const& params_transport,
       ReactionParameters const& params_reaction,
       Solvers::Parameters const& params_solvers)
      {
        return makeTransportTransitions<
        Geometry, Solvers::Steppers>(velocity_field,
                                     geometry,
                                     boundary,
                                     params_transport,
                                     params_reaction,
                                     params_solvers);
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry, double average_velocity_magnitude)
      {
        return makeLinearInterpolator(
                                      geometry,
                                      ptof::get_velocity_data_rescaled(geometry.mesh(),
                                                                       average_velocity_magnitude));
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry)
      {
        return makeLinearInterpolator(
                                      geometry,
                                      ptof::get_velocity_data(geometry.mesh()));
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
        "--------------------------------------------------\n"
        "Transport\n"
        "--------------------------------------------------\n"
        "Processes: Advection-diffusion\n"
        "Interpolation: Linear\n"
        "--------------------------------------------------\n";
      }
    };
    
    struct Reaction
    {
      using ReactionType = Reaction_DoNothing;
      
      struct Parameters
      {
        double damkohler{ 0. };
        double rate_constant{ 0. };
        double reaction_time{ std::numeric_limits<double>::infinity() };
        
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
       Parameters const& parameters,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        return Reaction_DoNothing{};
      }
      
      template
      <typename Geometry,
      typename TransportParameters,
      typename SolverParameters>
      static auto makeSurfaceReaction
      (Geometry const& geometry,
       Parameters const& parameters,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        return SurfaceReaction_DoNothing{};
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        Reaction_DoNothing::info(output);
        output << "\n";
        SurfaceReaction_DoNothing::info(output);
      }
    };
    
    struct InitialCondition
    {
      struct Parameters
      {
        InitialConditions::Type type;
        double distance_wall;
        std::vector<std::pair<double, double>> region_boundaries;
        std::string position_data;
        double initial_mass;
        std::size_t nr_particles;
        double time_min;
        double time_step_accuracy_adv;
        double time_step_accuracy_diff;
        double time_step_accuracy_react;
        double time_max;
        double time_step;
        
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
        {
          std::string ic_name;
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_initial_condition_"
                                         + name + ".dat");
          useful::read(input, ic_name);
          verify_initial_condition(ic_name,
                                   InitialConditions{});
          type = InitialConditions::type(ic_name);
          if (type == InitialConditions::Type::uniform_near_solid)
            useful::read(input, distance_wall);
          if (type == InitialConditions::Type::uniform_region_cartesian ||
              type == InitialConditions::Type::flux_weighted_region_cartesian)
          {
            input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::string line;
            std::getline(input, line);
            std::vector<std::string> split_line;
            boost::trim_if(line, boost::is_any_of("\t,|\r "));
            boost::algorithm::split(split_line, line, boost::is_any_of("\t,| "),
                                    boost::token_compress_on);
            if (split_line.size()%2 != 0)
              std::runtime_error{
                "Initial condition region boundaries must come in pairs" };
            for (std::size_t ii = 0; ii < split_line.size(); ii += 2)
            {
              region_boundaries.push_back( {
                std::stod(split_line[ii]), std::stod(split_line[ii+1]) } );
            }
          }
          if (type == InitialConditions::Type::prescribed_positions)
            useful::read(input, position_data);
          useful::read(input, initial_mass);
          useful::read(input, nr_particles);
          useful::read(input, time_min);
          if (type == InitialConditions::Type::uniform_inlet_continuous ||
              type == InitialConditions::Type::flux_weighted_inlet_continuous)
          {
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
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
          "--------------------------------------------------\n"
          "Initial condition parameters\n"
          "--------------------------------------------------\n"
          "- Initial condition type:\n"
          "\tuniform: Homogeneous throughout the domain\n"
          "\tflux_weighted: Flux-weighted throughout the domain\n"
          "\tuniform_inlet: Homogeneous at the inlet\n"
          "\tflux_weighted_inlet: Flux-weighted at the inlet\n"
          "\tuniform_solid: Homogeneous at the solid surface\n"
          "\tuniform_near_solid: Homogeneous at a fixed distance to the solid interface\n"
          "\tuniform_region_cartesian: Homogeneous in a prescribed cartesian region\n"
          "\tflux_weighted_region_cartesian: Flux-weighted in a prescribed cartesian region\n"
          "\tprescribed_positions: Prescribed positions\n"
          "\tuniform_inlet_continuous: Continuous injection homogeneous at the inlet\n"
          "\tflux_weighted_inlet_continuous: Continuous injection flux-weighted at the inlet\n"
          "- Initial distance from solid phase (with full path) for prescribed positions (pass only for type prescribed_positions)\n"
          "- Region boundaries (pass only for types uniform_region_cartesian or flux_weighted_region_cartesian)\n"
          "- Filename (with full path) for prescribed positions (pass only for type prescribed_positions)\n"
          "- Total transported mass in each injection step\n"
          "- Number of Lagrangian particles in each injection step\n"
          "- Initial injection time\n"
          "- Final injection time (pass only for types uniform_inlet_continuous or flux_weighted_inlet_continuous)\n"
          "- Maximum timestep for continuous injection discretization in units of advection time (pass only for types uniform_inlet_continuous or flux_weighted_inlet_continuous)\n"
          "- Maximum timestep for continuous injection discretization in units of diffusion time (pass only for types uniform_inlet_continuous or flux_weighted_inlet_continuous)\n"
          "- Maximum timestep for continuous injection discretization in units of reaction time (pass only for types uniform_inlet_continuous or flux_weighted_inlet_continuous)\n"
          "- Number of Lagrangian particles in each injection discretization\n"
          "--------------------------------------------------\n";
        }
      };
      
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters };
      }
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Mask>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters,
       Mask const& mask,
       double threshold = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters,
          mask, threshold };
      }
    };
    
    
    using VelocityField = VectorField_LinearInterpolation_OF
    <Foam::volVectorField, Locator_Cell_Parallel const&, 1, 1>;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    template
    <typename Transport,
    typename VelocityField,
    typename Geometry,
    typename Boundary,
    typename TransportParameters,
    typename ReactionParameters,
    typename SolverParameters,
    typename Reaction>
    struct TransitionMaker
    {
      TransitionMaker
      (VelocityField const& velocity_field,
       Geometry const& geometry,
       Boundary& boundary,
       TransportParameters const& params_transport,
       ReactionParameters const& params_reaction,
       SolverParameters const& params_solvers,
       Reaction const& reaction,
       std::size_t num_threads)
      : _velocity_field{ velocity_field }
      , _geometry{ geometry }
      , _boundary{ boundary }
      , _params_transport{ params_transport }
      , _params_reaction{ params_reaction }
      , _params_solvers{ params_solvers }
      , _reaction{ reaction }
      , _num_threads{ num_threads }
      {}
      
      TransitionMaker
      (VelocityField const& velocity_field,
       Geometry const& geometry,
       Boundary& boundary,
       TransportParameters const& params_transport,
       ReactionParameters const& params_reaction,
       SolverParameters const& params_solvers,
       Reaction const& reaction,
       std::size_t num_threads,
       useful::Selector_t<Transport>)
      : TransitionMaker(velocity_field,
                        geometry,
                        boundary,
                        params_transport,
                        params_reaction,
                        params_solvers,
                        reaction,
                        num_threads)
      {}
      
      auto operator()()
      {
        std::vector<Transitions> transitions;
        for (std::size_t thread = 0; thread < _num_threads; ++thread)
          transitions.emplace_back(Transport::makeTransitions(_velocity_field,
                                                              _geometry,
                                                              _boundary,
                                                              _params_transport,
                                                              _params_reaction,
                                                              _params_solvers),
                                   _reaction);
        
        return transitions;
      }
      
    private:
      VelocityField const& _velocity_field;
      Geometry const& _geometry;
      Boundary& _boundary;
      TransportParameters const& _params_transport;
      ReactionParameters const& _params_reaction;
      SolverParameters const& _params_solvers;
      Reaction const& _reaction;
      std::size_t _num_threads;
      
      using Transitions = decltype(ctrw::Transitions_CTRW_Transport_Reaction{
          Transport::makeTransitions(_velocity_field,
                                     _geometry, _boundary,
                                     _params_transport, _params_reaction, _params_solvers),
          _reaction });
    };
    template
    <typename Transport,
    typename VelocityField,
    typename Geometry,
    typename Boundary,
    typename TransportParameters,
    typename ReactionParameters,
    typename SolverParameters,
    typename Reaction>
    TransitionMaker
    (VelocityField const&,
     Geometry const&,
     Boundary&,
     TransportParameters const&,
     ReactionParameters const&,
     SolverParameters const&,
     Reaction const&,
     std::size_t,
     useful::Selector_t<Transport>) ->
    TransitionMaker<Transport, VelocityField, Geometry, Boundary,
    TransportParameters, ReactionParameters, SolverParameters,
    Reaction>;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return
      model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                             geometry,
                                                             boundary,
                                                             params_transport,
                                                             params_reaction,
                                                             params_solvers,
                                                             reaction,
                                                             num_threads,
                                                             useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_advection_2d_parallel
   Definitions for 2D advective transport, in parallel. */
  namespace model_advection_2d_parallel
  {
    struct Model
    {
      inline static const std::string name{ "advection_2d_parallel" };
      
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
    
    using model_advection_diffusion_2d_parallel::Geometry;
    using model_advection_diffusion_2d_parallel::Info;
    using model_advection_diffusion_2d_parallel::State;
    using model_advection_diffusion_2d_parallel::CTRW;
    using model_advection_diffusion_2d_parallel::Reaction;
    using model_advection_diffusion_2d_parallel::InitialCondition;
    using model_advection_diffusion_2d_parallel::VelocityField;
    using model_advection_diffusion_2d_parallel::Output;
    
    struct Solvers
    {
      using Steppers = Steppers_Advection_Euler_Diffusion_Euler;
      
      struct Parameters
      {
        double local_time_step_adv;
        double local_time_step_diff = std::numeric_limits<double>::infinity();
        double local_time_step_react = std::numeric_limits<double>::infinity();
        double global_time_step_adv;
        double global_time_step_diff = std::numeric_limits<double>::infinity();
        double global_time_step_react = std::numeric_limits<double>::infinity();
        
        template
        <typename TransportParameters,
        typename ReactionParameters>
        Parameters
        (Directories const& directories,
         std::string const& name,
         TransportParameters const& params_transport,
         ReactionParameters const& params_reaction)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_solvers_"
                                         + name + ".dat");
          useful::read(input, local_time_step_adv);
          useful::read(input, global_time_step_adv);
          input.close();
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
          "--------------------------------------------------\n"
          "Solver parameters\n"
          "--------------------------------------------------\n"
          "- Local timestep accuracy in terms of cell-based advection time\n"
          "- Global timestep accuracy in terms of characteristic advection time\n"
          "  (Note: Minimum between processes and maximum between local and global is used)\n"
          "  (Note: Initial values (e.g., of flow) are used for global quantities)\n"
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
        "Advection: Euler\n"
        "--------------------------------------------------\n";
      }
    };
    
    struct Transport
    {
      struct Parameters
      {
        double lengthscale;
        const double diff_coeff{ 0. };
        const double diffusion_time{ std::numeric_limits<double>::infinity() };
        double advection_time;
        std::string rescale_velocity_field;
        const double peclet{ std::numeric_limits<double>::infinity() };
        double mean_velocity;
        double velocity_rescaling_factor{ 1. };
        
        Parameters
        (Directories const& directories, std::string const& name)
        {
          auto input = useful::open_read(directories.dir_parameters
                                         + "/parameters_transport_"
                                         + name + ".dat");
          useful::read(input, lengthscale);
          useful::read(input, rescale_velocity_field);
          if (rescale_velocity_field == "rescale_to_mean")
          {
            useful::read(input, mean_velocity);
            advection_time = lengthscale/mean_velocity;
          }
          else if (rescale_velocity_field == "rescale_to_advection_time")
          {
            useful::read(input, advection_time);
            mean_velocity = lengthscale/advection_time;
          }
          else if (rescale_velocity_field == "no_rescale")
          {}
          else
            throw std::runtime_error{ "Flow field rescaling option "
              + rescale_velocity_field
              + " not supported" };
          input.close();
        }
        
        template <typename VelocityField, typename Mesh>
        void rescale
        (VelocityField& velocity_field, Mesh const& mesh)
        {
          double current_mean = ptof::magnitude_of_average(velocity_field.field(), mesh);
          if (rescale_velocity_field == "rescale_to_mean")
          {
            velocity_rescaling_factor = mean_velocity/current_mean;
            velocity_field.rescale(velocity_rescaling_factor);
          }
          else if (rescale_velocity_field == "rescale_to_advection_time")
          {
            velocity_rescaling_factor = mean_velocity/current_mean;
            velocity_field.rescale(velocity_rescaling_factor);
          }
          else if (rescale_velocity_field == "no_rescale")
          {
            advection_time = lengthscale/current_mean;
            mean_velocity = lengthscale/advection_time;
          }
          else
            throw std::runtime_error{ "Flow field rescaling option "
              + rescale_velocity_field
              + " not supported" };
        }
        
        template <typename OStream>
        static void info(OStream& output)
        {
          output <<
          "--------------------------------------------------\n"
          "Transport parameters\n"
          "--------------------------------------------------\n"
          "- Reference lengthscale\n"
          "- Diffusion coefficient\n"
          "- Whether and how to rescale velocity field:\n"
          "  (Note: Rescaling of the velocity field is not handled when transport parameters are set;\n"
          "         rescale method must be called)\n"
          "\tno_rescale: Do not rescale\n"
          "\t            (Note: rescale method should still be called to compute the Peclet number,\n"
          "\t                   advection time, and mean velocity)\n"
          "\trescale_to_mean: Rescale according to imposed mean flow velocity\n"
          "\trescale_to_advection_time: Rescale according to imposed advection time\n"
          "- Mean flow velocity (pass only if rescaling with rescale_to_mean)\n"
          "- Advection time (pass only if rescaling with rescale_to_advection_time)\n"
          "--------------------------------------------------\n";
        }
      };
      
      template
      <typename VelocityField,
      typename Geometry,
      typename Boundary,
      typename ReactionParameters>
      static auto makeTransitions
      (VelocityField const& velocity_field,
       Geometry const& geometry,
       Boundary& boundary,
       Parameters const& params_transport,
       ReactionParameters const& params_reaction,
       Solvers::Parameters const& params_solvers)
      {
        return makeTransportTransitions_Advection<
        Geometry, Solvers::Steppers>(velocity_field,
                                     geometry,
                                     boundary,
                                     params_transport,
                                     params_reaction,
                                     params_solvers);
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry, double average_velocity_magnitude)
      {
        return makeLinearInterpolator(
                                      geometry,
                                      ptof::get_velocity_data_rescaled(geometry.mesh(),
                                                                       average_velocity_magnitude));
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry)
      {
        return makeLinearInterpolator(
                                      geometry,
                                      ptof::get_velocity_data(geometry.mesh()));
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
        "--------------------------------------------------\n"
        "Transport\n"
        "--------------------------------------------------\n"
        "Process: Advection\n"
        "Interpolation: Linear\n"
        "--------------------------------------------------\n";
      }
    };
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  };
  
  /** \namespace ptof::model_advection_diffusion_fpt_2d_parallel
   Definitions for first-passage times under 2D advective--diffusive transport, in parallel. */
  namespace model_advection_diffusion_fpt_2d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "advection_diffusion_fpt_2d_parallel" };
      
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
    
    using Geometry = Geometry_Parallel<2,
      BoundaryConditionSet::Type::firstpassage>;
    using Info = ptof::Info_Absorbed_Reinjections;
    using State = State<Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW_Parallel<State>;
    using model_advection_diffusion_2d_parallel::Solvers;
    using model_advection_diffusion_2d_parallel::Transport;
    using model_advection_diffusion_2d_parallel::Reaction;
    using model_advection_diffusion_2d_parallel::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    struct InitialCondition final
    : model_advection_diffusion_2d_parallel::InitialCondition
    {
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters };
      }
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Mask>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters,
       Mask const& mask,
       double threshold = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters,
          mask, threshold };
      }
    };
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_advection_diffusion_decay_catalytic_2d_parallel
   Definitions for 2D advective--diffusive transport with surface reaction \f$ A_F + B_S \to B_S\f$, in parallel. */
  namespace model_advection_diffusion_decay_catalytic_2d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "advection_diffusion_decay_catalytic_2d_parallel" };
      
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
    
    using model_advection_diffusion_2d_parallel::Geometry;
    using model_advection_diffusion_2d_parallel::Info;
    using model_advection_diffusion_2d_parallel::State;
    using model_advection_diffusion_2d_parallel::CTRW;
    using model_advection_diffusion_2d_parallel::Solvers;
    using model_advection_diffusion_2d_parallel::Transport;
    using model_advection_diffusion_2d_parallel::InitialCondition;
    using model_advection_diffusion_2d_parallel::VelocityField;
    using model_advection_diffusion_2d_parallel::Output;
    
    struct Reaction
    {
      using ReactionType = Reaction_DoNothing;
    
      struct Parameters
      {
        double damkohler;
        std::string initial_distribution;
        double surface_concentration;
        double rate_constant;
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
          rate_constant = params_transport.lengthscale*damkohler/
          (surface_concentration*params_transport.diffusion_time);
          reaction_time = params_transport.diffusion_time/damkohler;
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
       Parameters const& parameters,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        return Reaction_DoNothing{};
      }
      
      template
      <typename Geometry,
      typename TransportParameters,
      typename SolverParameters>
      static auto makeSurfaceReaction
      (Geometry const& geometry,
       Parameters const& parameters,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        if (parameters.initial_distribution == "uniform")
          return SurfaceReaction_AFluidPlusASolidtoASolid{
            parameters.rate_constant,
            params_transport.diff_coeff,
            uniform_solid_reactant_patches(parameters.surface_concentration,
                                           { "wallFluidSolid" },
                                           geometry.mesh() ) };
        throw std::runtime_error{
          "Initial reactant distribution "
          + parameters.initial_distribution
          + " not supported"
        };
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        Reaction_DoNothing::info(output);
        output << "\n";
        SurfaceReaction_AFluidPlusASolidtoASolid::info(output);
      }
    };
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_periodic_cartesian_advection_diffusion_2d_parallel
   Definitions for 2D advective--diffusive transport with some periodic boundaries aligned with the Cartesian axes, in parallel. */
  namespace model_periodic_cartesian_advection_diffusion_2d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_2d_parallel" };
      
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
    
    using Geometry = Geometry_Periodic_Cartesian_Parallel<2,
      BoundaryConditionSet::Type::transport>;
    using model_advection_diffusion_2d_parallel::Info;
    using State = State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW_Parallel<State>;
    using model_advection_diffusion_2d_parallel::Solvers;
    using model_advection_diffusion_2d_parallel::Transport;
    using model_advection_diffusion_2d_parallel::Reaction;
    using model_advection_diffusion_2d_parallel::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    struct InitialCondition final
    : model_advection_diffusion_2d_parallel::InitialCondition
    {
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker_Periodic{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            geometry.boundary_periodic,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters };
      }
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Mask>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters,
       Mask const& mask,
       double threshold = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker_Periodic{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            geometry.boundary_periodic,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters,
          mask, threshold };
      }
    };
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_periodic_cartesian_advection_2d_parallel
   Definitions for 2D advective transport with some periodic boundaries aligned with the Cartesian axes. */
  namespace model_periodic_cartesian_advection_2d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_2d_parallel" };
      
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
    
    using model_periodic_cartesian_advection_diffusion_2d_parallel::Geometry;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::Info;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::State;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::CTRW;
    using model_advection_2d_parallel::Solvers;
    using model_advection_2d_parallel::Transport;
    using model_advection_2d_parallel::Reaction;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::InitialCondition;
    using model_advection_2d_parallel::VelocityField;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::Output;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  };
  
  /** \namespace ptof::model_periodic_cartesian_advection_diffusion_fpt_2d_parallel
   Definitions for first-passage times under 2D advective--diffusive transport with some periodic boundaries aligned with the Cartesian axes, in parallel. */
  namespace model_periodic_cartesian_advection_diffusion_fpt_2d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_fpt_2d_parallel" };
      
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
    
    using Geometry = Geometry_Periodic_Cartesian_Parallel<2,
      BoundaryConditionSet::Type::firstpassage>;
    using model_advection_diffusion_fpt_2d_parallel::Info;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::State;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::CTRW;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::Solvers;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::Transport;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::Reaction;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::InitialCondition;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_periodic_cartesian_advection_diffusion_decay_catalytic_2d_parallel
   Definitions for first-passage times under 2D advective--diffusive transport with surface reaction \f$ A_F + B_S \to B_S\f$, with some periodic boundaries aligned with the Cartesian axes, in parallel. */
  namespace model_periodic_cartesian_advection_diffusion_decay_catalytic_2d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_decay_catalytic_2d_parallel" };
      
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
    
    using model_periodic_cartesian_advection_diffusion_2d_parallel::Geometry;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::Info;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::State;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::CTRW;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::Solvers;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::Transport;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::InitialCondition;
    using model_advection_diffusion_2d_parallel::VelocityField;
    using model_periodic_cartesian_advection_diffusion_2d_parallel::Output;
    using model_advection_diffusion_decay_catalytic_2d_parallel::Reaction;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_advection_diffusion_3d_parallel
   Definitions for 2D advective--diffusive transport, in parallel. */
  namespace model_advection_diffusion_3d_parallel
  {
    struct Model
    {
      inline static const std::string name{ "advection_diffusion_3d_parallel" };
      
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
    
    using Geometry = Geometry_Parallel<3,
      BoundaryConditionSet::Type::transport>;
    using Info = ptof::Info_Absorbed;
    using State = State<Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW_Parallel<State>;
    using Solvers = model_advection_diffusion_2d_parallel::Solvers;
    
    struct Transport
    {
      using Parameters = model_advection_diffusion_2d_parallel::Transport::Parameters;
      
      template
      <typename VelocityField,
      typename Geometry,
      typename Boundary,
      typename ReactionParameters>
      static auto makeTransitions
      (VelocityField const& velocity_field,
       Geometry const& geometry,
       Boundary& boundary,
       Parameters const& params_transport,
       ReactionParameters const& params_reaction,
       Solvers::Parameters const& params_solvers)
      {
        return makeTransportTransitions<
          Geometry, Solvers::Steppers>(velocity_field,
                                       geometry,
                                       boundary,
                                       params_transport,
                                       params_reaction,
                                       params_solvers);
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry, double average_velocity_magnitude)
      {
        return makeLinearInterpolator(
                                      geometry,
                                      ptof::get_velocity_data_rescaled(geometry.mesh(),
                                                                       average_velocity_magnitude));
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry)
      {
        return makeLinearInterpolator(
                                      geometry,
                                      ptof::get_velocity_data(geometry.mesh()));
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
        "--------------------------------------------------\n"
        "Transport\n"
        "--------------------------------------------------\n"
        "Processes: Advection-diffusion\n"
        "Interpolation: Linear\n"
        "--------------------------------------------------\n";
      }
    };
    
    struct Reaction
    {
      using ReactionType = Reaction_DoNothing;
      
      struct Parameters
      {
        double damkohler{ 0. };
        double rate_constant{ 0. };
        double reaction_time{ std::numeric_limits<double>::infinity() };
        
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
       Parameters const& parameters,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        return Reaction_DoNothing{};
      }
      
      template
      <typename Geometry,
      typename TransportParameters,
      typename SolverParameters>
      static auto makeSurfaceReaction
      (Geometry const& geometry,
       Parameters const& parameters,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        return SurfaceReaction_DoNothing{};
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        Reaction_DoNothing::info(output);
        output << "\n";
        SurfaceReaction_DoNothing::info(output);
      }
    };
    
    struct InitialCondition final
    : model_advection_diffusion_2d_parallel::InitialCondition
    {
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters };
      }
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Mask>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters,
       Mask const& mask,
       double threshold = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters,
          mask, threshold };
      }
    };
    
    using VelocityField = VectorField_LinearInterpolation_OF
    <Foam::volVectorField, Locator_Cell_Parallel const&, 1, 1>;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_advection_3d_parallel
   Definitions for 3D advective--diffusive transport, in parallel. */
  namespace model_advection_3d_parallel
  {
    struct Model
    {
      inline static const std::string name{ "advection_3d_parallel" };
      
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
    
    using model_advection_diffusion_3d_parallel::Geometry;
    using model_advection_diffusion_3d_parallel::Info;
    using model_advection_diffusion_3d_parallel::State;
    using model_advection_diffusion_3d_parallel::CTRW;
    using model_advection_diffusion_3d_parallel::Reaction;
    using model_advection_diffusion_3d_parallel::InitialCondition;
    using model_advection_diffusion_3d_parallel::VelocityField;
    using model_advection_diffusion_3d_parallel::Output;
    using model_advection_2d_parallel::Solvers;
    using model_advection_2d_parallel::Transport;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  };
  
  /** \namespace ptof::model_advection_diffusion_fpt_3d_parallel
   Definitions for first-passage times under 3D advective--diffusive transport, in parallel. */
  namespace model_advection_diffusion_fpt_3d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "advection_diffusion_fpt_3d_parallel" };
      
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
    
    using Geometry = Geometry_Parallel<3,
      BoundaryConditionSet::Type::firstpassage>;
    using Info = ptof::Info_Absorbed_Reinjections;
    using State = State<Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW_Parallel<State>;
    using model_advection_diffusion_3d_parallel::Solvers;
    using model_advection_diffusion_3d_parallel::Transport;
    using model_advection_diffusion_3d_parallel::Reaction;
    using model_advection_diffusion_3d_parallel::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    struct InitialCondition final
    : model_advection_diffusion_2d_parallel::InitialCondition
    {
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters };
      }
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Mask>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters,
       Mask const& mask,
       double threshold = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters,
          mask, threshold };
      }
    };
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_advection_diffusion_decay_catalytic_3d_parallel
   Definitions for 3D advective--diffusive transport with surface reaction \f$ A_F + B_S \to B_S\f$, in parallel. */
  namespace model_advection_diffusion_decay_catalytic_3d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "advection_diffusion_decay_catalytic_3d_parallel" };
      
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
    
    using model_advection_diffusion_3d_parallel::Geometry;
    using model_advection_diffusion_3d_parallel::Info;
    using model_advection_diffusion_3d_parallel::State;
    using model_advection_diffusion_3d_parallel::CTRW;
    using model_advection_diffusion_3d_parallel::Solvers;
    using model_advection_diffusion_3d_parallel::Transport;
    using model_advection_diffusion_3d_parallel::InitialCondition;
    using model_advection_diffusion_2d_parallel::VelocityField;
    using model_advection_diffusion_3d_parallel::Output;
    
    struct Reaction
    {
      using ReactionType = Reaction_DoNothing;
      
      struct Parameters
      {
        double damkohler;
        std::string initial_distribution;
        double surface_concentration;
        double rate_constant;
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
          rate_constant = params_transport.lengthscale*damkohler/
          (surface_concentration*params_transport.diffusion_time);
          reaction_time = params_transport.diffusion_time/damkohler;
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
       Parameters const& parameters,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        return Reaction_DoNothing{};
      }
      
      template
      <typename Geometry,
      typename TransportParameters,
      typename SolverParameters>
      static auto makeSurfaceReaction
      (Geometry const& geometry,
       Parameters const& parameters,
       TransportParameters const& params_transport,
       SolverParameters const& params_solvers)
      {
        if (parameters.initial_distribution == "uniform")
          return SurfaceReaction_AFluidPlusASolidtoASolid{
            parameters.rate_constant,
            params_transport.diff_coeff,
            uniform_solid_reactant_patches(parameters.surface_concentration,
                                           { "wallFluidSolid" },
                                           geometry.mesh() ) };
        throw std::runtime_error{
          "Initial reactant distribution "
          + parameters.initial_distribution
          + " not supported" };
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        Reaction_DoNothing::info(output);
        output << "\n";
        SurfaceReaction_AFluidPlusASolidtoASolid::info(output);
      }
    };
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_periodic_cartesian_advection_diffusion_3d_parallel
   Definitions for 3D advective--diffusive transport with some periodic boundaries aligned with the Cartesian axes. */
  namespace model_periodic_cartesian_advection_diffusion_3d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_3d_parallels" };
      
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
    
    using Geometry = Geometry_Periodic_Cartesian_Parallel<3,
      BoundaryConditionSet::Type::transport>;
    using model_advection_diffusion_3d_parallel::Info;
    using State = State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW_Parallel<State>;
    using model_advection_diffusion_3d_parallel::Solvers;
    using model_advection_diffusion_3d_parallel::Transport;
    using model_advection_diffusion_3d_parallel::Reaction;
    using model_advection_diffusion_3d_parallel::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    struct InitialCondition final
    : model_advection_diffusion_2d_parallel::InitialCondition
    {
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker_Periodic{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            geometry.boundary_periodic,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters };
      }
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Mask>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters,
       Mask const& mask,
       double threshold = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker_Periodic{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            geometry.boundary_periodic,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters,
          mask, threshold };
      }
    };
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_periodic_cartesian_advection_3d_parallel
   Definitions for 3D advective transport with some periodic boundaries aligned with the Cartesian axes. */
  namespace model_periodic_cartesian_advection_3d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_3d_parallel" };
      
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
    
    using model_periodic_cartesian_advection_diffusion_3d_parallel::Geometry;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::Info;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::State;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::CTRW;
    using model_advection_3d_parallel::Solvers;
    using model_advection_3d_parallel::Transport;
    using model_advection_3d_parallel::Reaction;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::InitialCondition;
    using model_advection_3d_parallel::VelocityField;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::Output;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  };
  
  /** \namespace ptof::model_periodic_cartesian_advection_diffusion_fpt_3d_parallel
   Definitions for first-passage times under 3D advective--diffusive transport with some periodic boundaries aligned with the Cartesian axes, in parallel. */
  namespace model_periodic_cartesian_advection_diffusion_fpt_3d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_fpt_3d_parallel" };
      
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
    
    using Geometry = Geometry_Periodic_Cartesian_Parallel<3,
      BoundaryConditionSet::Type::firstpassage>;
    using model_advection_diffusion_fpt_3d_parallel::Info;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::State;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::CTRW;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::Solvers;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::Transport;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::Reaction;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::InitialCondition;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_periodic_cartesian_advection_diffusion_decay_catalytic_3d_parallel
   Definitions for 3D advective--diffusive transport with surface reaction \f$ A_F + B_S \to B_S\f$, with some periodic boundaries aligned with the Cartesian axes, in parallel. */
  namespace model_periodic_cartesian_advection_diffusion_decay_catalytic_3d_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "periodic_cartesian_advection_diffusion_decay_catalytic_3d_parallel" };
      
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
    
    using model_periodic_cartesian_advection_diffusion_3d_parallel::Geometry;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::Info;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::State;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::CTRW;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::Solvers;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::Transport;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::InitialCondition;
    using model_advection_diffusion_2d_parallel::VelocityField;
    using model_periodic_cartesian_advection_diffusion_3d_parallel::Output;
    using model_advection_diffusion_decay_catalytic_3d_parallel::Reaction;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
   /** \namespace ptof::model_bcc_cartesian_advection_diffusion_parallel
    Definitions for 3D advective--diffusive transport in a body centered cubic pack, based on the primitive unit cel, in parallel. */
  namespace model_bcc_cartesian_advection_diffusion_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_cartesian_advection_diffusion_parallel" };
      
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
    
    using Geometry = Geometry_Bcc_Parallel<3,
      BoundaryConditionSet::Type::transport,
      BoundaryConditionSet::Type::cartesian>;
    using model_advection_diffusion_2d_parallel::Info;
    using State = State_Periodic<Geometry::dim, Info, double, double, std::size_t>;
    using CTRW = ctrw::CTRW_Parallel<State>;
    using model_advection_diffusion_2d_parallel::Solvers;
    using model_advection_diffusion_2d_parallel::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    using model_advection_diffusion_2d_parallel::Reaction;
    
    struct InitialCondition final
    : model_advection_diffusion_2d_parallel::InitialCondition
    {
      template
      <typename Geometry,
      typename VelocityField>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker_Periodic{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            geometry.boundary_periodic,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters };
      }
      
      template
      <typename Geometry,
      typename VelocityField,
      typename Mask>
      static auto makeInitialCondition
      (Geometry const& geometry,
       VelocityField const& velocity_field,
       Parameters const& parameters,
       Mask const& mask,
       double threshold = 0.)
      {
        return InitialCondition_Cases{
          geometry, velocity_field,
          ParticleMaker_Periodic{
            useful::Selector_t<CTRW::Particle>{},
            geometry.locator,
            geometry.boundary_periodic,
            parameters.time_min,
            parameters.initial_mass/parameters.nr_particles },
          parameters,
          mask, threshold };
      }
    };
    
    struct Transport
    {
      struct Parameters
      {
        double radius;
        double lengthscale;
        double diff_coeff;
        double diffusion_time;
        std::string rescale_velocity_field;
        double peclet;
        double mean_velocity;
        double velocity_rescaling_factor{ 1. };
        
        double cell_side;
        std::vector<std::pair<double, double>> primitive_cell_boundaries;
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
              "Lengthscale definition"
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
            throw std::runtime_error{ "Lengthscale definition"
              + lengthscale_definition
              + "not supported" };
          useful::read(input, diff_coeff);
          diffusion_time = lengthscale*lengthscale/(2.*diff_coeff);
          useful::read(input, rescale_velocity_field);
          if (rescale_velocity_field == "rescale_to_peclet")
          {
            useful::read(input, peclet);
            advection_time = 2.*diffusion_time/peclet;
            mean_velocity = lengthscale/advection_time;
          }
          else if (rescale_velocity_field == "rescale_to_mean")
          {
            useful::read(input, mean_velocity);
            peclet = lengthscale*mean_velocity/diff_coeff;
            advection_time = lengthscale/mean_velocity;
          }
          else if (rescale_velocity_field == "rescale_to_advection_time")
          {
            useful::read(input, advection_time);
            peclet = 2.*diffusion_time/advection_time;
            mean_velocity = lengthscale/advection_time;
          }
          else if (rescale_velocity_field == "no_rescale")
          {}
          else
            throw std::runtime_error{ "Flow field rescaling option "
              + rescale_velocity_field
              + " not supported" };
          input.close();
        }
        
        template <typename VelocityField, typename Mesh>
        void rescale
        (VelocityField& velocity_field, Mesh const& mesh)
        {
          double current_mean = ptof::magnitude_of_average(velocity_field.field(), mesh);
          if (rescale_velocity_field == "rescale_to_peclet")
          {
            velocity_rescaling_factor = mean_velocity/current_mean;
            velocity_field.rescale(velocity_rescaling_factor);
          }
          else if (rescale_velocity_field == "rescale_to_mean")
          {
            velocity_rescaling_factor = mean_velocity/current_mean;
            velocity_field.rescale(velocity_rescaling_factor);
          }
          else if (rescale_velocity_field == "rescale_to_advection_time")
          {
            velocity_rescaling_factor = mean_velocity/current_mean;
            velocity_field.rescale(velocity_rescaling_factor);
          }
          else if (rescale_velocity_field == "no_rescale")
          {
            peclet = lengthscale*current_mean/diff_coeff;
            advection_time = lengthscale/current_mean;
            mean_velocity = lengthscale/advection_time;
          }
          else
            throw std::runtime_error{ "Flow field rescaling option "
              + rescale_velocity_field
              + " not supported" };
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
          "- Diffusion coefficient\n"
          "- Whether and how to rescale velocity field:\n"
          "  (Note: Rescaling of the velocity field is not handled when transport parameters are set;\n"
          "         rescale method must be called)\n"
          "\tno_rescale: Do not rescale\n"
          "\t            (Note: rescale method should still be called to compute the Peclet number,\n"
          "\t                   advection time, and mean velocity)\n"
          "\trescale_to_peclet: Rescale according to imposed peclet number\n"
          "\trescale_to_mean: Rescale according to imposed mean flow velocity\n"
          "\trescale_to_advection_time: Rescale according to imposed advection time\n"
          "- Peclet number (pass only if rescaling with rescale_to_peclet)\n"
          "- Mean flow velocity (pass only if rescaling with rescale_to_mean)\n"
          "- Advection time (pass only if rescaling with rescale_to_advection_time)\n"
          "--------------------------------------------------\n";
        }
      };
      
      template
      <typename VelocityField,
      typename Geometry,
      typename Boundary,
      typename ReactionParameters>
      static auto makeTransitions
      (VelocityField const& velocity_field,
       Geometry const& geometry,
       Boundary& boundary,
       Parameters const& params_transport,
       ReactionParameters const& params_reaction,
       Solvers::Parameters const& params_solvers)
      {
        return makeTransportTransitions<
        Geometry, Solvers::Steppers>(velocity_field,
                                     geometry,
                                     boundary,
                                     params_transport,
                                     params_reaction,
                                     params_solvers);
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry, double average_velocity_magnitude)
      {
        return makeLinearInterpolator(
                                      geometry,
                                      ptof::get_velocity_data_rescaled(geometry.mesh(),
                                                                       average_velocity_magnitude));
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry)
      {
        return makeLinearInterpolator(
                                      geometry,
                                      ptof::get_velocity_data(geometry.mesh()));
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
        "--------------------------------------------------\n"
        "Transport\n"
        "--------------------------------------------------\n"
        "Process: Advection-diffusion\n"
        "Interpolation: Linear\n"
        "--------------------------------------------------\n";
      }
    };
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_bcc_cartesian_advection_diffusion_fpt_parallel
   Definitions for advective--diffusive transport in a body centered cubic beadpack, based on the primitive unit cell, in parallel. */
  namespace model_bcc_cartesian_advection_diffusion_fpt_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_cartesian_advection_diffusion_fpt_parallel" };
      
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
    
    using Geometry = Geometry_Bcc_Parallel<3,
      BoundaryConditionSet::Type::firstpassage,
      BoundaryConditionSet::Type::cartesian>;
    using model_advection_diffusion_fpt_2d_parallel::Info;
    using model_bcc_cartesian_advection_diffusion_parallel::State;
    using model_bcc_cartesian_advection_diffusion_parallel::CTRW;
    using model_bcc_cartesian_advection_diffusion_parallel::Solvers;
    using model_bcc_cartesian_advection_diffusion_parallel::VelocityField;
    using model_bcc_cartesian_advection_diffusion_parallel::Reaction;
    using model_bcc_cartesian_advection_diffusion_parallel::Transport;
    using model_advection_diffusion_2d_parallel::VelocityField;
    using model_bcc_cartesian_advection_diffusion_parallel::InitialCondition;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_bcc_cartesian_advection_parallel
   Definitions for first-passage times under advective--diffusive transport in a body centered cubic beadpack, based on the primitive unit cell. */
  namespace model_bcc_cartesian_advection_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_cartesian_advectio_parallel" };
      
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
    
    using model_bcc_cartesian_advection_diffusion_parallel::Geometry;
    using model_bcc_cartesian_advection_diffusion_parallel::Info;
    using model_bcc_cartesian_advection_diffusion_parallel::State;
    using model_bcc_cartesian_advection_diffusion_parallel::CTRW;
    using model_advection_2d_parallel::Solvers;
    using model_bcc_cartesian_advection_diffusion_parallel::InitialCondition;
    using model_advection_2d_parallel::VelocityField;
    using model_bcc_cartesian_advection_diffusion_parallel::Output;
    using model_advection_2d_parallel::Reaction;
    
    struct Transport
    {
      struct Parameters
      {
        double radius;
        double lengthscale;
        const double diff_coeff{ 0. };
        const double diffusion_time { std::numeric_limits<double>::infinity() };
        std::string rescale_velocity_field;
        const double peclet{ std::numeric_limits<double>::infinity() };
        double mean_velocity;
        double velocity_rescaling_factor{ 1. };
        
        double cell_side;
        std::vector<std::pair<double, double>> primitive_cell_boundaries;
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
            throw std::runtime_error{ "Lengthscale definition"
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
            throw std::runtime_error{ "Lengthscale definition"
              + lengthscale_definition
              + "not supported" };
          useful::read(input, rescale_velocity_field);
          if (rescale_velocity_field == "rescale_to_mean")
          {
            useful::read(input, mean_velocity);
            advection_time = lengthscale/mean_velocity;
          }
          else if (rescale_velocity_field == "rescale_to_advection_time")
          {
            useful::read(input, advection_time);
            mean_velocity = lengthscale/advection_time;
          }
          else if (rescale_velocity_field == "no_rescale")
          {}
          else
            throw std::runtime_error{ "Flow field rescaling option "
              + rescale_velocity_field
              + " not supported" };
          input.close();
        }
        
        template <typename VelocityField, typename Mesh>
        void rescale
        (VelocityField& velocity_field, Mesh const& mesh)
        {
          double current_mean = ptof::magnitude_of_average(velocity_field.field(), mesh);
          if (rescale_velocity_field == "rescale_to_mean")
          {
            velocity_rescaling_factor = mean_velocity/current_mean;
            velocity_field.rescale(velocity_rescaling_factor);
          }
          else if (rescale_velocity_field == "rescale_to_advection_time")
          {
            velocity_rescaling_factor = mean_velocity/current_mean;
            velocity_field.rescale(velocity_rescaling_factor);
          }
          else if (rescale_velocity_field == "no_rescale")
          {
            advection_time = lengthscale/current_mean;
            mean_velocity = lengthscale/advection_time;
          }
          else
            throw std::runtime_error{ "Flow field rescaling option "
              + rescale_velocity_field
              + " not supported" };
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
          "- Whether and how to rescale velocity field:\n"
          "  (Note: Rescaling of the velocity field is not handled when transport parameters are set;\n"
          "         rescale method must be called)\n"
          "\tno_rescale: Do not rescale\n"
          "\t            (Note: rescale method should still be called to compute the Peclet number,\n"
          "\t                   advection time, and mean velocity)\n"
          "\trescale_to_mean: Rescale according to imposed mean flow velocity\n"
          "\trescale_to_advection_time: Rescale according to imposed advection time\n"
          "- Mean flow velocity (pass only if rescaling with rescale_to_mean)\n"
          "- Advection time (pass only if rescaling with rescale_to_advection_time)\n"
          "--------------------------------------------------\n";
        }
      };
      
      template
      <typename VelocityField,
      typename Geometry,
      typename Boundary,
      typename ReactionParameters>
      static auto makeTransitions
      (VelocityField const& velocity_field,
       Geometry const& geometry,
       Boundary& boundary,
       Parameters const& params_transport,
       ReactionParameters const& params_reaction,
       Solvers::Parameters const& params_solvers)
      {
        return makeTransportTransitions_Advection<
        Geometry, Solvers::Steppers>(velocity_field,
                                     geometry,
                                     boundary,
                                     params_transport,
                                     params_reaction,
                                     params_solvers);
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry, double average_velocity_magnitude)
      {
        return makeLinearInterpolator(
                                      geometry,
                                      ptof::get_velocity_data_rescaled(geometry.mesh(),
                                                                       average_velocity_magnitude));
      }
      
      template <typename Geometry>
      static auto makeVelocityInterpolator
      (Geometry const& geometry)
      {
        return makeLinearInterpolator(
                                      geometry,
                                      ptof::get_velocity_data(geometry.mesh()));
      }
      
      template <typename OStream>
      static void info(OStream& output)
      {
        output <<
        "--------------------------------------------------\n"
        "Transport\n"
        "--------------------------------------------------\n"
        "Process: Advection"
        "Interpolation: Linear\n"
        "--------------------------------------------------\n";
      }
    };
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_bcc_cartesian_advection_diffusion_decay_catalytic_parallel
   Definitions for advective transport in a body centered cubic beadpack, based on the primitive unit cell. */
  namespace model_bcc_cartesian_advection_diffusion_decay_catalytic_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_cartesian_advection_diffusion_decay_catalytic_parallel" };
      
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
    
    using model_bcc_cartesian_advection_diffusion_parallel::Geometry;
    using model_bcc_cartesian_advection_diffusion_parallel::Info;
    using model_bcc_cartesian_advection_diffusion_parallel::State;
    using model_bcc_cartesian_advection_diffusion_parallel::CTRW;
    using model_bcc_cartesian_advection_diffusion_parallel::Solvers;
    using model_bcc_cartesian_advection_diffusion_parallel::Transport;
    using model_bcc_cartesian_advection_diffusion_parallel::InitialCondition;
    using model_bcc_cartesian_advection_diffusion_parallel::VelocityField;
    using model_bcc_cartesian_advection_diffusion_parallel::Output;
    using model_advection_diffusion_decay_catalytic_2d_parallel::Reaction;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_bcc_symmetryplanes_advection_diffusion_parallel
   Definitions for advective--diffusive transport with surface reaction \f$ A_F + B_S \to B_S\f$ in a body centered cubic beadpack, based on the primitive unit cell, in parallel. */
  namespace model_bcc_symmetryplanes_advection_diffusion_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion_parallel" };
      
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
    
    using Geometry = Geometry_Bcc_Parallel<3,
      BoundaryConditionSet::Type::transport,
      BoundaryConditionSet::Type::symmetryplanes>;
    using model_bcc_cartesian_advection_diffusion_parallel::Info;
    using model_bcc_cartesian_advection_diffusion_parallel::State;
    using model_bcc_cartesian_advection_diffusion_parallel::CTRW;
    using model_bcc_cartesian_advection_diffusion_parallel::Solvers;
    using model_bcc_cartesian_advection_diffusion_parallel::Transport;
    using model_bcc_cartesian_advection_diffusion_parallel::InitialCondition;
    using model_bcc_cartesian_advection_diffusion_parallel::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    using model_bcc_cartesian_advection_diffusion_parallel::Reaction;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_bcc_symmetryplanes_advection_diffusion_fpt_parallel
   Definitions for first-passage times under advective--diffusive transport in a body centered cubic beadpack, based on the minimal periodic unit cel, in parallell. */
  namespace model_bcc_symmetryplanes_advection_diffusion_fpt_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_cartesian_advection_diffusion_fpt_parallel" };
      
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
    
    using Geometry = Geometry_Bcc_Parallel<3,
      BoundaryConditionSet::Type::firstpassage,
      BoundaryConditionSet::Type::symmetryplanes>;
    using model_advection_diffusion_fpt_2d_parallel::Info;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::State;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::CTRW;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::Solvers;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::VelocityField;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::Reaction;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::Transport;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::InitialCondition;
    using model_advection_diffusion_2d_parallel::VelocityField;
    using Output = ptof::Output_Cases<CTRW, VelocityField, Geometry>;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_bcc_symmetryplanes_advection_parallel
  Definitions for advective transport in a body centered cubic beadpack, based on the minimal periodic unit cell, in parallel. */
  namespace model_bcc_symmetryplanes_advection_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_symmetryplanes_advection_parallel" };
      
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
    
    using model_bcc_symmetryplanes_advection_diffusion_parallel::Geometry;
    using model_bcc_cartesian_advection_parallel::Info;
    using model_bcc_cartesian_advection_parallel::State;
    using model_bcc_cartesian_advection_parallel::CTRW;
    using model_bcc_cartesian_advection_parallel::Solvers;
    using model_bcc_cartesian_advection_parallel::Transport;
    using model_bcc_cartesian_advection_parallel::InitialCondition;
    using model_bcc_cartesian_advection_parallel::VelocityField;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::Output;
    using model_bcc_cartesian_advection_parallel::Reaction;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
  
  /** \namespace ptof::model_bcc_symmetryplanes_advection_diffusion_decay_catalytic
  Definitions for advective--diffusive transport with surface reaction \f$ A_F + B_S \to B_S\f$ in a body centered cubic beadpack, based on the minimal periodic unit cell, in parallel. */
  namespace model_bcc_symmetryplanes_advection_diffusion_decay_catalytic_parallel
  {
    struct Model
    {
      inline static const std::string name{
        "bcc_symmetryplanes_advection_diffusion_decay_catalytic_parallel" };
      
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
    
    using model_bcc_symmetryplanes_advection_diffusion_parallel::Geometry;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::Info;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::State;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::CTRW;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::Solvers;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::Transport;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::InitialCondition;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::VelocityField;
    using model_bcc_symmetryplanes_advection_diffusion_parallel::Output;
    using model_advection_diffusion_decay_catalytic_2d_parallel::Reaction;
    
    template <typename Boundary>
    auto makeTransitions
    (VelocityField const& velocity_field,
     Geometry const& geometry,
     Boundary& boundary,
     Transport::Parameters const& params_transport,
     Reaction::Parameters const& params_reaction,
     Solvers::Parameters const& params_solvers,
     Reaction::ReactionType const& reaction,
     std::size_t num_threads)
    {
      return model_advection_diffusion_2d_parallel::TransitionMaker(velocity_field,
                                                                    geometry,
                                                                    boundary,
                                                                    params_transport,
                                                                    params_reaction,
                                                                    params_solvers,
                                                                    reaction,
                                                                    num_threads,
                                                                    useful::Selector_t<Transport>{})();
    }
  }
}

#endif /* PTOF_MODELS_PARALLEL_H */
