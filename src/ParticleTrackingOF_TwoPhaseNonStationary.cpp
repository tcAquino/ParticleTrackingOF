//
//  ParticleTrackingOF_TwoPhaseNonStationary.cpp
//  PTOF
//
//  Created by Tomás Aquino on 07/02/2024.
//

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <fvcGrad.H>
#include "General/Useful.h"
#include "PTOF/Directories.h"
#include "PTOF/Models.h"
#include "PTOF/Phase.h"

int main(int argc, char * argv[])
{
  using namespace ptof::model_periodic_cartesian_advection_diffusion_2d;
  using Phase = ptof::Phase;
  
  if (useful::check_options_help(argc, argv))
  {
    std::cout <<
      "--------------------------------------------------\n"
      "ParticleTrackingOF_TwoPhaseNonStationary\n"
      "--------------------------------------------------\n";
    std::cout << std::endl;
    std::cout <<
      "--------------------------------------------------\n"
      "Executable parameters (pass '' for default in []):\n"
      "--------------------------------------------------\n"
      "- Cases directory [../cases]\n"
      "- Name of case\n"
      "- Name of transport parameter set\n"
      "- Name of phase parameter set\n"
      "- Name of reaction parameter set\n"
      "- Name of solver parameter set\n"
      "- Name of initial condition parameter set\n"
      "- Name of output parameter set\n"
      "- Output directory [<Case directory>/output]\n"
      "--------------------------------------------------\n";
    std::cout << std::endl;
    Model::info(std::cout);
    std::cout << std::endl;
    Geometry::info(std::cout);
    std::cout << std::endl;
    ptof::DirectoriesOF::info(std::cout);
    std::cout << std::endl;
    Transport::info(std::cout);
    std::cout << std::endl;
    Reaction::info(std::cout);
    std::cout << std::endl;
    Solvers::info(std::cout);
    std::cout << std::endl;
    Transport::Parameters::info(std::cout);
    std::cout << std::endl;
    Phase::Parameters::info(std::cout);
    std::cout << std::endl;
    Reaction::Parameters::info(std::cout);
    std::cout << std::endl;
    Solvers::Parameters::info(std::cout);
    std::cout << std::endl;
    InitialCondition::Parameters::info(std::cout);
    std::cout << std::endl;
    Output::Parameters::info(std::cout);
    std::cout << std::endl;
    return 0;
  }
  if (argc != 10)
    throw useful::bad_parameters_help();
  
  std::size_t arg = 1;
  std::string dir = argv[arg++];
  std::string case_name = argv[arg++];
  std::string parameters_transport_name = argv[arg++];
  std::string parameters_phase_name = argv[arg++];
  std::string parameters_reaction_name = argv[arg++];
  std::string parameters_solvers_name = argv[arg++];
  std::string parameters_initial_condition_name = argv[arg++];
  std::string parameters_output_name = argv[arg++];
  std::string dir_output = argv[arg++];
  
  ptof::Directories directories{ dir, case_name, dir_output };
  directories.info_runtime(std::cout);
  std::cout << std::endl;
  ptof::DirectoriesOF directories_of{ directories };
  directories_of.info_runtime(std::cout);
  
  std::cout << std::setprecision(2) << std::scientific;
  
  std::cout << "\n" << "Importing transport parameters..." << std::endl;
  auto execution_begin = std::chrono::high_resolution_clock::now();
  Transport::Parameters params_transport{ directories,
    parameters_transport_name };
  auto execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Importing phase parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Phase::Parameters params_phase{ directories,
    parameters_phase_name };
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Importing reaction parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Reaction::Parameters params_reaction{ directories,
    parameters_reaction_name,
    params_transport };
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n"  << "Importing solver parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Solvers::Parameters params_solvers{ directories,
    parameters_solvers_name, params_transport, params_reaction };
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n"  << "Importing initial condition parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  InitialCondition::Parameters params_initial_condition{ directories,
    parameters_initial_condition_name,
    params_transport, params_reaction, params_solvers,
    
  };
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Setting up geometry...\n";
  execution_begin = std::chrono::high_resolution_clock::now();
  Geometry geometry{
    directories_of, directories, params_transport };
  geometry.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Setting up phases..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  auto excluded_phase_field = Phase::get_excluded_phase_data(geometry.mesh(),
                                                             params_phase);
  auto grad_excluded_phase_field = Phase::grad_excluded_phase(geometry.mesh(),
                                                              params_phase,
                                                              excluded_phase_field);
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Setting up velocity interpolation..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  auto velocity_field = Transport::makeVelocityInterpolator(geometry);
  params_transport.rescale(velocity_field, geometry.mesh());
  velocity_field.sum(Foam::dimensionedScalar("",
                                             Foam::dimensionSet(0, 2, -1, 0, 0, 0, 0),
                                             params_transport.diff_coeff)
                     * params_phase.leakage_coefficient
                     * grad_excluded_phase_field);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Setting up reaction..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  auto reaction = Reaction::makeReaction(geometry,
                                         params_reaction,
                                         params_transport,
                                         params_solvers);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Setting up initial condition...\n";
  execution_begin = std::chrono::high_resolution_clock::now();
  auto initial_condition
    = InitialCondition::makeInitialCondition(geometry,
                                             velocity_field,
                                             params_initial_condition,
                                             excluded_phase_field,
                                             params_phase.phase_threshold);
  initial_condition.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Setting up boundary conditions...\n";
  execution_begin = std::chrono::high_resolution_clock::now();
  auto boundary = geometry.makeBoundary(directories,
                                        params_transport,
                                        params_reaction,
                                        params_solvers,
                                        Reaction::makeSurfaceReaction(
                                          geometry,
                                          params_reaction,
                                          params_transport,
                                          params_solvers),
                                        initial_condition);
  boundary.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Setting up dynamics..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  CTRW ctrw{ initial_condition(), CTRW::Tag{} };
  ctrw::Transitions_CTRW_Transport_Reaction transitions{
    Transport::makeTransitions(velocity_field,
                               geometry, boundary,
                               params_transport, params_reaction, params_solvers),
    reaction };
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n"  << "Importing output parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Output::Parameters params_output{ directories,
    parameters_output_name,
    params_transport, params_reaction, params_solvers,
    params_transport.velocity_rescaling_factor };
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Setting up output..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Output measurer{
    ctrw,
    velocity_field,
    geometry,
    directories,
    params_output,
    std::string("M_") + "TwoPhaseNonStationary_" + Model::name
      + "_C_" + case_name
      + "_OF_" + directories_of.case_name
      + "_T_" + parameters_transport_name
      + "_P_" + parameters_phase_name
      + "_R_" + parameters_reaction_name
      + "_S_" + parameters_solvers_name
      + "_I_" + parameters_initial_condition_name
      + "_O_" + parameters_output_name,
    std::vector<Phase::PhaseField const*>{ &excluded_phase_field },
    { 1. - params_phase.phase_threshold }
  };
  measurer.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n" << "Starting dynamics..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  double current_time = directories_of.time.value();
  ptof::info_time(std::cout, params_output, current_time);
  auto const& flow_times = directories_of.time.times();
  auto closest_time_index = [&flow_times](auto value)
  {
    return std::distance(flow_times.cbegin(),
                         std::min_element(flow_times.cbegin(), flow_times.cend(),
                                          [value](auto xx, auto yy)
                                          {
      if (xx.name() == "constant" || yy.name() == "constant")
        return true;
      else
        return std::abs(xx.value() - value) < std::abs(yy.value() - value); }));
  };
  Foam::label time_index = closest_time_index(current_time);
  double next_flow_time = std::numeric_limits<double>::infinity();
  if (time_index < flow_times.size()-1)
  {
    next_flow_time = flow_times[time_index+1].value();
  }
  
  while (!measurer.done(current_time))
  {
    if (time_index < flow_times.size()-1)
    {
      bool velocity_field_needs_update = false;
      while (time_index < flow_times.size()-1
             && current_time >= next_flow_time)
      {
        ++time_index;
        next_flow_time = directories_of.time.times()[time_index].value();
        directories_of.time.setTime(next_flow_time, 0);
        velocity_field_needs_update = true;
      }
      if (velocity_field_needs_update)
      {
        std::cout << "Field updates required...\n";
        ptof::info_time(std::cout, params_output, current_time);
        std::cout << "Updating phase field...\n";
        excluded_phase_field = Phase::get_excluded_phase_data(geometry.mesh(),
                                                              params_phase);
        grad_excluded_phase_field = Phase::grad_excluded_phase(geometry.mesh(),
                                                               params_phase,
                                                               excluded_phase_field);
        std::cout << "Done!\n";
        std::cout << "Updating velocity field...\n";
        velocity_field.set(ptof::get_velocity_data(geometry.mesh()));
        if (params_transport.velocity_rescaling_factor != 1.)
          velocity_field.rescale(params_transport.velocity_rescaling_factor);
        velocity_field.sum(Foam::dimensionedScalar("",
                                                   Foam::dimensionSet(0, 2, -1, 0, 0, 0, 0),
                                                   params_transport.diff_coeff)
                           * params_phase.leakage_coefficient
                           * grad_excluded_phase_field);
        std::cout << "Done!\n";
        std::cout << "Done!\n";
      }
    }
    while (measurer.next_measure_time() <= current_time)
    {
      std::cout << "Measurement required...\n";
      ptof::info_time(std::cout, params_output, measurer.next_measure_time());
      ptof::info_fraction_not_absorbed(std::cout, ctrw, measurer.next_measure_time());
      measurer(measurer.next_measure_time());
      std::cout << "Done!\n";
    }
    current_time = next_flow_time > current_time
      ? std::min(measurer.next_measure_time(), next_flow_time)
      : measurer.next_measure_time();
    ctrw.evolve([current_time](CTRW::Particle const& part)
                { return part.state_new().time < current_time
                    && !part.state_new().info.absorbed; },
                transitions);
  }
  if (measurer.next_measure_time() <= current_time)
  {
    std::cout << "Measurement required...\n";
    ptof::info_time(std::cout, params_output, measurer.next_measure_time());
    ptof::info_fraction_not_absorbed(std::cout, ctrw, measurer.next_measure_time());
    measurer(measurer.next_measure_time());
    std::cout << "Done!\n";
  }
  measurer();
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  return 0;
}
