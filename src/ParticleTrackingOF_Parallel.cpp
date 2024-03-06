//
//  ParticleTrackingOF.cpp
//  PTOF
//
//  Created by Tomás Aquino on 16/02/2022.
//

#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <limits>
#include <string>
#include "General/Useful.h"
#include "PTOF/Directories.h"
#include "PTOF/Models_Parallel.h"

int main(int argc, char * argv[])
{
  using namespace ptof::model_advection_diffusion_decay_catalytic_2d_parallel;
  
  if (useful::check_options_help(argc, argv))
  {
    std::cout <<
      "--------------------------------------------------\n"
      "ParticleTrackingOF_Parallel\n"
      "--------------------------------------------------\n";
    std::cout << std::endl;
    std::cout <<
      "--------------------------------------------------\n"
      "Executable parameters (pass '' for default in []):\n"
      "--------------------------------------------------\n"
      "- Cases directory [../cases]\n"
      "- Name of case\n"
      "- Name of transport parameter set\n"
      "- Name of reaction parameter set\n"
      "- Name of solver parameter set\n"
      "- Name of initial condition parameter set\n"
      "- Name of output parameter set\n"
      "- Output directory [<Case directory>/output]\n"
      "- Number of parallel threads\n"
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
  std::string parameters_reaction_name = argv[arg++];
  std::string parameters_solvers_name = argv[arg++];
  std::string parameters_initial_condition_name = argv[arg++];
  std::string parameters_output_name = argv[arg++];
  std::string dir_output = argv[arg++];
  std::size_t num_threads = strtoul(argv[arg++], NULL, 0);
  
  if (num_threads < 1)
    throw useful::bad_parameters_help();
  omp_set_num_threads(int(num_threads));
  
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
    params_transport, params_reaction, params_solvers };
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Setting up geometry...\n";
  execution_begin = std::chrono::high_resolution_clock::now();
  Geometry geometry{
    directories_of, directories, num_threads, params_transport };
  geometry.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Setting up velocity interpolation..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  auto velocity_field = Transport::makeVelocityInterpolator(geometry);
  params_transport.rescale(velocity_field, geometry.mesh());
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
                                             params_initial_condition);
  initial_condition.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  
  std::cout << "\n" << "Setting up boundary conditions...\n";
  execution_begin = std::chrono::high_resolution_clock::now();
  auto surface_reaction = Reaction::makeSurfaceReaction(geometry,
                                                        params_reaction,
                                                        params_transport,
                                                        params_solvers);
  auto boundary = geometry.makeBoundary(directories,
                                        params_transport,
                                        params_reaction,
                                        params_solvers,
                                        surface_reaction,
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
  auto transitions = makeTransitions(velocity_field,
                                     geometry, boundary,
                                     params_transport, params_reaction, params_solvers,
                                     reaction,
                                     num_threads);
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
    "M_" + Model::name
      + "_C_" + case_name
      + "_OF_" + directories_of.case_name
      + "_T_" + parameters_transport_name
      + "_R_" + parameters_reaction_name
      + "_S_" + parameters_solvers_name
      + "_I_" + parameters_initial_condition_name
      + "_O_" + parameters_output_name };
  measurer.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n" << "Starting dynamics..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  double current_time = 0.;
  ptof::info_time(std::cout, params_output, current_time);
  while (!measurer.done(current_time))
  {
    while (measurer.next_measure_time() <= current_time)
    {
      std::cout << "Measurement required...\n";
      ptof::info_time(std::cout, params_output, measurer.next_measure_time());
      ptof::info_fraction_not_absorbed(std::cout, ctrw, measurer.next_measure_time());
      measurer(measurer.next_measure_time());
      std::cout << "Done!\n";
    }
    current_time = measurer.next_measure_time();
    
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
