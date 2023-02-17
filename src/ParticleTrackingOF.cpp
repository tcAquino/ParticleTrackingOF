// File: ParticleTrackingOF.cpp
// Author: Tomás Aquino
// Date: 16/02/2022

#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include "General/Useful.h"
#include "PTOF/Directories.h"
#include "PTOF/Models.h"

int main(int argc, char * argv[])
{
  using namespace ptof::model_advection_diffusion_fpt_2d;
  
  if (useful::check_options_help(argc, argv))
  {
    std::cout <<
      "--------------------------------------------------\n"
      "ParticleTrackingOF\n"
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
    Transport::Parameters::info(std::cout);
    std::cout << std::endl;
    Reaction::Parameters::info(std::cout);
    std::cout << std::endl;
    Solvers::info(std::cout);
    std::cout << std::endl;
    Solvers::Parameters::info(std::cout);
    std::cout << std::endl;
    InitialCondition::Parameters::info(std::cout);
    std::cout << std::endl;
    Output::Parameters::info(std::cout);
    std::cout << std::endl;
    return 0;
  }
  if (argc != 9)
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
  
  ptof::Directories directories{ dir, case_name, dir_output };
  directories.info_runtime(std::cout);
  std::cout << std::endl;
  ptof::DirectoriesOF directories_of{ directories };
  directories_of.info_runtime(std::cout);
  
  std::cout << std::setprecision(2) << std::scientific;
  
  std::cout << "\n" << "Importing transport parameters..." << std::endl;
  Transport::Parameters params_transport{ directories,
    parameters_transport_name };
  std::cout << "Done!" << std::endl;
  std::cout << "\n" << "Importing reaction parameters..." << std::endl;
  Reaction::Parameters params_reaction{ directories,
    parameters_reaction_name,
    params_transport };
  std::cout << "Done!" << std::endl;
  std::cout << "\n"  << "Importing solver parameters..." << std::endl;
  Solvers::Parameters params_solvers{ directories,
    parameters_solvers_name, params_transport, params_reaction };
  std::cout << "Done!" << std::endl;
  std::cout << "\n"  << "Importing initial condition parameters..." << std::endl;
  InitialCondition::Parameters params_initial_condition{ directories,
    parameters_initial_condition_name,
    params_transport, params_reaction, params_solvers };
  std::cout << "Done!" << std::endl;
  
  std::cout << "\n" << "Setting up geometry...\n";
  Geometry geometry{
    directories_of, directories, params_transport };
  geometry.info_runtime(std::cout);
  std::cout << "Done!" << std::endl;
  
  std::cout << "\n" << "Setting up velocity interpolation..." << std::endl;
  double velocity_rescaling = 1.;
  auto velocity_field = Transport::makeVelocityInterpolator(geometry);
  if (params_transport.peclet
      == std::numeric_limits<double>::infinity())
  {
    params_transport.advection_time = params_transport.lengthscale/
    ptof::magnitude_of_average(velocity_field.get_field(), geometry.mesh);
  }
  else
  {
    double average_velocity_magnitude =
      params_transport.diff_coeff*params_transport.peclet/
      params_transport.lengthscale;
    velocity_rescaling = average_velocity_magnitude/
      ptof::magnitude_of_average(velocity_field.get_field(),
                                 geometry.mesh);
    velocity_field.rescale(velocity_rescaling);
  }
  std::cout << "Done!" << std::endl;
  
  std::cout << "\n" << "Setting up reaction..." << std::endl;
  auto reaction = Reaction::makeReaction(geometry,
                                         params_reaction,
                                         params_transport,
                                         params_solvers);
  std::cout << "Done!" << std::endl;
  
  std::cout << "\n" << "Setting up initial condition...\n";
  auto initial_condition
    = InitialCondition::makeInitialCondition(geometry,
                                             velocity_field,
                                             params_initial_condition);
  initial_condition.info_runtime(std::cout);
  std::cout << "Done!" << std::endl;
  
  std::cout << "\n" << "Setting up boundary conditions...\n";
  auto boundary = geometry.makeBoundary(directories,
                                        params_transport,
                                        params_reaction,
                                        params_solvers,
                                        initial_condition);
  boundary.info_runtime(std::cout);
  std::cout << "Done!" << std::endl;
  
  std::cout << "\n" << "Setting up dynamics..." << std::endl;
  CTRW ctrw{ initial_condition(), CTRW::Tag{} };
  ctrw::Transitions_CTRW_Transport_Reaction transitions{
    Transport::makeTransitions(velocity_field,
                               geometry, boundary,
                               params_transport, params_solvers),
    reaction };
  std::cout << "Done!" << std::endl;
  
  std::cout << "\n"  << "Importing output parameters..." << std::endl;
  Output::Parameters params_output{ directories,
    parameters_output_name,
    params_transport, params_reaction, params_solvers,
    velocity_rescaling };
  std::cout << "Done!" << std::endl;
  
  std::cout << "\n" << "Setting up output..." << std::endl;
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
  std::cout << "Done!" << std::endl;

  std::cout << "\n" << "Starting dynamics..." << std::endl;
  double current_time = 0.;
  const bool step_in_time = (params_solvers.time_step != 0);
  while (!measurer.done(current_time))
  {
    while (measurer.next_measure_time() <= current_time)
    {
      std::cout << "Time "
                << "[" << params_output.time_units << " times]: "
                << measurer.next_measure_time()/params_output.time_unit_factor
                << "\n"
                << "Fraction not absorbed: "
                << 1. - double(ptof::nr_absorbed(ctrw, current_time))/ctrw.size()
                << "\n";
      measurer(measurer.next_measure_time());
    }
    if (step_in_time)
    {
      current_time += params_solvers.time_step;
      ctrw.step([current_time](CTRW::Particle& part)
                { return part.state_new().time < current_time
                    && !part.state_new().info.absorbed; },
                transitions);
    }
    else
    {
      current_time = measurer.next_measure_time();
      ctrw.evolve([current_time](CTRW::Particle& part)
      { return part.state_new().time < current_time
          && !part.state_new().info.absorbed; },
                transitions);
    }
  }
  if (measurer.next_measure_time() <= current_time)
  {
    std::cout << "Time "
              << "[" << params_output.time_units << " times]: "
              << measurer.next_measure_time()/params_output.time_unit_factor
              << "\n"
              << "Fraction not absorbed: "
              << 1. - double(ptof::nr_absorbed(ctrw, current_time))/ctrw.size()
              << "\n";
    measurer(measurer.next_measure_time());
  }
  measurer();
  
  std::cout << "Done!" << std::endl;
  
  return 0;
}
