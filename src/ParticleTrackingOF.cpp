//
//  ParticleTrackingOF.cpp
//
//  Created by Tomás Aquino on 16/02/2022.
//

#include "General/IO.h"
#include "General/Parallel.h"
#include "General/Useful.h"
#include "PTOF/Directories.h"
#include "PTOF/Models.h"
#include "PTOF/Useful.h"
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <string>

template <typename ParallelOption> struct ExecutableInfo {
  template <typename OStream> static void banner(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "ParticleTrackingOF\n"
           "--------------------------------------------------------------\n";
  }

  static int constexpr nr_parameters =
      std::is_same_v<ParallelOption, par::ParallelOptions::Serial> ? 10 : 11;

  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------------------\n"
              "Executable parameters (pass '' for default in []):\n"
              "--------------------------------------------------------------\n"
              "Number of parameters: "
           << nr_parameters << "\n"
           << "- Cases directory [../cases]\n"
              "- Name of case\n"
              "- Name of transport parameter set\n"
              "- Name of reaction parameter set\n"
              "- Name of solver parameter set\n"
              "- Name of initial condition parameter set\n"
              "- Name of output parameter set\n"
              "- Output directory [<Case directory>/output]\n"
              "- Output file identifier [Based on parameter set names]\n"
              "- Run number (nonnegative integer to index output) [None]\n";
    if constexpr (!std::is_same_v<ParallelOption, par::ParallelOptions::Serial>)
      output << "- Number of parallel threads\n";
    output
        << "--------------------------------------------------------------\n";
  }

  template <typename OStream> static void help(OStream &output) {
    output
        << "--------------------------------------------------------------\n"
           "Help options (pass any number after -h or --help):\n"
           "--------------------------------------------------------------\n"
           "-a / --all : All available info\n"
           "-e / --executable : Main executable info\n"
           "-g / --geometry : Geometry info\n"
           "-d / --directories-of : OpenFOAM directories info\n"
           "-t / --transport : Transport info\n"
           "-r / --reaction : Reaction info\n"
           "-s / --solvers : Solvers info\n"
           "-i / --initial-condition : Initial condition info\n"
           "-o / --output : Output info\n"
           "--------------------------------------------------------------\n";
  }

  ExecutableInfo() = delete;
};

int main(int argc, char *argv[]) {
  using ParallelOption = par::ParallelOptions::Parallel;
  namespace model = ptof::model_periodic_cartesian_advection_diffusion_fpt_3d;

  std::cout << std::setprecision(2) << std::scientific;
  ExecutableInfo<ParallelOption>::banner(std::cout);
  std::cout << "\n";
  model::Definitions<ParallelOption>::Model::info(std::cout);

  if (ptof::options_help<
          ExecutableInfo<ParallelOption>,
          model::Definitions<ParallelOption>::Geometry, ptof::DirectoriesOF,
          model::Definitions<ParallelOption>::Transport, useful::Empty,
          model::Definitions<ParallelOption>::Reaction,
          model::Definitions<ParallelOption>::Solvers,
          model::Definitions<ParallelOption>::InitialCondition,
          model::Definitions<ParallelOption>::Output>(std::cout, argc, argv))
    return 0;
  if (argc != ExecutableInfo<ParallelOption>::nr_parameters + 1)
    throw io::bad_parameters_help();

  std::size_t arg = 1;
  std::string dir = argv[arg++];
  std::string case_name = argv[arg++];
  std::string params_transport_name = argv[arg++];
  std::string params_reaction_name = argv[arg++];
  std::string params_solvers_name = argv[arg++];
  std::string params_initial_condition_name = argv[arg++];
  std::string params_output_name = argv[arg++];
  std::string dir_output = argv[arg++];
  std::string filename_output_identifier = argv[arg++];
  std::string run_nr = argv[arg++];
  if constexpr (!std::is_same_v<ParallelOption, par::ParallelOptions::Serial>) {
    omp_set_num_threads(atoi(argv[arg++]));
    std::cout
        << "\n"
        << "--------------------------------------------------------------\n"
           "Number of threads: "
        << par::get_num_threads(ParallelOption{})
        << "\n"
           "--------------------------------------------------------------\n";
  }
  if (!io::is_empty(run_nr))
    std::cout
        << "\n"
        << "--------------------------------------------------------------\n"
           "Run number: "
        << std::stoul(run_nr)
        << "\n"
           "--------------------------------------------------------------\n";

  ptof::Directories directories{dir, case_name, dir_output};
  std::cout << "\n";
  directories.info_runtime(std::cout);

  ptof::DirectoriesOF directories_of{directories};
  std::cout << "\n";
  directories_of.info_runtime(std::cout);

  std::cout << "\n"
            << "Setting up geometry...\n";
  auto execution_begin = std::chrono::high_resolution_clock::now();
  model::Definitions<ParallelOption>::Geometry geometry{directories_of,
                                                        directories};
  geometry.info_runtime(std::cout);
  auto execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up velocity interpolation..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  auto velocity_field =
      model::Definitions<ParallelOption>::Transport::makeVelocityInterpolator(
          geometry);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing transport parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  model::Definitions<ParallelOption>::Transport::Parameters params_transport{
    directories, params_transport_name, geometry, velocity_field};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing reaction parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  model::Definitions<ParallelOption>::Reaction::Parameters params_reaction{
      directories, params_reaction_name, geometry, params_transport};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing solver parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  model::Definitions<ParallelOption>::Solvers::Parameters params_solvers{
      directories, params_solvers_name, geometry, params_transport,
      params_reaction};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up reaction..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  auto bulk_reaction =
      model::Definitions<ParallelOption>::Reaction::makeBulkReaction(
          geometry, params_reaction, params_transport, params_solvers);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing initial condition parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  model::Definitions<ParallelOption>::InitialCondition::Parameters
      params_initial_condition{directories, params_initial_condition_name,
                               geometry, params_transport, params_reaction};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up initial condition...\n";
  execution_begin = std::chrono::high_resolution_clock::now();
  auto initial_condition = model::Definitions<ParallelOption>::
      InitialCondition::makeInitialCondition<
          model::Definitions<ParallelOption>::CTRW::Particle>(
          geometry, velocity_field, params_initial_condition, params_solvers);
  initial_condition.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up boundary conditions...\n";
  execution_begin = std::chrono::high_resolution_clock::now();
  auto boundary = geometry.makeBoundary(
      directories, params_transport, params_reaction, params_solvers,
      velocity_field,
      model::Definitions<ParallelOption>::Reaction::makeSurfaceReaction(
          geometry, params_reaction, params_transport, params_solvers),
      initial_condition);
  boundary.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up dynamics..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  model::Definitions<ParallelOption>::CTRW ctrw{initial_condition()};
  auto transitions =
      ptof::makeTransitions<model::Definitions<ParallelOption>::Transport,
                            model::Definitions<ParallelOption>::Solvers>(
          velocity_field, geometry, boundary, params_transport, params_reaction,
          params_solvers, bulk_reaction);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing output parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  model::Definitions<ParallelOption>::Output::Parameters params_output{
      directories,      params_output_name, geometry,
      params_transport, params_reaction,    params_solvers};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up output..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  if (io::is_empty(filename_output_identifier))
    filename_output_identifier = ptof::identifier(
        model::Definitions<ParallelOption>::Model::name, case_name,
        directories_of.case_name, params_transport_name, params_reaction_name,
        params_solvers_name, params_initial_condition_name, params_output_name);
  filename_output_identifier +=
      (io::is_empty(run_nr) ? "" : "_RUN_" + run_nr);
  auto output = model::Definitions<ParallelOption>::Output::makeOutput(
      ctrw, velocity_field, geometry, directories, params_output,
      filename_output_identifier);
  output.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Starting dynamics..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  double current_time = 0.;
  ptof::info_time(std::cout, params_output, current_time);
  while (!output.done(current_time)) {
    while (output.next_measurement_time() <= current_time) {
      std::cout << "Measurement required...\n";
      ptof::info_time(std::cout, params_output, output.next_measurement_time());
      ptof::info_fraction_not_absorbed(std::cout, ctrw,
                                       output.next_measurement_time());
      output(output.next_measurement_time());
      std::cout << "Done!\n";
    }
    current_time = output.next_measurement_time();
    ctrw.evolve(
        [current_time](
            model::Definitions<ParallelOption>::CTRW::Particle const &part) {
          return part.state_new().time < current_time &&
                 !part.state_new().info.absorbed;
        },
        transitions);
  }
  if (output.next_measurement_time() <= current_time) {
    std::cout << "Measurement required...\n";
    ptof::info_time(std::cout, params_output, output.next_measurement_time());
    ptof::info_fraction_not_absorbed(std::cout, ctrw,
                                     output.next_measurement_time());
    output(output.next_measurement_time());
    std::cout << "Done!\n";
  }
  output();
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  return 0;
}
