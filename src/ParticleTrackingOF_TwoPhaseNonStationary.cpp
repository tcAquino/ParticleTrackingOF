//
//  ParticleTrackingOF_TwoPhaseNonStationary.cpp
//  PTOF
//
//  Created by Tomás Aquino on 07/02/2024.
//

#include "General/IO.h"
#include "General/Parallel.h"
#include "General/Useful.h"
#include "PTOF/Advection.h"
#include "PTOF/Directories.h"
#include "PTOF/Model.h"
#include "PTOF/Phase.h"
#include "PTOF/Transitions.h"
#include "PTOF/Useful.h"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <exception>
#include <fvcGrad.H>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>

template <typename ParallelOption> struct ExecutableInfo {
  inline static std::string banner() {
    return io::line() + "ParticleTrackingOF_TwoPhaseNonStationary\n" +
           io::line();
  };

  static int constexpr nr_parameters =
      std::is_same_v<ParallelOption, par::ParallelOptions::Serial> ? 11 : 12;

  inline static std::ostream &info(std::ostream &output) {
    output << io::line()
           << "Executable parameters (pass '' for default in []):\n"
           << io::line() << "Number of parameters: " << nr_parameters << "\n"
           << "- Cases directory [../cases]\n"
              "- Name of case\n"
              "- Name of transport parameter set\n"
              "- Name of phase parameter set\n"
              "- Name of reaction parameter set\n"
              "- Name of solver parameter set\n"
              "- Name of initial condition parameter set\n"
              "- Name of output parameter set\n"
              "- Output directory [<Case directory>/output]\n"
              "- Output file identifier [Based on parameter set names]\n"
              "- Run number (nonnegative integer to index output) [None]\n";
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Parallel>)
      output << "- Number of parallel threads\n";
    output << io::line();
    return output;
  }

  inline static std::ostream &help(std::ostream &output) {
    output << io::line()
           << "Help options (pass any number after -h or --help):\n"
           << io::line()
           << "-a / --all : All available info\n"
              "-e / --executable : Main executable info\n"
              "-g / --geometry : Geometry info\n"
              "-d / --directories-of : OpenFOAM directories info\n"
              "-t / --transport : Transport info\n"
              "-p / --phase : Phase info\n"
              "-r / --reaction : Reaction info\n"
              "-s / --solvers : Solver info\n"
              "-i / --initial-condition : Initial condition info\n"
              "-o / --output : Output info\n"
           << io::line();
    return output;
  }

  ExecutableInfo() = delete;
};

int main(int argc, char *argv[]) {
  using ParallelOption = par::ParallelOptions::Parallel;
  using Model = ptof::Model::advection_diffusion_2d;
  using Phase = ptof::Phase;
  using Definitions = Model::Definitions<ParallelOption>;

  std::cout << std::setprecision(2) << std::scientific;
  std::cout << ExecutableInfo<ParallelOption>::banner() << "\n";
  Model::info(std::cout);

  if (ptof::options_help<ExecutableInfo<ParallelOption>, Definitions::Geometry,
                         ptof::DirectoriesOF, Definitions::Transport, Phase,
                         Definitions::Reaction, Definitions::Solvers,
                         Definitions::InitialCondition, Definitions::Output>(
          std::cout, argc, argv))
    return 0;
  if (argc != ExecutableInfo<ParallelOption>::nr_parameters + 1)
    throw io::bad_parameters_help();

  std::size_t arg = 1;
  std::string dir = argv[arg++];
  std::string case_name = argv[arg++];
  std::string params_transport_name = argv[arg++];
  std::string params_phase_name = argv[arg++];
  std::string params_reaction_name = argv[arg++];
  std::string params_solvers_name = argv[arg++];
  std::string params_initial_condition_name = argv[arg++];
  std::string params_output_name = argv[arg++];
  std::string dir_output = argv[arg++];
  std::string filename_output_identifier = argv[arg++];
  std::string run_nr = argv[arg++];
  if constexpr (!std::is_same_v<ParallelOption, par::ParallelOptions::Serial>) {
    omp_set_num_threads(atoi(argv[arg++]));
    std::cout << io::line()
              << "Number of threads: " << par::get_num_threads(ParallelOption{})
              << "\n"
              << io::line();
  }
  if (!io::is_empty(run_nr))
    std::cout << "\n"
              << io::line() << "Run number: " << std::stoul(run_nr) << "\n"
              << io::line();

  ptof::Directories directories{dir, case_name, dir_output};
  std::cout << "\n";
  directories.info_runtime(std::cout);

  ptof::DirectoriesOF directories_of{directories};
  std::cout << "\n";
  directories_of.info_runtime(std::cout);

  std::cout << "\n"
            << "Setting up geometry...\n";
  auto execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::Geometry geometry{directories_of, directories, io::CoutLogger{}};
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
      Definitions::Transport::makeVelocityInterpolator_WithUninterpolated(
          geometry);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing transport parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::Transport::Parameters params_transport{
      directories, params_transport_name, geometry, velocity_field};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing phase parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Phase::Parameters params_phase{directories, params_phase_name, geometry};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up phases..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  auto carrier_phase_field = Phase::makePhaseField(geometry, params_phase);
  Phase::update_effective_velocity_field(carrier_phase_field, velocity_field,
                                         geometry, params_transport,
                                         params_phase);
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing reaction parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::Reaction::Parameters params_reaction{
      directories, params_reaction_name, geometry, params_transport};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing solver parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::Solvers::Parameters params_solvers{
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
  auto bulk_reaction = Definitions::Reaction::makeBulkReaction(
      geometry, params_reaction, params_transport, params_solvers);
  auto surface_reaction = Definitions::Reaction::makeSurfaceReaction(
      geometry, params_reaction, params_transport, params_solvers);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing initial condition parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::InitialCondition::Parameters params_initial_condition{
      directories, params_initial_condition_name, geometry, params_transport,
      params_reaction};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up initial condition...\n";
  execution_begin = std::chrono::high_resolution_clock::now();
  auto initial_condition = Definitions::InitialCondition::makeInitialCondition<
      Definitions::CTRW::Particle>(
      geometry, velocity_field, params_initial_condition, params_solvers,
      carrier_phase_field, 1. - params_phase.phase_tolerance);
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
      velocity_field, surface_reaction, initial_condition);
  boundary.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up dynamics..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::CTRW ctrw{initial_condition()};
  auto transitions = ptof::makeTransitions<Definitions::Transport>(
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
  Definitions::Output::Parameters params_output{
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
        "TwoPhaseNonStationary_" + Model::name, case_name,
        directories_of.case_name, params_transport_name, params_phase_name,
        params_reaction_name, params_solvers_name,
        params_initial_condition_name, params_output_name);
  filename_output_identifier += (io::is_empty(run_nr) ? "" : "_RUN_" + run_nr);
  auto output = Definitions::Output::makeOutput(
      ctrw, velocity_field, geometry, directories, params_output,
      filename_output_identifier, {std::cref(carrier_phase_field)},
      {params_phase.phase_tolerance});
  output.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  useful::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Starting dynamics..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  double current_time = directories_of.time.value();
  ptof::info_time(std::cout, params_output, current_time);
  auto const &flow_times = directories_of.time.times();
  auto closest_time_index = [&flow_times](auto value) {
    return std::distance(
        flow_times.cbegin(),
        std::min_element(
            flow_times.cbegin(), flow_times.cend(), [value](auto xx, auto yy) {
              if (xx.name() == "constant" || yy.name() == "constant")
                return true;
              else
                return std::abs(xx.value() - value) <
                       std::abs(yy.value() - value);
            }));
  };
  Foam::label time_index = closest_time_index(current_time);
  double next_flow_time = std::numeric_limits<double>::infinity();
  if (time_index < flow_times.size() - 1)
    next_flow_time = flow_times[time_index + 1].value();
  while (!output.done(current_time)) {
    if (time_index < flow_times.size() - 1) {
      bool velocity_field_needs_update = false;
      while (time_index < flow_times.size() - 1 &&
             current_time >= next_flow_time) {
        ++time_index;
        next_flow_time = flow_times[time_index].value();
        directories_of.time.setTime(next_flow_time, 0);
        velocity_field_needs_update = true;
      }
      if (velocity_field_needs_update) {
        std::cout << "Field updates required...\n";
        ptof::info_time(std::cout, params_output, current_time);
        std::cout << "Updating velocity and phase fields...\n";
        Phase::update_phase_and_velocity_field(carrier_phase_field,
                                               velocity_field, geometry,
                                               params_transport, params_phase);
        std::cout << "Done!\n";
        std::cout << "Updating output handlers...\n";
        output.update(current_time, next_flow_time);
        std::cout << "Done!\n";
        std::cout << "Done!\n";
      }
    }
    while (output.next_measurement_time() <= current_time) {
      std::cout << "Measurement required...\n";
      ptof::info_time(std::cout, params_output, output.next_measurement_time());
      ptof::info_fraction_not_absorbed(std::cout, ctrw,
                                       output.next_measurement_time());
      output(output.next_measurement_time());
      std::cout << "Done!\n";
    }
    double previous_time = current_time;
    current_time =
        next_flow_time > current_time
            ? std::min(output.next_measurement_time(), next_flow_time)
            : output.next_measurement_time();
    Definitions::Solvers::evolve(
        ctrw, transitions, params_solvers, current_time, previous_time,
        [&bulk_reaction, &surface_reaction](Definitions::CTRW ctrw,
                                            Definitions::State::Time new_time,
                                            Definitions::State::Time old_time) {
          Definitions::Reaction::update(bulk_reaction, surface_reaction, ctrw,
                                        new_time, old_time);
        });
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
