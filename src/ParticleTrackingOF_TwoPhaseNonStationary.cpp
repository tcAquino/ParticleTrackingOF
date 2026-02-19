//
//  ParticleTrackingOF_TwoPhaseNonStationary.cpp
//  PTOF
//
//

#include "General/Chrono.h"
#include "General/IO.h"
#include "General/Meta.h"
#include "General/Parallel.h"
#include "PTOF/Advection.h"
#include "PTOF/Directories.h"
#include "PTOF/Field.h"
#include "PTOF/Model.h"
#include "PTOF/Phase.h"
#include "PTOF/Solvers.h"
#include "PTOF/Transitions.h"
#include "PTOF/TransportHandler.h"
#include "PTOF/Useful.h"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fvcGrad.H>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>

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
                                 par::ParallelOptions::Parallel>) {
      output << "- Number of parallel threads\n";
    }
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
  using ChemicalPotentialModel = ptof::ChemicalPotentialModels::AGG;
  using TimeInterpolationType = ptof::InterpolationTypes::Linear;

  constexpr bool chemical_potential =
      !std::is_same_v<ChemicalPotentialModel,
                      ptof::ChemicalPotentialModels::None>;
  constexpr bool hard_reflection = !chemical_potential;
  using Phase = ptof::Phase<ChemicalPotentialModel, TimeInterpolationType,
                            hard_reflection>;
  using Definitions = Model::Definitions<ParallelOption>;

  constexpr bool advection = Definitions::Solvers::Parameters::advection;
  using Solvers = std::conditional_t<
      !advection && chemical_potential,
      ptof::Solvers_Generic<
          ptof::Steppers::Euler, Definitions::Solvers::Stepper_Diffusion,
          Definitions::Solvers::reaction, Definitions::Solvers::Stepper_CTRW>,
      Definitions::Solvers>;
  std::cout << std::setprecision(2) << std::scientific;
  std::cout << ExecutableInfo<ParallelOption>::banner();
  Model::info(std::cout);
  if (!advection && chemical_potential) {
    std::cout
        << io::line()
        << "Note: Turning on advection with forward Euler stepping to apply\n"
           "      chemical potential\n"
        << io::line();
  }

  if (ptof::options_help<ExecutableInfo<ParallelOption>, Definitions::Geometry,
                         ptof::DirectoriesOF, Definitions::TransportHandler,
                         Phase, Definitions::ReactionHandler, Solvers,
                         Definitions::InitialConditionHandler,
                         Definitions::OutputHandler>(std::cout, argc, argv)) {
    return 0;
  }
  if (argc != ExecutableInfo<ParallelOption>::nr_parameters + 1) {
    throw io::bad_parameters_help();
  }

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
  if (!io::is_empty(run_nr)) {
    std::cout << io::line() << "Run number: " << std::stoul(run_nr) << "\n"
              << io::line();
  }

  ptof::Directories directories{dir, case_name, dir_output};
  directories.info_runtime(std::cout);

  ptof::DirectoriesOF directories_of{directories};
  directories_of.info_runtime(std::cout);

  std::cout << "\n"
            << "Setting up geometry...\n";
  auto execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::Geometry geometry{directories_of, directories, io::CoutLogger{}};
  geometry.info_runtime(std::cout);
  auto execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing velocity data..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Foam::scalar start_time_value = directories_of.time.value();
  Foam::word start_time_name = directories_of.time.timeName();
  auto velocity_data_old = ptof::get_velocity_data(
      geometry.mesh(), meta::Selector<bool, advection>{});
  auto const &flow_times = directories_of.time.times();
  Foam::label next_time_index =
      ptof::closest_time_index(flow_times, start_time_value) + 1;
  Foam::instant next_flow_time{std::numeric_limits<double>::infinity()};
  if (next_time_index < flow_times.size()) {
    next_flow_time = flow_times[next_time_index];
  }
  auto velocity_data_new = ptof::get_velocity_data(
      geometry.mesh(),
      next_flow_time.value() < std::numeric_limits<double>::infinity()
          ? next_flow_time.name()
          : directories_of.time.name(),
      meta::Selector<bool, advection>{});
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing transport parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::TransportHandler::Parameters params_transport{
      directories, params_transport_name, geometry, velocity_data_old};
  ptof::rescale(velocity_data_new, params_transport.velocity_rescaling_factor);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up velocity interpolation..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  auto velocity_field_base =
      Definitions::TransportHandler::makeVelocityTimeInterpolator<
          Definitions::Solvers>(
          geometry, std::move(velocity_data_new), std::move(velocity_data_old),
          next_flow_time.value(), start_time_value, TimeInterpolationType{});
  std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing phase parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Phase::Parameters params_phase{directories, params_phase_name, geometry};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up phases..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  auto carrier_phase_field = Phase::makeCarrierPhaseTimeInterpolator(
      geometry,
      Phase::get_carrier_phase_data(geometry.mesh(), start_time_name,
                                    params_phase),
      Phase::get_carrier_phase_data(
          geometry.mesh(),
          next_flow_time.value() < std::numeric_limits<double>::infinity()
              ? next_flow_time.name()
              : start_time_name,
          params_phase),
      next_flow_time.value(), start_time_value);
  auto velocity_field = Phase::EffectiveVelocity{
      velocity_field_base, carrier_phase_field, params_transport, params_phase};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing reaction parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::ReactionHandler::Parameters params_reaction{
      directories, params_reaction_name, geometry, params_transport};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing solver parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Solvers::Parameters params_solvers{directories, params_solvers_name, geometry,
                                     params_transport, params_reaction};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up reaction..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  auto bulk_reaction = Definitions::ReactionHandler::makeBulkReaction(
      geometry, params_reaction, params_transport, params_solvers);
  auto surface_reaction = Definitions::ReactionHandler::makeSurfaceReaction(
      geometry, params_reaction, params_transport, params_solvers);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing initial condition parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::InitialConditionHandler::Parameters params_initial_condition{
      directories, params_initial_condition_name, geometry, params_transport,
      params_reaction};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up initial condition...\n";
  execution_begin = std::chrono::high_resolution_clock::now();
  auto initial_condition =
      Definitions::InitialConditionHandler::makeInitialCondition<
          Definitions::CTRW::Particle>(
          geometry, velocity_field, params_initial_condition, params_solvers,
          carrier_phase_field, 1. - params_phase.phase_tolerance);
  initial_condition.info_runtime(std::cout);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up boundary conditions...\n";
  execution_begin = std::chrono::high_resolution_clock::now();
  auto boundary = geometry.makeBoundary<Definitions::State>(
      directories, params_transport, params_reaction, params_solvers,
      velocity_field, surface_reaction, initial_condition);
  boundary.info_runtime(std::cout);
  auto boundary_plus_optional_hard_interface =
      Phase::BoundaryPlusOptionalHardInterface(boundary, carrier_phase_field,
                                               velocity_field, params_phase);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up dynamics..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::CTRW ctrw{initial_condition()};
  auto transitions =
      ptof::makeTransitions<Definitions::TransportHandler, Solvers>(
          velocity_field, geometry, boundary_plus_optional_hard_interface,
          params_transport, params_reaction, params_solvers, bulk_reaction);
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Importing output parameters..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  Definitions::OutputHandler::Parameters params_output{
      directories,      params_output_name, geometry,
      params_transport, params_reaction,    params_solvers};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Setting up output..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  if (io::is_empty(filename_output_identifier)) {
    filename_output_identifier = ptof::identifier(
        "TwoPhaseNonStationary_" + Model::name, case_name,
        directories_of.case_name, params_transport_name, params_phase_name,
        params_reaction_name, params_solvers_name,
        params_initial_condition_name, params_output_name);
  }
  filename_output_identifier += (io::is_empty(run_nr) ? "" : "_RUN_" + run_nr);
  auto output = Definitions::OutputHandler::makeOutput(
      ctrw, velocity_field, geometry, boundary, directories, params_output,
      filename_output_identifier, {std::cref(carrier_phase_field)},
      {params_phase.phase_tolerance});
  output.info_runtime(std::cout);
  auto masks = std::vector{std::cref(carrier_phase_field)};
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;

  std::cout << "\n"
            << "Starting dynamics..." << std::endl;
  execution_begin = std::chrono::high_resolution_clock::now();
  double previous_time = start_time_value;
  double current_time = start_time_value;
  ptof::info_time(std::cout, params_output, current_time);
  while (output.next_measurement_time() <= current_time) {
    std::cout << "Measurement required...\n";
    ptof::info_time(std::cout, params_output, output.next_measurement_time());
    ptof::info_fraction_not_absorbed(std::cout, ctrw,
                                     output.next_measurement_time());
    output(output.next_measurement_time());
    std::cout << "Done!\n";
  }
  while (current_time < std::numeric_limits<double>::infinity() &&
         !output.done(current_time)) {
    bool fields_need_update = false;
    while (next_time_index < flow_times.size() &&
           current_time >= next_flow_time.value()) {
      directories_of.time.setTime(next_flow_time, 0);
      next_flow_time = flow_times[++next_time_index];
      fields_need_update = true;
    }
    if (next_time_index == flow_times.size()) {
      next_flow_time.value() = std::numeric_limits<double>::infinity();
    }
    if (fields_need_update) {
      std::cout << "Field updates required...\n";
      ptof::info_time(std::cout, params_output, current_time);
      Phase::update_effective_velocity_and_phase_field(
          velocity_field, carrier_phase_field,
          next_flow_time.value() < std::numeric_limits<double>::infinity()
              ? next_flow_time
              : flow_times[next_time_index - 1],
          params_transport, params_phase);
      output.update(current_time, flow_times[next_time_index - 1],
                    next_flow_time.value() <
                            std::numeric_limits<double>::infinity()
                        ? next_flow_time
                        : flow_times[next_time_index - 1]);
      std::cout << "Done!\n";
    }
    previous_time = current_time;
    current_time =
        next_flow_time.value() > current_time
            ? std::min(output.next_measurement_time(), next_flow_time.value())
            : output.next_measurement_time();
    Solvers::evolve(
        ctrw, transitions, params_solvers, current_time, previous_time,
        [&bulk_reaction, &surface_reaction](Definitions::CTRW ctrw,
                                            Definitions::State::Time new_time,
                                            Definitions::State::Time old_time) {
          Definitions::ReactionHandler::update(bulk_reaction, surface_reaction,
                                               ctrw, new_time, old_time);
        });
    if (current_time == output.next_measurement_time()) {
      std::cout << "Measurement required...\n";
      ptof::info_time(std::cout, params_output, current_time);
      ptof::info_fraction_not_absorbed(std::cout, ctrw, current_time);
      output(current_time);
      std::cout << "Done!\n";
    }
  }
  if (output.done(output.next_measurement_time())) {
    std::cout << "End criterion reached." << std::endl;
  } else {
    std::cout << "Maximum measurement time reached before end criterion."
              << std::endl;
  }
  output();
  execution_end = std::chrono::high_resolution_clock::now();
  std::cout << "Done!";
  std::cout << " (";
  chrono::display_duration(std::cout, execution_begin, execution_end);
  std::cout << ")" << std::endl;
  return 0;
}
