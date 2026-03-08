/**
 * @file   TransportParameters.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Sun Jan 19 00:00:00 2025
 *
 * @brief  Transport parameters and utilities.
 */

#ifndef PTOF_TRANSPORTPARAMETERS_H
#define PTOF_TRANSPORTPARAMETERS_H

#include "General/IO.h"
#include "General/Meta.h"
#include "PTOF/Directories.h"
#include <limits>
#include <string>
#include <utility>

namespace ptof {
/** @brief Transport parameters for advective-diffusive transport. */
struct TransportParameters_AdvectionDiffusion {
public:
  std::string peclet_option;
  double lengthscale;
  double diff_coeff;
  double diffusion_time;
  double advection_time;
  double peclet;
  double mean_velocity;
  double velocity_rescaling_factor;

  template <typename Geometry, typename VelocityData>
  TransportParameters_AdvectionDiffusion(Directories const &directories,
                                         std::string const &parameter_set_name,
                                         Geometry const &geometry,
                                         VelocityData &velocity_data) {
    std::string filename = directories.dir_parameters + "/transport_" +
                           parameter_set_name + ".param";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse reference lengthscale", lengthscale);

    split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse Peclet setting option", peclet_option);
    std::string for_peclet_option =
        std::string{"Peclet option "} + peclet_option + " : ";
    if (peclet_option == "rescale_velocity_to_peclet") {
      io::read(split_line, param_index,
               in_file + for_peclet_option +
                   "Could not parse Peclet number and diffusion coefficient",
               peclet, diff_coeff);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      advection_time = 2. * diffusion_time / peclet;
      mean_velocity = lengthscale / advection_time;
    } else if (peclet_option == "rescale_velocity_to_mean") {
      io::read(split_line, param_index,
               in_file + for_peclet_option +
                   "Could not parse mean velocity and diffusion coefficient",
               mean_velocity, diff_coeff);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      peclet = lengthscale * mean_velocity / diff_coeff;
      advection_time = lengthscale / mean_velocity;
    } else if (peclet_option == "rescale_velocity_to_advection_time") {
      io::read(split_line, param_index,
               in_file + for_peclet_option +
                   "Could not parse advection time and diffusion coefficient",
               advection_time, diff_coeff);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      peclet = 2. * diffusion_time / advection_time;
      mean_velocity = lengthscale / advection_time;
    } else if (peclet_option == "compute_from_diffusion_coefficient") {
      io::read(split_line, param_index,
               in_file + for_peclet_option +
                   "Could not parse diffusion coefficient",
               diff_coeff);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
    } else if (peclet_option == "compute_from_diffusion_time") {
      io::read(split_line, param_index,
               in_file + for_peclet_option + "Could not parse diffusion time",
               diffusion_time);
      diff_coeff = lengthscale * lengthscale / (2. * diffusion_time);
    } else if (peclet_option == "compute_diffusion_coefficient") {
      io::read(split_line, param_index,
               in_file + for_peclet_option + "Could not parse Peclet number",
               peclet);
    } else {
      throw std::runtime_error{in_file + for_peclet_option + "Not supported"};
    }

    rescale(velocity_data, geometry.mesh(), in_file + for_peclet_option);
  }

  template <typename VelocityData, typename Mesh>
  void rescale(VelocityData &velocity_data, Mesh const &mesh,
               std::string const &in_file_for_peclet_option) {
    double current_mean = magnitude_of_average(velocity_data, mesh);
    if (peclet_option == "rescale_velocity_to_peclet" ||
        peclet_option == "rescale_velocity_to_mean" ||
        peclet_option == "rescale_velocity_to_advection_time") {
      velocity_rescaling_factor = mean_velocity / current_mean;
      velocity_data *= velocity_rescaling_factor;
      velocity_data.boundaryFieldRef() ==
          velocity_rescaling_factor *velocity_data.boundaryField();
    } else if (peclet_option == "compute_from_diffusion_coefficient" ||
               peclet_option == "compute_from_diffusion_time") {
      velocity_rescaling_factor = 1.;
      mean_velocity = current_mean;
      advection_time = lengthscale / mean_velocity;
      peclet = 2. * diffusion_time / advection_time;
    } else if (peclet_option == "compute_diffusion_coefficient") {
      velocity_rescaling_factor = 1.;
      mean_velocity = current_mean;
      diff_coeff = lengthscale * mean_velocity / peclet;
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      advection_time = lengthscale / mean_velocity;
    } else {
      throw std::runtime_error{in_file_for_peclet_option + " Not supported"};
    }
  }

  /**
   * @brief Output generic information about object.
   *
   * @param output Output stream.
   */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Transport parameters\n"
           << io::line() << R"((Note:
  - Diffusion time is defined as <lengthscale>^2 / (2 * <diffusion
    coefficient>).
  - Advection time is defined as <lengthscale> / <mean velocity>.
  - Mean velocity is defined as the absolute value of mean velocity vector.
  - Peclet number is defined as <lengthscale> * <mean velocity> / <diffusion
    coefficient>.
- Reference lengthscale.
- How to set the Peclet number:
  (Note:
    - Options that involve rescaling the velocity field are expensive for
      simulations with time-dependent flow fields because they require
      recomputing the interpolator.)
  - compute_from_diffusion_coefficient
    - Compute from given diffusion coefficient (do not rescale velocity field).
    - Pass on same line:
      - Diffusion coefficient.
  - compute_diffusion_coefficient
    - Compute diffusion coefficient to impose given Peclet number (do not
      rescale velocity field).
    - Pass on same line:
      - Peclet number.
  - compute_from_diffusion_time
    - Compute from given diffusion time (do not rescale velocity field).
    - Pass on same line:
      - Diffusion coefficient.
  - rescale_velocity_to_peclet
    - Rescale velocity field according to given peclet number.
    - Pass on same line:
      - Peclet number.
      - Diffusion coefficient.
  - rescale_velocity_to_mean
    - Rescale according to given mean flow velocity.
    - Pass on same line:
      - Absolute value of mean velocity vector.
      - Diffusion coefficient.
  - rescale_velocity_to_advection_time
    - Rescale according to given advection time.
    - Pass on same line:
      - Advection time, based on absolute value of mean velocity vector.
      - Diffusion coefficient.
)" << io::line();
    return output;
  }
};

/** @brief Transport parameters for advective transport. */
struct TransportParameters_Advection {
public:
  std::string rescale_velocity_option;
  double lengthscale;
  const double diff_coeff{0.};
  const double diffusion_time{std::numeric_limits<double>::infinity()};
  double advection_time;
  const double peclet{std::numeric_limits<double>::infinity()};
  double mean_velocity;
  double velocity_rescaling_factor;

  template <typename Geometry, typename VelocityData>
  TransportParameters_Advection(Directories const &directories,
                                std::string const &parameter_set_name,
                                Geometry const &geometry,
                                VelocityData &velocity_data) {
    std::string filename = directories.dir_parameters + "/transport_" +
                           parameter_set_name + ".param";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse reference length scale", lengthscale);

    split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not flow velocity field rescaling option",
             rescale_velocity_option);
    std::string for_rescale_velocity_option =
        std::string{"Velocity rescaling option "} + rescale_velocity_option +
        " : ";
    if (rescale_velocity_option == "rescale_velocity_to_mean") {
      io::read(split_line, param_index,
               in_file + for_rescale_velocity_option +
                   "Could not parse mean velocity",
               mean_velocity);
      advection_time = lengthscale / mean_velocity;
    } else if (rescale_velocity_option ==
               "rescale_velocity_to_advection_time") {
      io::read(split_line, param_index,
               in_file + for_rescale_velocity_option +
                   "Could not parse advection time",
               advection_time);
      mean_velocity = lengthscale / advection_time;
    } else if (rescale_velocity_option == "no_rescale_velocity") {
    } else {
      throw std::runtime_error{in_file + for_rescale_velocity_option +
                               "Not supported"};
    }

    rescale(velocity_data, geometry.mesh(),
            in_file + for_rescale_velocity_option);
  }

  template <typename VelocityData, typename Mesh>
  void rescale(VelocityData &velocity_data, Mesh const &mesh,
               std::string const &in_file_for_rescale_velocity_option) {
    double current_mean = magnitude_of_average(velocity_data, mesh);
    if (rescale_velocity_option == "rescale_velocity_to_mean" ||
        rescale_velocity_option == "rescale_velocity_to_advection_time") {
      velocity_rescaling_factor = mean_velocity / current_mean;
      velocity_data *= velocity_rescaling_factor;
      velocity_data.boundaryFieldRef() ==
          velocity_rescaling_factor *velocity_data.boundaryField();
    } else if (rescale_velocity_option == "no_rescale_velocity") {
      velocity_rescaling_factor = 1.;
      advection_time = lengthscale / current_mean;
      mean_velocity = lengthscale / advection_time;
    } else {
      throw std::runtime_error{in_file_for_rescale_velocity_option +
                               "Not supported"};
    }
  }

  /**
   * @brief Output generic information about object.
   *
   * @param output Output stream.
   */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Transport parameters\n"
           << io::line() << R"((Note:
  - Advection time is defined as <lengthscale> / <mean velocity>.
  - Mean velocity is defined as the absolute value of mean velocity vector.)
- Reference lengthscale.
- Whether and how to rescale velocity field:
  (Note:
    - Options that involve rescaling the velocity field are expensive for
      simulations with time-dependent flow fields because they require
      recomputing the interpolator.)
  - rescale_velocity_to_mean
    - Rescale according to given mean.
    - Pass on same line:
      - Absolute value of mean velocity vector.
  - rescale_velocity_to_advection_time
    - Rescale according to given advection time.
    - Pass on same line:
      - Advection time, based on absolute value of mean velocity vector.
  - no_rescale_velocity
    - Do not rescale.
)" << io::line();
    return output;
  }
};

/** @brief Transport parameters for diffusive transport. */
struct TransportParameters_Diffusion {
public:
  double lengthscale;
  double diff_coeff;
  double diffusion_time;
  const double advection_time{std::numeric_limits<double>::infinity()};
  const double peclet{0.};
  const double mean_velocity{0.};
  const double velocity_rescaling_factor{1.};

  template <typename Geometry, typename VelocityData>
  TransportParameters_Diffusion(Directories const &directories,
                                std::string const &parameter_set_name,
                                Geometry const &geometry,
                                VelocityData &velocity_data) {
    std::string filename = directories.dir_parameters + "/transport_" +
                           parameter_set_name + ".param";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse reference length scale", lengthscale);

    split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    param_index = 0;
    auto diffusion_option = io::read<std::string>(
        split_line, param_index,
        in_file + "Could not parse diffusion coefficient setting option");
    std::string for_diffusion_option =
        std::string{"Diff option "} + diffusion_option + " : ";
    if (diffusion_option == "set_value") {
      io::read(split_line, param_index,
               in_file + for_diffusion_option +
                   "Could not parse diffusion coefficient",
               diff_coeff);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
    } else if (diffusion_option == "compute_from_diffusion_time") {
      io::read(split_line, param_index,
               in_file + for_diffusion_option +
                   "Could not parse diffusion time",
               diffusion_time);
      diff_coeff = lengthscale * lengthscale / (2. * diffusion_time);
    } else {
      throw std::runtime_error{in_file + for_diffusion_option +
                               "Not supported"};
    }
  }

  /**
   * @brief Output generic information about object.
   *
   * @param output Output stream.
   */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Transport parameters\n"
           << io::line() << R"((Note:
  - Diffusion time is defined as <lengthscale>^2 / (2 * <diffusion
    coefficient>).)
- Reference lengthscale.
- How to set the diffusion coefficient:
  - set_value
    - Set given value.
    - Pass on same line:
      - Diffusion coefficient.
  - compute_from_diffusion_time
    - Compute from given diffusion time.
    - Pass on same line:
      - Diffusion time.
)" << io::line();
    return output;
  }
};

/**
 * @brief Transport parameters for advective-diffusive transport in the
 *        primitive cell of a body centered cubic packing.
 */
struct TransportParameters_AdvectionDiffusion_Bcc {
public:
  std::string lengthscale_option;
  std::string peclet_option;
  double lengthscale;
  double peclet;
  double diff_coeff;
  double diffusion_time;
  double mean_velocity;
  double velocity_rescaling_factor;
  double cell_side;
  std::vector<std::pair<double, double>> primitive_cell_boundaries;
  double advection_time;

  template <typename Geometry, typename VelocityData>
  TransportParameters_AdvectionDiffusion_Bcc(
      Directories const &directories, std::string const &parameter_set_name,
      Geometry const &geometry, VelocityData &velocity_data) {
    std::string filename = directories.dir_parameters + "/transport_" +
                           parameter_set_name + ".param";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file +
                 "Could not parse reference lengthscale definition option",
             lengthscale_option);
    std::string for_lengthscale_option =
        std::string{"Lengthscale definition option "} + lengthscale_option +
        " : ";

    double radius = geometry.radius;
    cell_side = 4. / std::sqrt(3.) * radius;
    if (lengthscale_option == "radius") {
      lengthscale = radius;
    } else if (lengthscale_option == "diameter")
      lengthscale = 2. * radius;
    else if (lengthscale_option == "cell_side") {
      lengthscale = cell_side;
    } else if (lengthscale_option == "custom") {
      io::read(split_line, param_index,
               in_file + for_lengthscale_option +
                   "Could not parse reference lengthscale",
               lengthscale);
      lengthscale = std::stod(split_line[param_index++]);
    } else {
      throw std::runtime_error{in_file + for_lengthscale_option +
                               "Not supported"};
    }

    split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse Peclet setting option", peclet_option);
    std::string for_peclet_option =
        std::string{"Peclet option "} + peclet_option + " : ";
    if (peclet_option == "rescale_velocity_to_peclet") {
      io::read(split_line, param_index,
               in_file + for_peclet_option +
                   "Could not parse Peclet number and diffusion coefficient",
               peclet, diff_coeff);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      advection_time = 2. * diffusion_time / peclet;
      mean_velocity = lengthscale / advection_time;
    } else if (peclet_option == "rescale_velocity_to_mean") {
      io::read(split_line, param_index,
               in_file + for_peclet_option +
                   "Could not parse mean velocity and diffusion coefficient",
               mean_velocity, mean_velocity);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      peclet = lengthscale * mean_velocity / diff_coeff;
      advection_time = lengthscale / mean_velocity;
    } else if (peclet_option == "rescale_velocity_to_advection_time") {
      io::read(split_line, param_index,
               in_file + for_peclet_option +
                   "Could not parse advection time and diffusion coefficient",
               advection_time, diff_coeff);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      peclet = 2. * diffusion_time / advection_time;
      mean_velocity = lengthscale / advection_time;
    } else if (peclet_option == "compute_from_diffusion_coefficient") {
      io::read(split_line, param_index,
               in_file + for_peclet_option +
                   "Could not parse diffusion coefficient",
               diff_coeff);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
    } else if (peclet_option == "compute_from_diffusion_time") {
      io::read(split_line, param_index,
               in_file + for_peclet_option + "Could not parse diffusion time",
               diffusion_time);
      diff_coeff = lengthscale * lengthscale / (2. * diffusion_time);
    } else if (peclet_option == "compute_diffusion_coefficient") {
      io::read(split_line, param_index,
               in_file + for_peclet_option + "Could not parse Peclet number",
               peclet);
    } else {
      throw std::runtime_error{in_file + for_peclet_option + "Not supported"};
    }

    rescale(velocity_data, geometry.mesh(), in_file + for_peclet_option);
  }

  template <typename VelocityData, typename Mesh>
  void rescale(VelocityData &velocity_data, Mesh const &mesh,
               std::string const &in_file_for_peclet_option) {
    double current_mean = magnitude_of_average(velocity_data, mesh);
    if (peclet_option == "rescale_velocity_to_peclet" ||
        peclet_option == "rescale_velocity_to_mean" ||
        peclet_option == "rescale_velocity_to_advection_time") {
      velocity_rescaling_factor = mean_velocity / current_mean;
      velocity_data *= velocity_rescaling_factor;
      velocity_data.boundaryFieldRef() ==
          velocity_rescaling_factor *velocity_data.boundaryField();
    } else if (peclet_option == "compute_from_diffusion_coefficient" ||
               peclet_option == "compute_from_diffusion_time") {
      velocity_rescaling_factor = 1.;
      mean_velocity = current_mean;
      advection_time = lengthscale / mean_velocity;
      peclet = 2. * diffusion_time / advection_time;
    } else if (peclet_option == "compute_diffusion_coefficient") {
      velocity_rescaling_factor = 1.;
      mean_velocity = current_mean;
      diff_coeff = lengthscale * mean_velocity / peclet;
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      advection_time = lengthscale / mean_velocity;
    } else {
      throw std::runtime_error{in_file_for_peclet_option + "Not supported"};
    }
  }

  /**
   * @brief Output generic information about object.
   *
   * @param output Output stream.
   */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Transport parameters\n"
           << io::line() << R"((Note:
  - Diffusion time is defined as <lengthscale>^2 / (2 * <diffusion
    coefficient>).
  - Advection time is defined as <lengthscale> / <mean velocity>.
  - Mean velocity is defined as the absolute value of mean velocity vector.
  - Peclet number is defined as <lengthscale> * <mean velocity> / <diffusion
    coefficient>.)
- How to set the reference lengthscale:
  - radius
    - Equal to bead radius.
  - diameter
    - Equal to bead diameter.
  - cell_side
    - Equal to primitive cell side.
  - custom
    - Equal to specified value.
    - Pass on same line:
      - Reference lengthscale value.
- How to set the Peclet number:
  (Note:
    - Options that involve rescaling the velocity field are expensive for
      simulations with time-dependent flow fields because they require
      recomputing the interpolator.)
  - compute_from_diffusion_coefficient
    - Compute from given diffusion coefficient (do not rescale velocity field).
    - Pass on same line:
      - Diffusion coefficient.
  - compute_diffusion_coefficient
    - Compute diffusion coefficient to impose given Peclet number (do not
      rescale velocity field).
    - Pass on same line:
      - Peclet number.
  - compute_from_diffusion_time
    - Compute from given diffusion time (do not rescale velocity field).
    - Pass on same line:
      - Diffusion coefficient.
  - rescale_velocity_to_peclet
    - Rescale velocity field according to given peclet number.
    - Pass on same line:
      - Peclet number.
      - Diffusion coefficient.
  - rescale_velocity_to_mean
    - Rescale according to given mean flow velocity.
    - Pass on same line:
      - Absolute value of mean velocity vector.
      - Diffusion coefficient.
  - rescale_velocity_to_advection_time
    - Rescale according to given advection time.
    - Pass on same line:
      - Advection time, based on absolute value of mean velocity vector.
      - Diffusion coefficient.
)" << io::line();
    return output;
  }
};

/**
 * @brief Transport parameters for advective transport in the primitive cell of
 *        a body centered cubic packing.
 */
struct TransportParameters_Advection_Bcc {
public:
  std::string lengthscale_option;
  std::string rescale_velocity_option;
  double lengthscale;
  const double diff_coeff{0.};
  const double diffusion_time{std::numeric_limits<double>::infinity()};
  const double peclet{std::numeric_limits<double>::infinity()};
  double mean_velocity;
  double velocity_rescaling_factor;
  double cell_side;
  std::vector<std::pair<double, double>> primitive_cell_boundaries;
  double advection_time;

  template <typename Geometry, typename VelocityData>
  TransportParameters_Advection_Bcc(Directories const &directories,
                                    std::string const &parameter_set_name,
                                    Geometry const &geometry,
                                    VelocityData &velocity_data) {
    std::string filename = directories.dir_parameters + "/transport_" +
                           parameter_set_name + ".param";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse lengthscale definition option",
             lengthscale_option);
    std::string for_lengthscale_option =
        std::string{"Lengthscale definition option "} + lengthscale_option +
        " : ";

    double radius = geometry.radius;
    cell_side = 4. / std::sqrt(3.) * radius;
    if (lengthscale_option == "radius") {
      lengthscale = radius;
    } else if (lengthscale_option == "diameter") {
      lengthscale = 2. * radius;
    } else if (lengthscale_option == "cell_side") {
      lengthscale = cell_side;
    } else if (lengthscale_option == "custom") {
      io::read(split_line, param_index,
               in_file + for_lengthscale_option +
                   "Could not parse reference lengthscale",
               lengthscale);
      lengthscale = std::stod(split_line[param_index++]);
    } else {
      throw std::runtime_error{in_file + for_lengthscale_option +
                               "Not supported"};
    }

    split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not flow velocity field rescaling option",
             rescale_velocity_option);
    std::string for_rescale_velocity_option =
        std::string{"Velocity rescaling option "} + rescale_velocity_option +
        " : ";
    if (rescale_velocity_option == "rescale_velocity_to_mean") {
      io::read(split_line, param_index,
               in_file + for_rescale_velocity_option +
                   "Could not parse mean velocity",
               mean_velocity);
      advection_time = lengthscale / mean_velocity;
    } else if (rescale_velocity_option ==
               "rescale_velocity_to_advection_time") {
      io::read(split_line, param_index,
               in_file + for_rescale_velocity_option +
                   "Could not parse advection time",
               advection_time);
      mean_velocity = lengthscale / advection_time;
    } else if (rescale_velocity_option == "no_rescale_velocity") {
    } else {
      throw std::runtime_error{in_file + for_rescale_velocity_option +
                               "Not supported"};
    }

    rescale(velocity_data, geometry.mesh(),
            in_file + for_rescale_velocity_option);
  }

  template <typename VelocityData, typename Mesh>
  void rescale(VelocityData &velocity_data, Mesh const &mesh,
               std::string const &in_file_for_rescale_velocity_option) {
    double current_mean = magnitude_of_average(velocity_data, mesh);
    if (rescale_velocity_option == "rescale_velocity_to_mean" ||
        rescale_velocity_option == "rescale_velocity_to_advection_time") {
      velocity_rescaling_factor = mean_velocity / current_mean;
      velocity_data *= velocity_rescaling_factor;
      velocity_data.boundaryFieldRef() ==
          velocity_rescaling_factor *velocity_data.boundaryField();
    } else if (rescale_velocity_option == "no_rescale_velocity") {
      velocity_rescaling_factor = 1.;
      advection_time = lengthscale / current_mean;
      mean_velocity = lengthscale / advection_time;
    } else {
      throw std::runtime_error{in_file_for_rescale_velocity_option +
                               "Not supported"};
    }
  }

  /**
   * @brief Output generic information about object.
   *
   * @param output Output stream.
   */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Transport parameters\n"
           << io::line() << R"((Note:
  - Advection time is defined as <lenghtscale> / <mean velocity>.
  - Mean velocity is defined as the absolute value of mean velocity vector.)
- How to set the reference lengthscale:
  - radius
    - Equal to bead radius.
  - diameter
    - Equal to bead diameter.
  - cell_side
    - Equal to primitive cell side.
  - custom
    - Equal to specified value.
    - Pass on same line:
      - Reference lengthscale value.
- Whether and how to rescale velocity field:
  (Note:
    - Options that involve rescaling the velocity field are expensive for
      simulations with time-dependent flow fields because they require
      recomputing the interpolator.)
  - rescale_velocity_to_mean
    - Rescale according to given mean.
    - Pass on same line:
      - Absolute value of mean velocity vector.
  - rescale_velocity_to_advection_time
    - Rescale according to given advection time.
    - Pass on same line:
      - Advection time, based on absolute value of mean velocity vector.
  - no_rescale_velocity
    - Do not rescale.
)" << io::line();
    return output;
  }
};

/**
 * @brief Transport parameters for diffusive transport in the primitive cell of
 *        a body centered cubic packing.
 */
struct TransportParameters_Diffusion_Bcc {
public:
  std::string lengthscale_option;
  double lengthscale;
  double diff_coeff;
  double diffusion_time;
  double cell_side;
  std::vector<std::pair<double, double>> primitive_cell_boundaries;
  const double peclet{0.};
  const double mean_velocity{0.};
  const double velocity_rescaling_factor{1.};
  const double advection_time{std::numeric_limits<double>::infinity()};

  template <typename Geometry, typename VelocityData>
  TransportParameters_Diffusion_Bcc(
      Directories const &directories, std::string const &parameter_set_name,
      Geometry const &geometry, VelocityData &&velocity_data = meta::Empty{}) {
    std::string filename = directories.dir_parameters + "/transport_" +
                           parameter_set_name + ".param";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse lengthscale definition option",
             lengthscale_option);
    std::string for_lengthscale_option =
        std::string{"Lengthscale definition option "} + lengthscale_option +
        " : ";

    double radius = geometry.radius;
    cell_side = 4. / std::sqrt(3.) * radius;
    if (lengthscale_option == "radius") {
      lengthscale = radius;
    } else if (lengthscale_option == "diameter")
      lengthscale = 2. * radius;
    else if (lengthscale_option == "cell_side") {
      lengthscale = cell_side;
    } else if (lengthscale_option == "custom") {
      io::read(split_line, param_index,
               in_file + for_lengthscale_option +
                   "Could not parse reference lengthscale",
               lengthscale);
      lengthscale = std::stod(split_line[param_index++]);
    } else {
      throw std::runtime_error{in_file + for_lengthscale_option +
                               "Not supported"};
    }

    split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
    param_index = 0;
    auto diffusion_option = io::read<std::string>(
        split_line, param_index,
        in_file + "Could not parse diffusion coefficient setting option");
    std::string for_diffusion_option =
        std::string{"Diffusion option "} + diffusion_option + " : ";
    if (diffusion_option == "set_value") {
      io::read(split_line, param_index,
               in_file + for_diffusion_option +
                   "Could not parse diffusion coefficient",
               diff_coeff);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
    } else if (diffusion_option == "compute_from_diffusion_time") {
      io::read(split_line, param_index,
               in_file + for_diffusion_option +
                   "Could not parse diffusion time",
               diffusion_time);
      diff_coeff = lengthscale * lengthscale / (2. * diffusion_time);
    } else {
      throw std::runtime_error{in_file + for_diffusion_option +
                               "Not supported"};
    }
  }

  /**
   * @brief Output generic information about object.
   *
   * @param output Output stream.
   */
  inline static std::ostream &info(std::ostream &output) {
    output << io::line() << "Transport parameters\n"
           << io::line() << R"((Note:
  - Diffusion time is defined as <lengthscale>^2 / (2 * <diffusion
    coefficient>).)
- How to set the reference lengthscale:
  - radius
    - Equal to bead radius.
  - diameter
    - Equal to bead diameter.
  - cell_side
    - Equal to primitive cell side.
  - custom
    - Equal to specified value.
    - Pass on same line:
      - Reference lengthscale value.
- How to set the diffusion coefficient:
  - set_value
    - Set given value.
    - Pass on same line:
      - Diffusion coefficient.
  - compute_from_diffusion_time
    - Compute from given diffusion time.
    - Pass on same line:
      - Diffusion time.
)" << io::line();
    return output;
  }
};

template <bool advection, bool diffusion>
using TransportParameters_Generic = std::conditional_t<
    advection && diffusion, TransportParameters_AdvectionDiffusion,
    std::conditional_t<advection, TransportParameters_Advection,
                       TransportParameters_Diffusion>>;

template <bool advection, bool diffusion>
using TransportParameters_Bcc = std::conditional_t<
    advection && diffusion, TransportParameters_AdvectionDiffusion_Bcc,
    std::conditional_t<advection, TransportParameters_Advection_Bcc,
                       TransportParameters_Diffusion_Bcc>>;
} // namespace ptof

#endif /* PTOF_TRANSPORTPARAMETERS_H */
