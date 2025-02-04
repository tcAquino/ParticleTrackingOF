/**
   \file PTOF/Transport.h
   \author Tomás Aquino
   \date 19/01/2025
   \brief Transport parameters and utilities.
*/

#ifndef PTOF_TRANSPORT_H
#define PTOF_TRANSPORT_H

#include "General/IO.h"
#include "General/Meta.h"
#include "PTOF/Directories.h"
#include "PTOF/Transitions.h"
#include <limits>
#include <string>
#include <type_traits>
#include <utility>

namespace ptof {
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

  template <typename Geometry, typename VelocityField>
  TransportParameters_AdvectionDiffusion(Directories const &directories,
                                         std::string const &parameter_set_name,
                                         Geometry const &geometry,
                                         VelocityField &velocity_field) {
    std::string filename = directories.dir_parameters +
                           "/parameters_transport_" + parameter_set_name +
                           ".dat";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse reference lengthscale", lengthscale);

    split_line = io::split_line(input);
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
    } else if (peclet_option == "compute_from_diff_coeff") {
      io::read(split_line, param_index,
               in_file + for_peclet_option +
                   "Could not parse diffusion coefficient",
               diff_coeff);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
    } else if (peclet_option == "compute_from_diff_time") {
      io::read(split_line, param_index,
               in_file + for_peclet_option + "Could not parse diffusion time",
               diffusion_time);
      diff_coeff = lengthscale * lengthscale / (2. * diffusion_time);
    } else if (peclet_option == "set_diff_coeff") {
      io::read(split_line, param_index,
               in_file + for_peclet_option + "Could not parse Peclet number",
               peclet);
    } else {
      throw std::runtime_error{in_file + for_peclet_option + "Not supported"};
    }

    rescale(velocity_field, geometry.mesh());
  }

  template <typename VelocityField, typename Mesh>
  void rescale(VelocityField &velocity_field, Mesh const &mesh) {
    double current_mean = magnitude_of_average(velocity_field.field(), mesh);
    std::cout << "Mean velocity = " << current_mean << std::endl;
    if (peclet_option == "rescale_velocity_to_peclet" ||
        peclet_option == "rescale_velocity_to_mean" ||
        peclet_option == "rescale_velocity_to_advection_time") {
      velocity_rescaling_factor = mean_velocity / current_mean;
      velocity_field.rescale(velocity_rescaling_factor);
    } else if (peclet_option == "compute_from_diff_coeff" ||
               peclet_option == "compute_from_diff_time") {
      velocity_rescaling_factor = 1.;
      mean_velocity = current_mean;
      advection_time = lengthscale / mean_velocity;
      peclet = 2. * diffusion_time / advection_time;
    } else if (peclet_option == "set_diff_coeff") {
      velocity_rescaling_factor = 1.;
      mean_velocity = current_mean;
      diff_coeff = lengthscale * mean_velocity / peclet;
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      advection_time = lengthscale / mean_velocity;
    } else
      throw std::runtime_error{"Peclet number setting option " + peclet_option +
                               " not supported"};
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << "--------------------------------------------------------------"
              "\n"
              "Transport parameters\n"
              "--------------------------------------------------------------"
              "\n"
              "- Reference lengthscale\n"
              "- How to set the Peclet number:\n"
              "  - compute_from_diff_coeff\n"
              "    - Compute from given diffusion coefficient (do not\n"
              "      rescale velocity field)\n"
              "    - Pass on same line:\n"
              "      - Diffusion coefficient\n"
              "  - set_diff_coeff\n"
              "    - Set diffusion coefficient to impose given Peclet number\n"
              "      (do not rescale velocity field)\n"
              "    - Pass on same line:\n"
              "      - Peclet number\n"
              "  - compute_from_diff_time\n"
              "    - Compute from given diffusion time (do not rescale\n"
              "      velocity field)\n"
              "    - Pass on same line:\n"
              "      - Diffusion coefficient\n"
              "  - rescale_velocity_to_peclet\n"
              "    - Rescale velocity field according to given peclet number\n"
              "    - Pass on same line:\n"
              "      - Peclet number\n"
              "      - Diffusion coefficient\n"
              "  - rescale_velocity_to_mean\n"
              "    - Rescale according to given mean flow velocity\n"
              "    - Pass on same line:\n"
              "      - Absolute value of mean velocity vector\n"
              "      - Diffusion coefficient\n"
              "  - rescale_velocity_to_advection_time\n"
              "    - Rescale according to given advection time\n"
              "    - Pass on same line:\n"
              "      - Advection time, based on absolute value of mean\n"
              "        velocity vector\n"
              "      - Diffusion coefficient\n"
              "--------------------------------------------------------------"
              "\n";
    return output;
  }
};

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

  template <typename Geometry, typename VelocityField>
  TransportParameters_Advection(Directories const &directories,
                                std::string const &parameter_set_name,
                                Geometry const &geometry,
                                VelocityField &velocity_field) {
    std::string filename = directories.dir_parameters +
                           "/parameters_transport_" + parameter_set_name +
                           ".dat";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse reference length scale", lengthscale);

    split_line = io::split_line(input);
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

    rescale(velocity_field, geometry.mesh());
  }

  template <typename VelocityField, typename Mesh>
  void rescale(VelocityField &velocity_field, Mesh const &mesh) {
    double current_mean = magnitude_of_average(velocity_field.field(), mesh);
    if (rescale_velocity_option == "rescale_velocity_to_mean" ||
        rescale_velocity_option == "rescale_velocity_to_advection_time") {
      velocity_rescaling_factor = mean_velocity / current_mean;
      velocity_field.rescale(velocity_rescaling_factor);
    } else if (rescale_velocity_option == "no_rescale_velocity") {
      velocity_rescaling_factor = 1.;
      advection_time = lengthscale / current_mean;
      mean_velocity = lengthscale / advection_time;
    } else
      throw std::runtime_error{"Flow velocity field rescaling option " +
                               rescale_velocity_option + " not supported"};
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << "--------------------------------------------------------------"
              "\n"
              "Transport parameters\n"
              "--------------------------------------------------------------"
              "\n"
              "- Reference lengthscale\n"
              "- Whether and how to rescale velocity field:\n"
              "  - rescale_velocity_to_mean\n"
              "    - Rescale according to given mean\n"
              "    - Pass on same line:\n"
              "      - Absolute value of mean velocity vector\n"
              "  - rescale_velocity_to_advection_time\n"
              "    - Rescale according to given advection time\n"
              "    - Pass on same line:\n"
              "      - Advection time, based on absolute value of mean "
              "        velocity vector\n"
              "  - no_rescale_velocity\n"
              "    - Do not rescale\n"
              "--------------------------------------------------------------"
              "\n";
    return output;
  }
};

struct TransportParameters_Diffusion {
public:
  double lengthscale;
  const double diff_coeff;
  const double diffusion_time{std::numeric_limits<double>::infinity()};
  double advection_time{std::numeric_limits<double>::infinity()};
  const double peclet{0.};
  double mean_velocity{0.};

  template <typename Geometry, typename VelocityField>
  TransportParameters_Diffusion(
      Directories const &directories, std::string const &parameter_set_name,
      Geometry const &geometry,
      VelocityField &&velocity_field = meta::Empty{}) {
    std::string filename = directories.dir_parameters +
                           "/parameters_transport_" + parameter_set_name +
                           ".dat";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input);
    std::size_t param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse reference length scale", lengthscale);

    split_line = io::split_line(input);
    param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse diffusion coefficient", diff_coeff);
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << "--------------------------------------------------------------"
              "\n"
              "Transport parameters\n"
              "--------------------------------------------------------------"
              "\n"
              "- Reference lengthscale\n"
              "- Diffusion coefficient\n"
              "--------------------------------------------------------------"
              "\n";
    return output;
  }
};

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

  template <typename Geometry, typename VelocityField>
  TransportParameters_AdvectionDiffusion_Bcc(
      Directories const &directories, std::string const &parameter_set_name,
      Geometry const &geometry, VelocityField &velocity_field) {
    std::string filename = directories.dir_parameters +
                           "/parameters_transport_" + parameter_set_name +
                           ".dat";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input);
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

    split_line = io::split_line(input);
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
    } else if (peclet_option == "compute_from_diff_coeff") {
      io::read(split_line, param_index,
               in_file + for_peclet_option +
                   "Could not parse diffusion coefficient",
               diff_coeff);
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
    } else if (peclet_option == "compute_from_diff_time") {
      io::read(split_line, param_index,
               in_file + for_peclet_option + "Could not parse diffusion time",
               diffusion_time);
      diff_coeff = lengthscale * lengthscale / (2. * diffusion_time);
    } else if (peclet_option == "set_diff_coeff") {
      io::read(split_line, param_index,
               in_file + for_peclet_option + "Could not parse Peclet number",
               peclet);
    } else {
      throw std::runtime_error{in_file + for_peclet_option + "Not supported"};
    }

    rescale(velocity_field, geometry.mesh());
  }

  template <typename VelocityField, typename Mesh>
  void rescale(VelocityField &velocity_field, Mesh const &mesh) {
    double current_mean = magnitude_of_average(velocity_field.field(), mesh);
    if (peclet_option == "rescale_velocity_to_peclet" ||
        peclet_option == "rescale_velocity_to_mean" ||
        peclet_option == "rescale_velocity_to_advection_time") {
      velocity_rescaling_factor = mean_velocity / current_mean;
      velocity_field.rescale(velocity_rescaling_factor);
    } else if (peclet_option == "compute_from_diff_coeff" ||
               peclet_option == "compute_from_diff_time") {
      velocity_rescaling_factor = 1.;
      mean_velocity = current_mean;
      advection_time = lengthscale / mean_velocity;
      peclet = 2. * diffusion_time / advection_time;
    } else if (peclet_option == "set_diff_coeff") {
      velocity_rescaling_factor = 1.;
      mean_velocity = current_mean;
      diff_coeff = lengthscale * mean_velocity / peclet;
      diffusion_time = lengthscale * lengthscale / (2. * diff_coeff);
      advection_time = lengthscale / mean_velocity;
    } else {
      throw std::runtime_error{"Peclet number setting option " + peclet_option +
                               " not supported"};
    }
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << "--------------------------------------------------------------"
              "\n"
              "Transport parameters\n"
              "--------------------------------------------------------------"
              "\n"
              "- How to set the reference lengthscale:\n"
              "  - radius\n"
              "    - Equal to bead radius\n"
              "  - diameter\n"
              "    - Equal to bead diameter\n"
              "  - cell_side\n"
              "    - Equal to primitive cell side\n"
              "  - custom\n"
              "    - Equal to specified value\n"
              "    - Pass on same line:\n"
              "      - Reference lengthscale value\n"
              "- How to set the Peclet number:\n"
              "  - compute_from_diff_coeff\n"
              "    - Compute from given diffusion coefficient (do not\n"
              "      rescale velocity field)\n"
              "    - Pass on same line:\n"
              "      - Diffusion coefficient\n"
              "  - set_diff_coeff\n"
              "    - Set diffusion coefficient to impose given Peclet number\n"
              "      (do not rescale velocity field)\n"
              "    - Pass on same line:\n"
              "      - Peclet number\n"
              "  - compute_from_diff_time\n"
              "    - Compute from given diffusion time (do not rescale\n"
              "      velocity field)\n"
              "    - Pass on same line:\n"
              "      - Diffusion coefficient\n"
              "  - rescale_velocity_to_peclet\n"
              "    - Rescale velocity field according to given peclet number\n"
              "    - Pass on same line:\n"
              "      - Peclet number\n"
              "      - Diffusion coefficient\n"
              "  - rescale_velocity_to_mean\n"
              "    - Rescale according to given mean flow velocity\n"
              "    - Pass on same line:\n"
              "      - Absolute value of mean velocity vector\n"
              "      - Diffusion coefficient\n"
              "  - rescale_velocity_to_advection_time\n"
              "    - Rescale according to given advection time\n"
              "    - Pass on same line:\n"
              "      - Advection time, based on absolute value of mean\n"
              "        velocity vector\n"
              "      - Diffusion coefficient\n"
              "--------------------------------------------------------------"
              "\n";
    return output;
  }
};

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

  template <typename Geometry, typename VelocityField>
  TransportParameters_Advection_Bcc(Directories const &directories,
                                    std::string const &parameter_set_name,
                                    Geometry const &geometry,
                                    VelocityField &velocity_field) {
    std::string filename = directories.dir_parameters +
                           "/parameters_transport_" + parameter_set_name +
                           ".dat";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input);
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
      throw std::runtime_error{"Lengthscale definition option " +
                               lengthscale_option + " not supported"};
    }

    split_line = io::split_line(input);
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

    rescale(velocity_field, geometry.mesh());
  }

  template <typename VelocityField, typename Mesh>
  void rescale(VelocityField &velocity_field, Mesh const &mesh) {
    double current_mean = magnitude_of_average(velocity_field.field(), mesh);
    if (rescale_velocity_option == "rescale_velocity_to_mean" ||
        rescale_velocity_option == "rescale_velocity_to_advection_time") {
      velocity_rescaling_factor = mean_velocity / current_mean;
      velocity_field.rescale(velocity_rescaling_factor);
    } else if (rescale_velocity_option == "no_rescale_velocity") {
      velocity_rescaling_factor = 1.;
      advection_time = lengthscale / current_mean;
      mean_velocity = lengthscale / advection_time;
    } else
      throw std::runtime_error{"Flow velocity field rescaling option " +
                               rescale_velocity_option + " not supported"};
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << "--------------------------------------------------------------"
              "\n"
              "Transport parameters\n"
              "--------------------------------------------------------------"
              "\n"
              "- How to set the reference lengthscale:\n"
              "  - radius\n"
              "    - Equal to bead radius\n"
              "  - diameter\n"
              "    - Equal to bead diameter\n"
              "  - cell_side\n"
              "    - Equal to primitive cell side\n"
              "  - custom\n"
              "    - Equal to specified value\n"
              "    - Pass on same line:\n"
              "      - Reference lengthscale value\n"
              "- Whether and how to rescale velocity field:\n"
              "  - rescale_velocity_to_mean\n"
              "    - Rescale according to given mean\n"
              "    - Pass on same line:\n"
              "      - Absolute value of mean velocity vector\n"
              "  - rescale_velocity_to_advection_time\n"
              "    - Rescale according to given advection time\n"
              "    - Pass on same line:\n"
              "      - Advection time, based on absolute value of mean "
              "        velocity vector\n"
              "  - no_rescale_velocity\n"
              "    - Do not rescale\n"
              "--------------------------------------------------------------"
              "\n";
    return output;
  }
};

struct TransportParameters_Diffusion_Bcc {
public:
  std::string lengthscale_option;
  double lengthscale;
  const double diff_coeff;
  const double diffusion_time;
  const double peclet{0.};
  double mean_velocity{0.};

  double cell_side;
  std::vector<std::pair<double, double>> primitive_cell_boundaries;
  double advection_time{std::numeric_limits<double>::infinity()};

  template <typename Geometry, typename VelocityField>
  TransportParameters_Diffusion_Bcc(
      Directories const &directories, std::string const &parameter_set_name,
      Geometry const &geometry,
      VelocityField &&velocity_field = meta::Empty{}) {
    std::string filename = directories.dir_parameters +
                           "/parameters_transport_" + parameter_set_name +
                           ".dat";
    auto input = io::open_read(filename);
    std::string in_file = std::string{"In file "} + filename + " : ";

    auto split_line = io::split_line(input);
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
      throw std::runtime_error{"Lengthscale definition option " +
                               lengthscale_option + " not supported"};
    }

    split_line = io::split_line(input);
    param_index = 0;
    io::read(split_line, param_index,
             in_file + "Could not parse diffusion coefficient", diff_coeff);
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << "--------------------------------------------------------------"
              "\n"
              "Transport parameters\n"
              "--------------------------------------------------------------"
              "\n"
              "- How to set the reference lengthscale:\n"
              "  - radius\n"
              "    - Equal to bead radius\n"
              "  - diameter\n"
              "    - Equal to bead diameter\n"
              "  - cell_side\n"
              "    - Equal to primitive cell side\n"
              "  - custom\n"
              "    - Equal to specified value\n"
              "    - Pass on same line:\n"
              "      - Reference lengthscale value\n"
              "- Diffusion coefficient"
              "--------------------------------------------------------------"
              "\n";
    return output;
  }
};

template <typename Solvers_t, typename Parameters_t>
struct Transport_LinearInterp {
  Transport_LinearInterp() = delete;

  using Solvers = Solvers_t;
  using Parameters = Parameters_t;

  template <typename VelocityField, typename Geometry, typename Boundary,
            typename ReactionParameters>
  static auto
  makeTransitions(VelocityField const &velocity_field, Geometry const &geometry,
                  Boundary &boundary, Parameters const &params_transport,
                  ReactionParameters const &params_reaction,
                  typename Solvers::Parameters const &params_solvers) {
    return makeTransportTransitions<Solvers>(velocity_field, geometry, boundary,
                                             params_transport, params_reaction,
                                             params_solvers);
  }

  template <typename Geometry>
  static auto makeVelocityInterpolator(Geometry const &geometry) {
    if constexpr (Solvers::Parameters::advection)
      return makeLinearVelocityInterpolator(geometry,
                                            get_velocity_data(geometry.mesh()));
    else
      return makeLinearVelocityInterpolator(
          geometry,
          Foam::dimensionedVector{"", Foam::dimVelocity, Foam::zero{}});
  }

  template <typename Geometry, typename Uninterpolated>
  static auto makeVelocityInterpolator(Geometry const &geometry,
                                       Uninterpolated &&uninterpolated) {
    if constexpr (Solvers::Parameters::advection)
      return makeLinearVelocityInterpolator(
          geometry, get_velocity_data(geometry.mesh()),
          std::forward<Uninterpolated>(uninterpolated));
    else
      return makeLinearVelocityInterpolator(
          geometry,
          Foam::dimensionedVector{"", Foam::dimVelocity, Foam::zero{}},
          Foam::dimensionedVector{"", Foam::dimVelocity, Foam::zero{}});
  }

  template <typename Geometry>
  static auto
  makeVelocityInterpolator_WithUninterpolated(Geometry const &geometry) {
    if constexpr (Solvers::Parameters::advection)
      return makeLinearVelocityInterpolator(
          geometry, get_velocity_data(geometry.mesh()),
          Foam::volVectorField{
              Foam::IOobject{"", geometry.mesh().time().timeName(),
                             geometry.mesh(), Foam::IOobject::NO_READ,
                             Foam::IOobject::NO_WRITE,
                             Foam::IOobject::NO_REGISTER},
              geometry.mesh(),
              Foam::dimensionedVector{"", Foam::dimVelocity, Foam::zero{}}});
    else
      return makeLinearVelocityInterpolator(
          geometry,
          Foam::dimensionedVector{"", Foam::dimVelocity, Foam::zero{}},
          Foam::dimensionedVector{"", Foam::dimVelocity, Foam::zero{}});
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    output << "--------------------------------------------------------------\n"
              "Transport\n"
              "--------------------------------------------------------------\n"
              "Advection: ";
    if constexpr (Solvers::Parameters::advection) {
      output << "Yes\n";
    } else {
      output << "No\n";
    }
    output << "Diffusion: ";
    if constexpr (Solvers::Parameters::diffusion) {
      output << "Yes\n";
    } else {
      output << "No\n";
    }
    if constexpr (Solvers::Parameters::advection) {
      output << "Velocity interpolation: Linear\n";
    }
    output
        << "--------------------------------------------------------------\n";
    return output;
  }
};

template <typename Solvers>
using Transport_Generic = Transport_LinearInterp<
    Solvers,
    std::conditional_t<Solvers::Parameters::advection &&
                           Solvers::Parameters::advection,
                       TransportParameters_AdvectionDiffusion,
                       std::conditional_t<Solvers::Parameters::advection,
                                          TransportParameters_Advection,
                                          TransportParameters_Diffusion>>>;

template <typename Solvers>
using Transport_Bcc = Transport_LinearInterp<
    Solvers,
    std::conditional_t<Solvers::Parameters::advection &&
                           Solvers::Parameters::advection,
                       TransportParameters_AdvectionDiffusion_Bcc,
                       std::conditional_t<Solvers::Parameters::advection,
                                          TransportParameters_Advection_Bcc,
                                          TransportParameters_Diffusion_Bcc>>>;

} // namespace ptof

#endif /* PTOF_TRANSPORT_H */
