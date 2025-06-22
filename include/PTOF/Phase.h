/**
   \file PTOF/Phase.h
   \author Tomas Aquino
   \date 08/02/2024
   \brief Definitions and utilities to handle VOF-type phase fields (saturation)
   and transport in two-phase flow.
*/

#ifndef PTOF_PHASE_H
#define PTOF_PHASE_H

#include "General/IO.h"
#include "PTOF/Advection.h"
#include "PTOF/Directories.h"
#include <IOobject.H>
#include <cmath>
#include <fieldTypes.H>
#include <fvcGrad.H>
#include <ostream>
#include <stdexcept>
#include <string>
#include <volFields.H>

namespace ptof {
/**
   \struct Phase PTOF/Phase.h "PTOF/Phase.h"
   \brief Phase field information.
*/
struct Phase {
  /** \brief Deleted constructor. */
  Phase() = delete;

  using PhaseField = Foam::volScalarField; /**< Saturation scalar field as a
                                              function of mesh cell index. */
  using GradPhaseField =
      Foam::volVectorField; /**< Gradient of saturation vector field as a
                               function of mesh cell index. */

  /**
     \struct Phase::Parameters PTOF/Phase.h "PTOF/Phase.h"
     \brief Parameters for phase-related quantities.
  */
  struct Parameters {
  public:
    std::string
        phase_name; /**< Name of phase saturation field to be read from file. */
    bool carrier_phase;    /**< Weather named phase is carrier phase. */
    bool compute_gradient; /**< Weather to compute or read phase gradient from
                              file. */
    double leakage_coefficient; /**< Log of ratio of concentrations between
                                   excluded and carrier phase at
                                   equilibrium. */
    double phase_tolerance; /**< Tolerance in local phase saturation to consider
                               pure phase. */

    /**
       \brief Constructor.
       \param directories Current case directory information.
       \param name Name of phase parameters set.
       \param geometry Domain geometry info and utilities.
    */
    template <typename Geometry>
    Parameters(Directories const &directories, std::string const &name,
               Geometry const &geometry) {
      std::string filename =
          directories.dir_parameters + "/parameters_phase_" + name + ".dat";
      auto input = io::open_read(filename);
      std::string in_file = std::string{"In file "} + filename + " : ";

      auto split_line = io::split_line(input);
      std::size_t param_index = 0;
      io::read(split_line, param_index, in_file + "Could not parse phase name",
               phase_name);
      auto carrier_option = io::read<std::string>(
          split_line, param_index,
          in_file + "Could not parse carrier or excluded phase option");
      std::string for_carrier_option =
          std::string{"Carrier or excluded phase option "} + carrier_option +
          " : ";
      if (carrier_option == "carrier") {
        carrier_phase = true;
      } else if (carrier_option == "excluded") {
        carrier_phase = false;
      } else
        throw std::invalid_argument{in_file + for_carrier_option +
                                    "Not supported"};

      split_line = io::split_line(input);
      param_index = 0;
      auto gradient_option = io::read<std::string>(
          split_line, param_index,
          in_file + "Could not parse compute or read phase gradient option");
      std::string for_gradient_option =
          std::string{"Compute or read phase gradient option "} +
          gradient_option + " : ";
      if (gradient_option == "compute") {
        compute_gradient = true;
      } else if (gradient_option == "read") {
        compute_gradient = false;
      } else
        throw std::invalid_argument{in_file + gradient_option +
                                    "Not supported"};

      split_line = io::split_line(input);
      param_index = 0;
      auto leakage_tolerance =
          io::read<double>(split_line, param_index,
                           in_file + "Could not parse leakage tolerance");
      leakage_coefficient = std::log(leakage_tolerance);
      io::read(split_line, param_index,
               in_file + "Could not parse pure phase tolerance",
               phase_tolerance);
    }

    /** \brief Output general information about object. */
    inline static std::ostream &info(std::ostream &output) {
      output << io::line() << "Phase parameters\n"
             << io::line()
             << "- Name of phase to be read from file\n"
                "- Pass on same line:\n"
                "  - Whether named phase is carrier or excluded phase:\n"
                "    - carrier\n"
                "      - Transport in this phase\n"
                "    - excluded\n"
                "      - No transport in this phase\n"
                "- How to compute gradient of carrier phase:\n"
                "  - read\n"
                "    - Based on reading gradient of named phase\n"
                "  - compute\n"
                "    - Computed from named phase data\n"
                "- Leakage tolerance (log of ratio of concentrations between\n"
                "  excluded phase and carrier phase)\n"
                "- Phase field tolerance to consider pure phase\n"
             << io::line();
      return output;
    }
  };

  /**
     \return Carrier phase saturation field.
     \details Obtained from OpenFOAM file data for either the carrier phase or
     the carrier phase. */
  template <typename Mesh>
  static PhaseField get_carrier_phase_data(Mesh const &mesh,
                                           Parameters const &parameters) {
    if (parameters.carrier_phase)
      return PhaseField{
          Foam::IOobject{std::string("alpha.") + parameters.phase_name,
                         mesh.time().timeName(), mesh,
                         Foam::IOobject::MUST_READ, Foam::IOobject::NO_WRITE},
          mesh};
    else
      return static_cast<PhaseField>(
          1. - PhaseField{
                   Foam::IOobject{std::string("alpha.") + parameters.phase_name,
                                  mesh.time().timeName(), mesh,
                                  Foam::IOobject::MUST_READ,
                                  Foam::IOobject::NO_WRITE},
                   mesh});
  }

  /**
     \return Carrier phase field gradient.
     \details Based on OpenFOAM file data.
  */
  template <typename Mesh>
  static auto get_grad_carrier_phase_data(Mesh const &mesh,
                                          Parameters const &parameters) {
    if (parameters.carrier_phase)
      return GradPhaseField{
          Foam::IOobject{std::string("gradAlpha.") + parameters.phase_name,
                         mesh.time().timeName(), mesh,
                         Foam::IOobject::MUST_READ, Foam::IOobject::NO_WRITE},
          mesh};
    else
      return static_cast<GradPhaseField>(-GradPhaseField{
          Foam::IOobject{std::string("gradAlpha.") + parameters.phase_name,
                         mesh.time().timeName(), mesh,
                         Foam::IOobject::MUST_READ, Foam::IOobject::NO_WRITE},
          mesh});
  }

  /**
     \return Carrier phase field gradient.
     \details Based on OpenFOAM file data or computed from carrier phase
     saturation field.
  */
  template <typename Mesh>
  static auto grad_carrier_phase(Mesh const &mesh, Parameters const &parameters,
                                 PhaseField const &carrier_phase_field) {
    if (parameters.compute_gradient)
      return static_cast<GradPhaseField>(Foam::fvc::grad(carrier_phase_field));
    else
      return get_grad_carrier_phase_data(mesh, parameters);
  }

  /**
     \brief Update phase field
  */
  template <typename Geometry>
  static void update_phase_field(PhaseField &carrier_phase_field,
                                 Geometry const &geometry,
                                 Parameters const &params_phase) {
    carrier_phase_field =
        Phase::get_carrier_phase_data(geometry.mesh(), params_phase);
  }

  /**
     \brief Update effective velocity field
  */
  template <typename VelocityField, typename Geometry,
            typename TransportParameters>
  static void update_effective_velocity_field(
      PhaseField &carrier_phase_field, VelocityField &velocity_field,
      Geometry const &geometry, TransportParameters const &params_transport,
      Parameters const &params_phase) {
    velocity_field.set_uninterpolated(
        Foam::dimensionedScalar("", Foam::dimViscosity,
                                params_transport.diff_coeff) *
        params_phase.leakage_coefficient *
        Phase::grad_carrier_phase(geometry.mesh(), params_phase,
                                  carrier_phase_field));
  }

  /**
     \brief Update phase field and velocity field
  */
  template <typename VelocityField, typename Geometry,
            typename TransportParameters>
  static void update_phase_and_velocity_field(
      PhaseField &carrier_phase_field, VelocityField &velocity_field,
      Geometry const &geometry, TransportParameters const &params_transport,
      Parameters const &params_phase) {
    update_phase_field(carrier_phase_field, geometry, params_phase);
    ptof::update_velocity_field(velocity_field, geometry, params_transport);
    update_effective_velocity_field(carrier_phase_field, velocity_field,
                                    geometry, params_transport, params_phase);
  }

  /**
     \brief Make phase field
  */
  template <typename Geometry>
  static PhaseField makePhaseField(Geometry const &geometry,
                                   Parameters const &params_phase) {
    return get_carrier_phase_data(geometry.mesh(), params_phase);
  }

  /**
     \brief Make phase field and update effective velocity field
  */
  template <typename VelocityField, typename Geometry,
            typename TransportParameters>
  static PhaseField
  makePhaseField(VelocityField &velocity_field, Geometry const &geometry,
                 Parameters const &params_phase,
                 TransportParameters const &params_transport) {
    PhaseField carrier_phase_field =
        get_carrier_phase_data(geometry.mesh(), params_phase);
    update_effective_velocity_field(carrier_phase_field, velocity_field,
                                    geometry, params_transport, params_phase);
    return carrier_phase_field;
  }
};
} // namespace ptof

#endif /* PTOF_PHASE_H */
