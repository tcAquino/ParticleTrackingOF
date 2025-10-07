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
#include "PTOF/Field.h"
#include "PTOF/Useful.h"
#include "vector.H"
#include "volFields.H"
#include "word.H"
#include <IOobject.H>
#include <cmath>
#include <fieldTypes.H>
#include <fvcGrad.H>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <volFields.H>

namespace ptof {
/**
   \struct Phase PTOF/Phase.h "PTOF/Phase.h"
   \brief Phase field information.
*/
struct Phase {
  /** \brief Deleted constructor. */
  Phase() = delete;

  /**
     \struct Phase::Parameters PTOF/Phase.h "PTOF/Phase.h"
     \brief Parameters for phase-related quantities.
  */
  struct Parameters {
  public:
    std::string
        phase_name; /**< Name of phase saturation field to be read from file. */
    bool carrier_phase;    /**< Weather named phase is carrier phase forl steep
                              chemical potential gradient. */
    bool compute_gradient; /**< Weather to compute or read phase gradient from
                              file. */
    double chemical_potential_difference; /**< Chemical potential difference
                                 between excluded and carrier phases at
                                 equilibrium (minus log of equilibrium ratio of
                                 concentrations) */
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

      auto split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
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
      } else {
        throw std::invalid_argument{in_file + for_carrier_option +
                                    "Not supported"};
      }

      split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
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
      } else {
        throw std::invalid_argument{in_file + gradient_option +
                                    "Not supported"};
      }

      split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
      param_index = 0;
      auto phase_concentration_ratio = io::read<double>(
          split_line, param_index,
          in_file + "Could not parse equilibrium phase concentration ratio");
      chemical_potential_difference = -std::log(phase_concentration_ratio);

      split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
      param_index = 0;
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
                "  - Whether named phase is carrier or excluded phase (in the\n"
                "    limit of steep chemical potential gradient)\n"
                "    - carrier\n"
                "      - Transport in this phase\n"
                "    - excluded\n"
                "      - No transport in this phase\n"
                "- How to compute gradient of carrier phase:\n"
                "  - read\n"
                "    - Based on reading gradient of named phase\n"
                "  - compute\n"
                "    - Computed from named phase data\n"
                "- Equilibrium ratio of concentrations between excluded phase\n"
                "  and carrier phase\n"
                "- Phase field tolerance to consider pure phase\n"
             << io::line();
      return output;
    }
  };

  /**
     \return Carrier phase saturation field.
     \details Obtained from OpenFOAM file data for either the carrier phase or
     the carrier phase.
  */
  template <typename Mesh>
  static auto get_carrier_phase_data(Mesh const &mesh,
                                     Foam::word const &timeName,
                                     Parameters const &params_phase) {
    if (params_phase.carrier_phase) {
      return Foam::volScalarField{
          Foam::IOobject{std::string("alpha.") + params_phase.phase_name,
                         timeName, mesh, Foam::IOobject::MUST_READ,
                         Foam::IOobject::NO_WRITE},
          mesh};
    } else {
      return static_cast<Foam::volScalarField>(
          1. -
          Foam::volScalarField{
              Foam::IOobject{std::string("alpha.") + params_phase.phase_name,
                             timeName, mesh, Foam::IOobject::MUST_READ,
                             Foam::IOobject::NO_WRITE},
              mesh});
    }
  }

  /**
     \return Carrier phase saturation field.
     \details Obtained from OpenFOAM file data for either the carrier phase or
     the carrier phase.
  */
  template <typename Mesh>
  static auto get_carrier_phase_data(Mesh const &mesh, Foam::scalar time,
                                     Parameters const &params_phase) {
    return get_carrier_phase_data(mesh, closest_time_name(mesh, time),
                                  params_phase);
  }

  /**
     \return Carrier phase saturation field.
     \details Obtained from OpenFOAM file data for either the carrier phase or
     the carrier phase.
  */
  template <typename Mesh>
  static auto get_carrier_phase_data(Mesh const &mesh,
                                     Parameters const &params_phase) {
    return get_carrier_phase_data(mesh, mesh.time().timeName(), params_phase);
  }

  /**
     \return Carrier phase field gradient.
     \details Based on OpenFOAM file data.
  */
  template <typename Mesh>
  static auto get_grad_carrier_phase_data(Mesh const &mesh,
                                          Foam::word const &timeName,
                                          Parameters const &params_phase) {
    if (params_phase.carrier_phase) {
      return Foam::volVectorField{
          Foam::IOobject{std::string("gradAlpha.") + params_phase.phase_name,
                         timeName, mesh, Foam::IOobject::MUST_READ,
                         Foam::IOobject::NO_WRITE},
          mesh};
    } else {
      return static_cast<Foam::volVectorField>(-Foam::volVectorField{
          Foam::IOobject{std::string("gradAlpha.") + params_phase.phase_name,
                         timeName, mesh, Foam::IOobject::MUST_READ,
                         Foam::IOobject::NO_WRITE},
          mesh});
    }
  }

  /**
     \return Carrier phase field gradient.
     \details Based on OpenFOAM file data.
  */
  template <typename Mesh>
  static auto get_grad_carrier_phase_data(Mesh const &mesh,
                                          Parameters const &params_phase) {
    return get_grad_carrier_phase_data(mesh, mesh.time().timeName(),
                                       params_phase);
  }

  /**
     \return Carrier phase field gradient.
     \details Based on OpenFOAM file data.
  */
  template <typename Mesh>
  static auto get_grad_carrier_phase_data(Mesh const &mesh, Foam::scalar time,
                                          Parameters const &params_phase) {
    return get_grad_carrier_phase_data(mesh, closest_time_name(mesh, time),
                                       params_phase);
  }

  /**
     \return Carrier phase field gradient.
     \details Based on OpenFOAM file data or computed from carrier phase
     saturation field.
  */
  template <typename Mesh>
  static auto
  grad_carrier_phase(Mesh const &mesh,
                     Foam::volScalarField const &carrier_phase_field,
                     Parameters const &params_phase) {
    if (params_phase.compute_gradient) {
      return static_cast<Foam::volVectorField>(
          Foam::fvc::grad(carrier_phase_field));
    } else {
      return get_grad_carrier_phase_data(mesh, params_phase);
    }
  }

  /**
     \return Carrier phase field gradient.
     \details Based on OpenFOAM file data or computed from carrier phase
     saturation field.
  */
  template <typename Mesh>
  static auto
  grad_carrier_phase(Mesh const &mesh, Foam::word const &timeName,
                     Foam::volScalarField const &carrier_phase_field,
                     Parameters const &params_phase) {
    if (params_phase.compute_gradient) {
      return static_cast<Foam::volVectorField>(
          Foam::fvc::grad(carrier_phase_field));
    } else {
      return get_grad_carrier_phase_data(mesh, timeName, params_phase);
    }
  }

  /**
     \return Carrier phase field gradient.
     \details Based on OpenFOAM file data or computed from carrier phase
     saturation field.
  */
  template <typename Mesh>
  static auto
  grad_carrier_phase(Mesh const &mesh, Foam::scalar time,
                     Foam::volScalarField const &carrier_phase_field,
                     Parameters const &params_phase) {
    if (params_phase.compute_gradient) {
      return static_cast<Foam::volVectorField>(
          Foam::fvc::grad(carrier_phase_field));
    } else {
      return get_grad_carrier_phase_data(mesh, closest_time_name(mesh, time),
                                         params_phase);
    }
  }

  /**
     \brief Update phase field and velocity field.
  */
  template <typename PhaseField, typename VelocityField, typename Geometry,
            typename TransportParameters>
  static void update_phase_and_velocity_field(
      PhaseField &carrier_phase_field, VelocityField &velocity_field,
      Geometry const &geometry, TransportParameters const &params_transport,
      Parameters const &params_phase) {
    update_phase_field(carrier_phase_field, geometry, params_phase);
    update_velocity_field(velocity_field, geometry, params_transport);
    update_effective_velocity_field(carrier_phase_field, velocity_field,
                                    geometry, params_transport, params_phase);
  }

  /**
     \brief Make phase field.
  */
  template <typename Geometry, typename PhaseFieldData>
  static auto makeCarrierPhase(Geometry const &geometry,
                               PhaseFieldData &&carrier_phase_data) {
    return ptof::ScalarField_Interpolation{
        std::forward<PhaseFieldData>(carrier_phase_data), geometry.locator,
        InterpolationTypes::Cell{}, CheckOptions::Check{}};
  }

  template <typename Geometry, typename PhaseFieldData>
  static auto makeCarrierPhaseTimeInterpolator(
      Geometry const &geometry, PhaseFieldData &&carrier_phase_data_new,
      PhaseFieldData &&carrier_phase_data_old, Foam::scalar time_new,
      Foam::scalar time_old) {

    using Type = decltype(makeCarrierPhase(
        geometry, std::forward<PhaseFieldData>(carrier_phase_data_new)));
    return ptof::Field_TimeInterpolation{
        std::unique_ptr<Type>{new Type{makeCarrierPhase(
            geometry, std::forward<PhaseFieldData>(carrier_phase_data_new))}},
        std::unique_ptr<Type>{new Type{makeCarrierPhase(
            geometry, std::forward<PhaseFieldData>(carrier_phase_data_old))}},
        time_new, time_old};
  }

  template <typename VelocityField, typename PhaseField>
  class EffectiveVelocity {
  public:
    using Point = Foam::point;                     /**> 3D point. */
    using Point2D = Foam::Vector2D<Foam::scalar>;  /**> 2D point. */
    using Vector = Foam::vector;                   /**> 3D vector. */
    using Vector2D = Foam::Vector2D<Foam::scalar>; /**> 2D vector. */
    using Scalar = Foam::scalar; /**> Scalar (also 1D point). */
    using Index = Foam::label;   /**> Cell index. */
    using Time = Foam::scalar;   /**> Time type. */

    template <typename TransportParameters>
    EffectiveVelocity(VelocityField &&velocity_field,
                      PhaseField &&carrier_phase_field,
                      TransportParameters const &params_transport,
                      Parameters const &params_phase)
        : _velocity_field(std::forward<VelocityField>(velocity_field)),
          _carrier_phase_field(std::forward<PhaseField>(carrier_phase_field)),
          _effective_drift{Field_TimeInterpolation{
              std::make_unique<DriftVelocityField>(
                  effective_drift_field(
                      closest_time_name(_carrier_phase_field.locator().mesh(),
                                        _carrier_phase_field.time_new()),
                      _carrier_phase_field.field_new(), params_transport,
                      params_phase),
                  carrier_phase_field.locator(), InterpolationTypes::Cell{},
                  CheckOptions::Check{}),
              std::make_unique<DriftVelocityField>(
                  effective_drift_field(
                      closest_time_name(_carrier_phase_field.locator().mesh(),
                                        _carrier_phase_field.time_old()),
                      _carrier_phase_field.field_old(), params_transport,
                      params_phase),
                  _carrier_phase_field.locator(), InterpolationTypes::Cell{},
                  CheckOptions::Check{}),
              _carrier_phase_field.time_new(),
              _carrier_phase_field.time_old()}} {}

    /**
      \brief Evaluate field.
      \param position Position.
      \param cell Mesh cell index position is in.
      \param time Interpolation time.
      \return Field value.
   */
    template <typename Position>
    auto operator()(Position const &position, Index cell, Time time) const {
      if constexpr (!std::is_same_v<std::remove_reference_t<
                                        std::remove_const_t<VelocityField>>,
                                    meta::Empty>) {
        std::cout << "OLA2!!!";
        io::print(std::cout, _effective_drift(position, cell, time))
            << std::endl;
        return _velocity_field(position, cell, time) +
               _effective_drift(position, cell, time);
      } else {
        std::cout << "OLA!!!";
        io::print(std::cout, _effective_drift(position, cell, time)) << std::endl;
        return _effective_drift(position, cell, time);
      }
    }

    /**
       \brief Field value.
       \param state Particle state to interpolate.
    */
    template <typename State> auto operator()(State const &state) const {
      return (*this)(state.position, state.cell, state.time);
    }

    /**
       \brief Evaluate field.
       \param position 3D position.
       \param time Interpolation time.
       \return Field value.
    */
    auto operator()(Point const &position, Time time) const {
      return (*this)(position, locate(position), time);
    }

    /**
       \brief Interpolate field.
       \param position 2D position.
       \param time Interpolation time.
       \return Field value.
    */
    auto operator()(Point2D const &position, Time time) const {
      return (*this)(position, locate(position), time);
    }

    /**
       \brief Evaluate field.
       \param position 1D position.
       \param time Interpolation time.
       \return Field value.
    */
    auto operator()(Scalar const &position, Time time) const {
      return (*this)(position, locate(position), time);
    }

    /** \brief Field at cell center. */
    auto operator()(Foam::label cell_id, Time time) {
      if constexpr (!std::is_same_v<std::remove_reference_t<
                                        std::remove_const_t<VelocityField>>,
                                    meta::Empty>) {
        return _velocity_field(cell_id, time) + _effective_drift(cell_id, time);
      } else {
        return _effective_drift(cell_id, time);
      }
    }

    /**
       \brief Locate a state in the mesh.
       \param state Particle state to locate.
       \return Mesh cell index.
    */
    template <typename State> auto locate(State const &state) const {
      return _carrier_phase_field.locate(state);
    }

    /** \return Locator object to find positions in mesh. */
    auto const &locator() const { return _carrier_phase_field->locator(); }

    /** \return Full field of underlying field type. */
    auto field(Time time) const {
      return static_cast<Foam::volVectorField>(_velocity_field.field(time) +
                                               _effective_drift.field(time));
    }

    /** \return Underlying boundary field. */
    auto boundaryField(Time time) const {
      using BoundaryField = decltype(_effective_drift.boundaryField(time));
      if constexpr (!std::is_same_v<std::remove_reference_t<
                                        std::remove_const_t<VelocityField>>,
                                    meta::Empty>) {
        return static_cast<BoundaryField>(_velocity_field.boundaryField(time) +
                                          _effective_drift.boundaryField(time));
      } else {
        return static_cast<BoundaryField>(_effective_drift.boundaryField(time));
      }
    }

    /** \return Underlying boundary field at patch. */
    auto boundaryField(Foam::label patch_id, Time time) const {
      using BoundaryField =
          decltype(_effective_drift.boundaryField(patch_id, time));
      if constexpr (!std::is_same_v<std::remove_reference_t<
                                        std::remove_const_t<VelocityField>>,
                                    meta::Empty>) {
        return static_cast<BoundaryField>(
            _velocity_field.boundaryField(patch_id, time) +
            _effective_drift.boundaryField(patch_id, time));
      } else {
        return static_cast<BoundaryField>(
            _effective_drift.boundaryField(patch_id, time));
      }
    }

    /**
       \return Underlying boundary field value at face of patch.
       \note Face id in patch is numbered relative to the patch, starting from
       0.
    */
    auto boundaryField(Foam::label patch_id, Foam::label face_in_patch,
                       Time time) const {
      using BoundaryField = decltype(_effective_drift.boundaryField(
          patch_id, face_in_patch, time));
      if constexpr (!std::is_same_v<std::remove_reference_t<
                                        std::remove_const_t<VelocityField>>,
                                    meta::Empty>) {
        return static_cast<BoundaryField>(
            _velocity_field.boundaryField(patch_id, face_in_patch, time) +
            _effective_drift.boundaryField(patch_id, face_in_patch, time));
      } else {
        return static_cast<BoundaryField>(
            _effective_drift.boundaryField(patch_id, face_in_patch, time));
      }
    }

    template <typename Geometry, typename TransportParameters>
    void update(Geometry const &geometry, Foam::instant const &instant,
                TransportParameters const &params_transport,
                Parameters const &params_phase) {
      if constexpr (!std::is_same_v<std::remove_reference_t<
                                        std::remove_const_t<VelocityField>>,
                                    meta::Empty>) {
        _velocity_field.update(
            get_velocity_data(geometry.mesh(), instant.name()),
            instant.value());
      }
      _carrier_phase_field.update(Phase::get_carrier_phase_data(geometry.mesh(),
                                                                instant.name(),
                                                                params_phase),
                                  instant.value());
      _effective_drift.update(
          effective_drift_field(instant.name(),
                                _carrier_phase_field.field_new(),
                                params_transport, params_phase),
          instant.value());
    }

  private:
    VelocityField _velocity_field;
    PhaseField _carrier_phase_field;
    using DriftVelocityField = ptof::VectorField_Interpolation<
        Foam::volVectorField, decltype(_carrier_phase_field.locator()),
        InterpolationTypes::Cell, CheckOptions::Check>;
    using DriftField = Field_TimeInterpolation<DriftVelocityField>;
    DriftField _effective_drift;

    template <typename TransportParameters>
    auto
    effective_drift_coefficient(TransportParameters const &params_transport,
                                Parameters const &params_phase) {
      return Foam::dimensionedScalar("", Foam::dimViscosity,
                                     params_transport.diff_coeff) *
             params_phase.chemical_potential_difference;
    }

    template <typename TransportParameters, typename PhaseFieldSpatialInterp>
    auto
    effective_drift_field(Foam::word const &timeName,
                          PhaseFieldSpatialInterp const &carrier_phase_field,
                          TransportParameters const &params_transport,
                          Parameters const &params_phase) {
      return static_cast<Foam::volVectorField>(
          effective_drift_coefficient(params_transport, params_phase) *
          grad_carrier_phase(carrier_phase_field.locator().mesh(), timeName,
                             carrier_phase_field.field(), params_phase));
    }
  };
  template <typename VelocityField, typename PhaseField,
            typename TransportParameters>
  EffectiveVelocity(VelocityField &&, PhaseField &&,
                    TransportParameters const &, Parameters const &)
      -> EffectiveVelocity<VelocityField, PhaseField>;
};
} // namespace ptof

#endif /* PTOF_PHASE_H */
