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
#include "General/Meta.h"
#include "General/Operation.h"
#include "PTOF/Advection.h"
#include "PTOF/CheckOptions.h"
#include "PTOF/Directories.h"
#include "PTOF/Field.h"
#include "PTOF/Useful.h"
#include "Roots/Roots.h"
#include "instant.H"
#include "vector.H"
#include "volFieldsFwd.H"
#include "word.H"
#include <IOobject.H>
#include <cmath>
#include <fieldTypes.H>
#include <fvcGrad.H>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <volFields.H>

namespace ptof {
/** \struct */
struct ChemicalPotentialModels {
  struct AGG;
  struct JM;
  struct None;
};

/**
   \class Phase PTOF/Phase.h "PTOF/Phase.h"
   \brief Phase field information.
*/
template <typename ChemicalPotentialModel,
          typename TimeInterpolationType = InterpolationTypes::Linear,
          bool hard_reflection = false>
class Phase {
public:
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
    double phase_concentration_ratio{1.}; /**< Equilibrium ratio of
                                             concentrations in excluded and
                                             carrier phases at equilibrium. */
    double chemical_potential_difference{-std::log(
        phase_concentration_ratio)}; /**< Chemical potential
                                        difference between excluded and carrier
                                        phases at equilibrium (minus log of
                                        equilibrium ratio of concentrations). */
    double phase_tolerance;          /**< Tolerance in local phase saturation to
                                        consider pure phase. */

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
          directories.dir_parameters + "/phase_" + name + ".param";
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

      if constexpr (!std::is_same_v<ChemicalPotentialModel,
                                    ChemicalPotentialModels::None>) {
        split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
        param_index = 0;
        io::read<double>(
            split_line, param_index,
            in_file + "Could not parse equilibrium phase concentration ratio",
            phase_concentration_ratio);
        chemical_potential_difference = -std::log(phase_concentration_ratio);
      }

      split_line = io::split_line(input, "#", "\t,|\r()[]{} ");
      param_index = 0;
      io::read(split_line, param_index,
               in_file + "Could not parse pure phase tolerance",
               phase_tolerance);
      if (phase_tolerance > 1. || phase_tolerance < 0.) {
        std::runtime_error{in_file + "Pure phase tolerance : " + "Value " +
                           std::to_string(phase_tolerance) +
                           " is not between 0 and 1"};
      }
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
                "    - Computed from named phase data\n";
      if constexpr (!std::is_same_v<ChemicalPotentialModel,
                                    ChemicalPotentialModels::None>) {
        output
            << "- Equilibrium ratio of concentrations between excluded phase\n"
               "  and carrier phase\n";
      }
      output << "- Phase field tolerance to consider pure phase\n"
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
        time_new, time_old, TimeInterpolationType{}};
  }

  template <typename EffectiveVelocityField, typename PhaseField,
            typename TransportParameters>
  static void update_effective_velocity_and_phase_field(
      EffectiveVelocityField &velocity_field, PhaseField &carrier_phase_field,
      Foam::instant const &instant, TransportParameters const &params_transport,
      Parameters const &params_phase) {
    if constexpr (!std::is_same_v<EffectiveVelocityField, meta::Empty>) {
      velocity_field.update(instant, params_transport, params_phase);
    } else {
      carrier_phase_field.update(
          Phase::get_carrier_phase_data(carrier_phase_field.locator().mesh(),
                                        instant.name(), params_phase),
          instant.value());
    }
  }

  // TransportParameters should be a template parameter of some member functions
  // only, but for some reason this crashes some versions of clangd. Retry
  // later.
  template <typename VelocityField, typename PhaseField,
            typename TransportParameters>
  class EffectiveVelocity {
  public:
    using Point = Foam::point;                     /**< 3D point. */
    using Point2D = Foam::Vector2D<Foam::scalar>;  /**< 2D point. */
    using Vector = Foam::vector;                   /**< 3D vector. */
    using Vector2D = Foam::Vector2D<Foam::scalar>; /**< 2D vector. */
    using Scalar = Foam::scalar; /**< Scalar (also 1D point). */
    using Index = Foam::label;   /**< Cell index. */
    using Time = Foam::scalar;   /**< Time type. */

    static_assert(
        std::is_same_v<ChemicalPotentialModel, ChemicalPotentialModels::AGG> ||
            std::is_same_v<ChemicalPotentialModel,
                           ChemicalPotentialModels::JM> ||
            std::is_same_v<ChemicalPotentialModel,
                           ChemicalPotentialModels::None>,
        "Phase : Effective drift: Only chemical potential models "
        "AGG, JM, or None are supported");

    static constexpr bool has_base_velocity = !std::is_same_v<
        std::remove_reference_t<std::remove_const_t<VelocityField>>,
        meta::Empty>;
    static constexpr bool has_chemical_potential =
        !std::is_same_v<ChemicalPotentialModel, ChemicalPotentialModels::None>;

    EffectiveVelocity(VelocityField &&velocity_field,
                      PhaseField &&carrier_phase_field,
                      TransportParameters const &params_transport,
                      Parameters const &params_phase)
        : _velocity_field(std::forward<VelocityField>(velocity_field)),
          _carrier_phase_field(std::forward<PhaseField>(carrier_phase_field)),
          _effective_drift{Field_TimeInterpolation{
              std::make_unique<DriftVelocityField>(
                  effective_drift_field(
                      closest_instant(_carrier_phase_field.locator().mesh(),
                                      _carrier_phase_field.time_new()),
                      _carrier_phase_field.field_new(), params_transport,
                      params_phase),
                  carrier_phase_field.locator(), InterpolationTypes::Cell{},
                  CheckOptions::Check{}),
              std::make_unique<DriftVelocityField>(
                  effective_drift_field(
                      closest_instant(_carrier_phase_field.locator().mesh(),
                                      _carrier_phase_field.time_old()),
                      _carrier_phase_field.field_old(), params_transport,
                      params_phase),
                  _carrier_phase_field.locator(), InterpolationTypes::Cell{},
                  CheckOptions::Check{}),
              _carrier_phase_field.time_new(), _carrier_phase_field.time_old(),
              TimeInterpolationType{}}} {}

    /**
       \brief Evaluate field.n
       \param position Position.
       \param cell Mesh cell index position is in.
       \param time Interpolation time.
       \return Field value.
    */
    template <typename Position>
    auto operator()(Position const &position, Index cell, Time time) const {
      if constexpr (has_base_velocity && has_chemical_potential) {
        return _velocity_field(position, cell, time) +
               _effective_drift(position, cell, time);
      }
      if constexpr (has_base_velocity) {
        return _velocity_field(position, cell, time);
      }
      return _effective_drift(position, cell, time);
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
      if constexpr (has_base_velocity && has_chemical_potential) {
        return _velocity_field(cell_id, time) + _effective_drift(cell_id, time);
      }
      if constexpr (has_base_velocity) {
        return _velocity_field(cell_id, time);
      }
      return _effective_drift(cell_id, time);
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
      if constexpr (has_base_velocity && has_chemical_potential) {
        return static_cast<Foam::volVectorField>(_velocity_field.field(time) +
                                                 _effective_drift.field(time));
      }
      if constexpr (has_base_velocity) {
        return _velocity_field.field(time);
      }
      return _effective_drift.field(time);
    }

    /** \return Underlying boundary field. */
    auto boundaryField(Time time) const {
      using BoundaryField = decltype(_effective_drift.boundaryField(time));
      if constexpr (has_base_velocity && has_chemical_potential) {
        return static_cast<BoundaryField>(_velocity_field.boundaryField(time) +
                                          _effective_drift.boundaryField(time));
      }
      if constexpr (has_base_velocity) {
        _velocity_field.boundaryField(time);
      }
      return _effective_drift.boundaryField(time);
    }

    /** \return Underlying boundary field at patch. */
    auto boundaryField(Foam::label patch_id, Time time) const {
      using BoundaryField =
          decltype(_effective_drift.boundaryField(patch_id, time));
      if constexpr (has_base_velocity && has_chemical_potential) {
        return static_cast<BoundaryField>(
            _velocity_field.boundaryField(patch_id, time) +
            _effective_drift.boundaryField(patch_id, time));
      }
      if constexpr (has_base_velocity) {
        return _velocity_field.boundaryField(patch_id, time);
      }
      return _effective_drift.boundaryField(patch_id, time);
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
      if constexpr (has_base_velocity && has_chemical_potential) {
        return static_cast<BoundaryField>(
            _velocity_field.boundaryField(patch_id, face_in_patch, time) +
            _effective_drift.boundaryField(patch_id, face_in_patch, time));
      }
      if constexpr (has_base_velocity) {
        return _velocity_field.boundaryField(patch_id, face_in_patch, time);
      }
      return _effective_drift.boundaryField(patch_id, face_in_patch, time);
    }

    void update(Foam::instant const &instant,
                TransportParameters const &params_transport,
                Parameters const &params_phase) {
      if constexpr (!has_chemical_potential) {
        return;
      }
      auto const &mesh = _carrier_phase_field.locator().mesh();
      if constexpr (has_base_velocity) {
        auto velocity_data = get_velocity_data(mesh, instant.name());
        rescale(velocity_data, params_transport.velocity_rescaling_factor);
        _velocity_field.update(std::move(velocity_data), instant.value());
      }
      _carrier_phase_field.update(
          Phase::get_carrier_phase_data(mesh, instant.name(), params_phase),
          instant.value());
      _effective_drift.update(
          effective_drift_field(instant, _carrier_phase_field.field_new(),
                                params_transport, params_phase),
          instant.value());
    }

  private:
    VelocityField _velocity_field;
    PhaseField _carrier_phase_field;
    using DriftVelocityField = ptof::VectorField_Interpolation<
        Foam::volVectorField, decltype(_carrier_phase_field.locator()),
        InterpolationTypes::Cell, CheckOptions::Check>;
    using DriftField =
        Field_TimeInterpolation<DriftVelocityField, TimeInterpolationType>;
    DriftField _effective_drift;

    template <typename PhaseFieldValues>
    auto chemical_force(PhaseFieldValues const &carrier_phase,
                        Parameters const &params_phase) {
      if constexpr (std::is_same_v<ChemicalPotentialModel,
                                   ChemicalPotentialModels::AGG>) {
        return params_phase.chemical_potential_difference;
      }

      if constexpr (std::is_same_v<ChemicalPotentialModel,
                                   ChemicalPotentialModels::JM>) {
        return (1. - params_phase.phase_concentration_ratio) / 2. /
               (carrier_phase +
                (1. - carrier_phase) * params_phase.phase_concentration_ratio);
      }

      if constexpr (std::is_same_v<ChemicalPotentialModel,
                                   ChemicalPotentialModels::None>) {
        return 0.;
      }
    }

    /** \brief Compute chemical potential force field.
        \note Returns a scalar if constant. */
    auto chemical_force_field(Foam::scalar time,
                              Parameters const &params_phase) {
      // Avoid computing carrier phase field at given time if not needed
      if constexpr (std::is_same_v<ChemicalPotentialModel,
                                   ChemicalPotentialModels::AGG> ||
                    std::is_same_v<ChemicalPotentialModel,
                                   ChemicalPotentialModels::None>) {
        return chemical_force(meta::Empty{}, params_phase);
      } else {
        // Otherwise compute field
        return chemical_force(_carrier_phase_field.field(time), params_phase);
      }
    }

    template <typename PhaseFieldSpatialInterp>
    auto
    effective_drift_field(Foam::instant const &instant,
                          PhaseFieldSpatialInterp const &carrier_phase_field,
                          TransportParameters const &params_transport,
                          Parameters const &params_phase) {
      return static_cast<Foam::volVectorField>(
          Foam::dimensionedScalar("", Foam::dimViscosity,
                                  params_transport.diff_coeff) *
          chemical_force_field(instant.value(), params_phase) *
          grad_carrier_phase(carrier_phase_field.locator().mesh(),
                             instant.name(), carrier_phase_field.field(),
                             params_phase));
    }
  };
  template <typename VelocityField, typename PhaseField,
            typename TransportParameters>
  EffectiveVelocity(VelocityField &&, PhaseField &&,
                    TransportParameters const &, Parameters const &)
      -> EffectiveVelocity<VelocityField, PhaseField, TransportParameters>;

  template <typename Boundary, typename PhaseField, typename VelocityField>
  struct BoundaryPlusOptionalHardInterface {
    BoundaryPlusOptionalHardInterface(Boundary &boundary,
                                      PhaseField const &carrier_phase_field,
                                      VelocityField const &velocity_field,
                                      Parameters const &params_phase)
        : _boundary{boundary}, _carrier_phase_field{carrier_phase_field},
          _phase_tolerance{params_phase.phase_tolerance},
          _velocity_field{velocity_field},
          _grad_carrier_phase_field{make_grad_carrier_phase_field(
              _carrier_phase_field, params_phase)} {}

    BoundaryPlusOptionalHardInterface(Boundary &boundary,
                                      PhaseField const &carrier_phase_field,
                                      VelocityField &&velocity_field,
                                      Parameters const &params_phase) = delete;

    BoundaryPlusOptionalHardInterface(Boundary &boundary,
                                      PhaseField &&carrier_phase_field,
                                      VelocityField const &velocity_field,
                                      Parameters const &params_phase) = delete;

    BoundaryPlusOptionalHardInterface(Boundary &boundary,
                                      PhaseField &&carrier_phase_field,
                                      VelocityField &&velocity_field,
                                      Parameters const &params_phase) = delete;

    template <typename State>
    auto operator()(State &state, State const &state_old) {
      if constexpr (hard_reflection) {
        bool any_collision = false;
        bool collision = true;
        while (collision) {
          auto boundary_collision = _boundary(state, state_old);
          bool interface_collision = phase_hard_reflect(
              state, state_old, boundary_collision.intersection.point(),
              boundary_collision.time, _carrier_phase_field,
              _grad_carrier_phase_field, _phase_tolerance, _velocity_field);
          collision = boundary_collision || interface_collision;
          any_collision = any_collision || collision;
        }
        return any_collision;
      } else {
        return _boundary(state, state_old);
      }
    }

    void update(Foam::instant const &instant, Parameters const &params_phase) {
      if constexpr (hard_reflection) {
        _grad_carrier_phase_field.update(
            grad_carrier_phase(_carrier_phase_field.locator().mesh(),
                               _carrier_phase_field.time_new(),
                               _carrier_phase_field.field_new(), params_phase));
      }
    }

    auto static make_grad_carrier_phase_field(
        PhaseField const &carrier_phase_field, Parameters const &params_phase) {
      if constexpr (hard_reflection) {
        return Field_TimeInterpolation{
            std::make_unique<GradPhaseFieldSpatial>(
                grad_carrier_phase(carrier_phase_field.locator().mesh(),
                                   carrier_phase_field.time_new(),
                                   carrier_phase_field.field_new().field(),
                                   params_phase),
                carrier_phase_field.locator()),
            std::make_unique<GradPhaseFieldSpatial>(
                grad_carrier_phase(carrier_phase_field.locator().mesh(),
                                   carrier_phase_field.time_old(),
                                   carrier_phase_field.field_old().field(),
                                   params_phase),
                carrier_phase_field.locator()),
            carrier_phase_field.time_new(), carrier_phase_field.time_old(),
            TimeInterpolationType{}};
      } else {
        return meta::Empty{};
      }
    }

  private:
    Boundary &_boundary;
    PhaseField const &_carrier_phase_field;
    double _phase_tolerance;
    VelocityField const &_velocity_field;
    using GradPhaseFieldData = Foam::volVectorField;
    using GradPhaseFieldSpatial = VectorField_Interpolation<
        GradPhaseFieldData, decltype(_carrier_phase_field.locator()) const &,
        InterpolationTypes::Cell, CheckOptions::Check>;
    using GradPhaseField = decltype(make_grad_carrier_phase_field(
        std::declval<PhaseField>(), std::declval<Parameters>()));
    GradPhaseField _grad_carrier_phase_field;

  public:
    using SurfaceReaction = decltype(_boundary.surface_reaction);
    SurfaceReaction &surface_reaction{_boundary.surface_reaction};

  private:
    /** \note Does not explicitly handle multiple interface collisions. */
    template <typename State, typename Position>
    auto phase_hard_reflect(State &state, State const &state_old,
                            Position const &last_contact_point,
                            double time_last_collision,
                            PhaseField const &carrier_phase_field,
                            GradPhaseField const &grad_carrier_phase_field,
                            double phase_tolerance,
                            VelocityField const &velocity_field) const {
      // If inside carrier phase, do nothing
      if (carrier_phase_field(state) >= (1. - phase_tolerance)) {
        return false;
      }

      auto cell_last_collision =
          carrier_phase_field.locator()(last_contact_point, state_old.cell);
      auto position_last_collision = State::make_position(last_contact_point);
      // If the last collision was already outside carrier phase, use old
      // state
      if (carrier_phase_field(position_last_collision, cell_last_collision,
                              time_last_collision) < (1. - phase_tolerance)) {
        position_last_collision = state_old.position;
        time_last_collision = state_old.time;
      }

      // If no time elapsed, do nothing
      double time_step_particle = state.time - time_last_collision;
      if (time_step_particle == 0.) {
        return false;
      }

      // Find collision position and time taking fields at collision time
      double time_step_fields =
          carrier_phase_field.time_new() - carrier_phase_field.time_old();
      double tolerance = 1e-15;
      std::size_t max_iterations = 10;
      double fraction_guess = 0.5;
      auto displacement = op::minus(state.position, position_last_collision);
      double fraction_of_displacement_to_collision = roots::NewtonRaphson{}(
          [&carrier_phase_field, &position_last_collision, &displacement,
           cell_last_collision, time_last_collision, time_step_particle,
           phase_tolerance](double frac) {
            auto contact_point = op::plus(position_last_collision,
                                          op::times_scalar(frac, displacement));
            auto cell = carrier_phase_field.locator()(contact_point,
                                                      cell_last_collision);
            return carrier_phase_field(contact_point, cell,
                                       time_last_collision +
                                           frac * time_step_particle) -
                   (1. - phase_tolerance);
          },
          [&carrier_phase_field, &grad_carrier_phase_field,
           &position_last_collision, &displacement, cell_last_collision,
           time_last_collision, time_step_particle,
           time_step_fields](double frac) {
            auto contact_point = op::plus(position_last_collision,
                                          op::times_scalar(frac, displacement));
            auto cell = carrier_phase_field.locator()(contact_point,
                                                      cell_last_collision);
            return op::dot(displacement,
                           grad_carrier_phase_field(contact_point, cell,
                                                    time_last_collision)) +
                   (carrier_phase_field(contact_point,
                                        carrier_phase_field.time_old()) -
                    carrier_phase_field(contact_point,
                                        carrier_phase_field.time_old())) *
                       time_step_particle / time_step_fields;
          },
          fraction_guess, tolerance, max_iterations);
      if (fraction_of_displacement_to_collision < 0. ||
          fraction_of_displacement_to_collision > 1.) {
        fraction_of_displacement_to_collision = 0.5;
      }
      auto contact_point =
          op::times_scalar(fraction_of_displacement_to_collision, displacement);
      double collision_time =
          time_last_collision +
          fraction_of_displacement_to_collision * time_step_particle;
      auto collision_cell =
          carrier_phase_field.locator()(contact_point, cell_last_collision);
      auto interface_normal = op::normalize(grad_carrier_phase_field(
          contact_point, collision_cell, time_last_collision));

      // If there is a velocity field, make displacement relative to interface
      if constexpr (!std::is_same_v<VelocityField, meta::Empty>) {
        op::minus(displacement, velocity_field(contact_point, collision_cell,
                                               collision_time));
      }
      auto reflection_displacement = op::times_scalar(
          2. * op::dot(interface_normal, displacement) *
              (state.time - collision_time) / time_step_particle,
          interface_normal);
      op::minus_inplace(state.position, reflection_displacement);
      state.cell = carrier_phase_field.locate(state);

      return true;
    }
  };
  // This is added because some versions of clangd complain, not needed with
  // gcc. Try to remove later.
  template <typename Boundary, typename PhaseField, typename VelocityField>
  BoundaryPlusOptionalHardInterface(Boundary &, PhaseField const &,
                                    VelocityField const &, Parameters const &)
      -> BoundaryPlusOptionalHardInterface<Boundary, PhaseField, VelocityField>;
};
} // namespace ptof

#endif /* PTOF_PHASE_H */
