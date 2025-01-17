/**
   \file PTOF/Field.h
   \author Tomás Aquino
   \date 16/02/2022
   \brief Field objects with interpolation and other utilities.
*/

#ifndef PTOF_FIELD_H
#define PTOF_FIELD_H

#include "General/Meta.h"
#include "PTOF/Useful.h"
#include <Vector2D.H>
#include <fieldTypes.H>
#include <interpolationCellPoint.H>
#include <point.H>
#include <type_traits>
#include <vector.H>

namespace ptof {
/**
   \class VectorField_LinearInterpolation_OF PTOF/Field.h "PTOF/Field.h"
   \brief Linear interpolation of vector field cell data based on OpenFOAM
   routines.
   \note See CheckOptions class for bounds checking options.
*/
template <typename Field, typename Locator, typename CheckOption,
          typename Uninterpolated = meta::Empty>
class VectorField_LinearInterpolation_OF {
public:
  using Point = Foam::point;                     /**> 3D point. */
  using Point2D = Foam::Vector2D<Foam::scalar>;  /**> 2D point. */
  using Vector = Foam::vector;                   /**> 3D vector. */
  using Vector2D = Foam::Vector2D<Foam::scalar>; /**> 2D vector. */
  using Scalar = Foam::scalar;                   /**> Scalar (also 1D point). */
  using Index = Foam::label;                     /**> Cell index. */

  /** Whether to check if requested positions are outside of mesh. */
  static constexpr bool check_if_outside =
      !std::is_same_v<CheckOption, CheckOptions::NoCheck>;

  /** Whether to warn when requested positions are outside of mesh. */
  static constexpr bool warn_if_outside =
      std::is_same_v<CheckOption, CheckOptions::Warn>;

  /**
     \brief Constructor.
     \param field Vector field cell data to interpolate.
     \param locator Mesh locator.
  */
  VectorField_LinearInterpolation_OF(Field &&field, Locator &&locator)
      : _field{std::forward<Field>(field)}, _locator{std::forward<Locator>(
                                                locator)} {}

  /**
     \brief Constructor.
     \param field Vector field cell data to interpolate.
     \param locator Mesh locator.
     \param uninterpolated Uninterpolated field to add to interpolated field.
  */
  VectorField_LinearInterpolation_OF(Field &&field, Locator &&locator,
                                     Uninterpolated &&uninterpolated)
      : _field{std::forward<Field>(field)}, _locator{std::forward<Locator>(
                                                locator)},
        _uninterpolated{std::forward<Uninterpolated>(uninterpolated)} {}

  /**
     \brief Constructor.
     \param field Vector field cell data to interpolate.
     \param locator Mesh locator.
  */
  VectorField_LinearInterpolation_OF(Field &&field, Locator &&locator,
                                     meta::Selector_t<CheckOption>)
      : VectorField_LinearInterpolation_OF{std::forward<Field>(field),
                                           std::forward<Locator>(locator)} {}

  /**
     \brief Constructor.
     \param field Vector field cell data to interpolate.
     \param locator Mesh locator.
     \param uninterpolated Uninterpolated field to add to interpolated field.
  */
  VectorField_LinearInterpolation_OF(Field &&field, Locator &&locator,
                                     Uninterpolated &&uninterpolated,
                                     meta::Selector_t<CheckOption>)
      : VectorField_LinearInterpolation_OF{
            std::forward<Field>(field), std::forward<Locator>(locator),
            std::forward<Uninterpolated>(uninterpolated)} {}

  /** \brief Deleted copy constructor. */
  VectorField_LinearInterpolation_OF(
      VectorField_LinearInterpolation_OF const &) = delete;

  /** \brief Move constructor. */
  VectorField_LinearInterpolation_OF(VectorField_LinearInterpolation_OF &&field)
      : _field{std::move(field._field)}, _locator{std::move(field._locator)},
        _interpolant{this->_field}, _uninterpolated{
                                        std::move(field._uninterpolated)} {}

  /**
     \brief Interpolate field.
     \param position 3D position.
     \param cell Mesh cell index position is in.
     \return interpolated field value.
  */
  auto operator()(Point const &position, Index cell) const {
    if constexpr (check_if_outside)
      if (outside<warn_if_outside>(
              cell, position,
              "Assigning vector field value from nearest cell")) {
        cell = _locator.nearest_cell(position);
        if constexpr (!std::is_same_v<Uninterpolated, meta::Empty>)
          return _field[cell] + _uninterpolated[cell];
        else
          return _field[cell];
      }

    if constexpr (!std::is_same_v<Uninterpolated, meta::Empty>)
      return _interpolant.interpolate(position, cell) + _uninterpolated[cell];

    return _interpolant.interpolate(position, cell);
  }

  /**
     \brief Interpolate field.
     \param position 2D position.
     \param cell Mesh cell index position is in.
     \return interpolated field value.
  */
  auto operator()(Point2D const &position, Index cell) const {
    auto interp = (*this)(make_point(position), cell);

    return Vector2D{interp[0], interp[1]};
  }

  /**
     \brief Interpolate field.
     \param position 1D position.
     \param cell Mesh cell index position is in.
     \return interpolated field value.
  */
  auto operator()(Scalar position, Index cell) const {
    return (*this)(make_point(position), cell)[0];
  }

  /**
     \brief Interpolate field.
     \param state Particle state to interpolate.
  */
  template <typename State> auto operator()(State const &state) const {
    return (*this)(state.position, state.cell);
  }

  /**
     \brief Interpolate field.
     \param position 3D position.
     \return interpolated field value.
  */
  auto operator()(Point const &position) const {
    return (*this)(position, _locator(position));
  }

  /**
     \brief Interpolate field.
     \param position 2D position.
     \return interpolated field value.
  */
  auto operator()(Point2D const &position) const {
    return (*this)(position, _locator(position));
  }

  /**
     \brief Interpolate field.
     \param position 1D position.
     \return interpolated field value.
  */
  auto operator()(Scalar const &position) const {
    return (*this)(position, _locator(position));
  }

  /**
     \brief Locate a state in the mesh.
     \param state Particle state to locate.
     \return Mesh cell index.
  */
  template <typename State> auto locate(State const &state) const {
    return _locator(state);
  }

  /** \return Underlying field data. */
  auto const &field() const { return _field; }

  /** \return Underlying uninterpolated data. */
  auto const &uninterpolated() const { return _uninterpolated; }

  /** \return Locator object to find positions in mesh. */
  auto const &locator() const { return _locator; }

  /** \return Underlying boundary field. */
  auto const &boundaryField() const { return _field.boundaryField(); }

  /** \return Underlying boundary field. */
  auto &boundaryFieldRef() const { return _field.boundaryFieldRef(); }

  /**
     \brief Rescale underlying field.
     \param factor Scaling factor.
  */
  void rescale(double factor) { _field *= factor; }

  /**
     \brief Rescale underlying uninterpolated field.
     \param factor Scaling factor.
  */
  void rescale_uninterpolated(double factor) { _uninterpolated *= factor; }

  /**
     \brief Sum to underlying field.
     \param field Field to sum.
  */
  void sum(Field const &field) { _field += field; }

  /**
     \brief Sum to underlying uninterpolated field.
     \param uninterpolated Field to sum.
  */
  void sum_uninterpolated(Uninterpolated const &uninterpolated) {
    _uninterpolated += uninterpolated;
  }

  /**
     \brief Change underlying field data.
     \param field Set field data to this field.
  */
  void set(Field const &field) { _field = field; }

  /**
     \brief Change underlying uninterpolated field data.
     \param uninterpolated Set uninterpolated field data to this field.
  */
  void set_uninterpolated(Uninterpolated const &uninterpolated) {
    _uninterpolated = uninterpolated;
  }

private:
  Field _field;     /**< Vector field cell data to interpolate. */
  Locator _locator; /**< Locator to find positions in mesh. */
  Foam::interpolationCellPoint<Vector> _interpolant{
      _field}; /**< Interpolation object. */
  Uninterpolated
      _uninterpolated; /** Uninterpolated field to add to interpolated field. */
};
template <typename Field, typename Locator, typename CheckOption>
VectorField_LinearInterpolation_OF(Field &&, Locator &&,
                                   meta::Selector_t<CheckOption>)
    -> VectorField_LinearInterpolation_OF<Field, Locator, CheckOption,
                                          meta::Empty>;
template <typename Field, typename Locator, typename Uninterpolated,
          typename CheckOption>
VectorField_LinearInterpolation_OF(Field &&, Locator &&, Uninterpolated &&,
                                   meta::Selector_t<CheckOption>)
    -> VectorField_LinearInterpolation_OF<Field, Locator, CheckOption,
                                          Uninterpolated>;

/**
   \class ScalarField_LinearInterpolation_OF PTOF/Field.h "PTOF/Field.h"
   \brief Linear interpolation of scalar field cell data based on OpenFOAM
   routines.
   \note See CheckOptions class for bounds checking options.
*/
template <typename Field, typename Locator, typename CheckOption,
          typename Uninterpolated = meta::Empty>
class ScalarField_LinearInterpolation_OF {
public:
  using Point = Foam::point;                    /**> 3D point. */
  using Point2D = Foam::Vector2D<Foam::scalar>; /**> 2D point. */
  using Scalar = Foam::scalar;                  /**> Scalar (also 1D point). */
  using Index = Foam::label;                    /**> Cell index. */

  /** Whether to check if requested positions are outside of mesh. */
  static constexpr bool check_if_outside =
      !std::is_same_v<CheckOption, CheckOptions::NoCheck>;

  /** Whether to warn when requested positions are outside of mesh. */
  static constexpr bool warn_if_outside =
      std::is_same_v<CheckOption, CheckOptions::NoCheck>;

  /**
     \brief Constructor.
     \param field Scalar field cell data to interpolate.
     \param locator Mesh locator.
   */
  ScalarField_LinearInterpolation_OF(Field &&field, Locator &&locator)
      : _field{std::forward<Field>(field)}, _locator{std::forward<Locator>(
                                                locator)} {}

  /**
     \brief Constructor.
     \param field Scalar field cell data to interpolate.
     \param locator Mesh locator.
     \param uninterpolated Uninterpolated field to add to interpolated field.
  */
  ScalarField_LinearInterpolation_OF(Field &&field, Locator &&locator,
                                     Uninterpolated &&uninterpolated)
      : _field{std::forward<Field>(field)}, _locator{std::forward<Locator>(
                                                locator)},
        _uninterpolated{std::forward<Uninterpolated>(uninterpolated)} {}

  /**
     \brief Constructor.
     \param field Scalar field cell data to interpolate.
     \param locator Mesh locator.
  */
  ScalarField_LinearInterpolation_OF(Field &&field, Locator &&locator,
                                     meta::Selector_t<CheckOption>)
      : ScalarField_LinearInterpolation_OF{std::forward<Field>(field),
                                           std::forward<Locator>(locator)} {}

  /**
     \brief Constructor.
     \param field Scalar field cell data to interpolate.
     \param locator Mesh locator.
     \param uninterpolated Uninterpolated field to add to interpolated field.
  */
  ScalarField_LinearInterpolation_OF(Field &&field, Locator &&locator,
                                     Uninterpolated &&uninterpolated,
                                     meta::Selector_t<CheckOption>)
      : ScalarField_LinearInterpolation_OF{
            std::forward<Field>(field), std::forward<Locator>(locator),
            std::forward<Uninterpolated>(uninterpolated)} {}

  /** \brief Deleted copy constructor. */
  ScalarField_LinearInterpolation_OF(
      ScalarField_LinearInterpolation_OF const &) = delete;

  /** \brief Move constructor. */
  ScalarField_LinearInterpolation_OF(ScalarField_LinearInterpolation_OF &&field)
      : _field{std::move(field._field)}, _locator{std::move(field._locator)},
        _interpolant{this->_field}, _uninterpolated{
                                        std::move(field._uninterpolated)} {}

  /**
     \brief Interpolate field.
     \param position 3D position.
     \param cell Mesh cell index position is in.
     \return interpolated field value.
  */
  auto operator()(Point const &position, Index cell) const {
    if constexpr (check_if_outside)
      if (outside<warn_if_outside>(
              cell, position,
              "Assigning scalar field value from nearest cell")) {
        cell = _locator.nearest_cell(position);
        if constexpr (!std::is_same_v<Uninterpolated, meta::Empty>)
          return _field[cell] + _uninterpolated[cell];
        else
          return _field[cell];
      }

    if constexpr (!std::is_same_v<Uninterpolated, meta::Empty>)
      return _interpolant.interpolate(position, cell) + _uninterpolated[cell];

    return _interpolant.interpolate(position, cell);
  }

  /**
     \brief Interpolate field.
     \param position 2D position.
     \param cell Mesh cell index position is in.
     \return interpolated field value.
  */
  auto operator()(Point2D const &position, Index cell) const {
    return (*this)(make_point(position), cell);
  }

  /**
     \brief Interpolate field.
     \param position 1D position.
     \param cell Mesh cell index position is in.
     \return interpolated field value.
  */
  auto operator()(Scalar position, Index cell) const {
    return (*this)(make_point(position), cell);
  }

  /**
     \brief Interpolate field.
     \param state Particle state to interpolate.
  */
  template <typename State> auto operator()(State const &state) const {
    return (*this)(state.position, state.cell);
  }

  /**
     \brief Interpolate field.
     \param position 3D position.
     \return interpolated field value.
  */
  auto operator()(Point const &position) const {
    return (*this)(position, _locator(position));
  }

  /**
     \brief Interpolate field.
     \param position 2D position.
     \return interpolated field value.
  */
  auto operator()(Point2D const &position) const {
    return (*this)(position, _locator(position));
  }

  /**
     \brief Interpolate field.
     \param position 1D position.
     \return interpolated field value.
  */
  auto operator()(Scalar const &position) const {
    return (*this)(position, _locator(position));
  }

  /**
     \brief Locate a state in the mesh.
     \param state Particle state to locate.
     \return Mesh cell index.
  */
  template <typename State> auto locate(State const &state) const {
    return _locator(state);
  }

  /** \return Underlying field data. */
  auto const &field() const { return _field; }

  /** \return Underlying uninterpolated data. */
  auto const &uninterpolated() const { return _uninterpolated; }

  /** \return Locator object to find positions in mesh. */
  auto const &locator() const { return _locator; }

  /** \return Underlying boundary field. */
  auto const &boundaryField() const { return _field.boundaryField(); }

  /** \return Underlying boundary field. */
  auto &boundaryFieldRef() const { return _field.boundaryFieldRef(); }

  /**
     \brief Rescale underlying field.
     \param factor Scaling factor.
  */
  void rescale(double factor) { _field *= factor; }

  /**
     \brief Rescale underlying uninterpolated field.
     \param factor Scaling factor.
  */
  void rescale_uninterpolated(double factor) { _uninterpolated *= factor; }

  /**
     \brief Sum to underlying field.
     \param field field to sum
  */
  void sum(Field const &field) { _field += field; }

  /**
     \brief Sum to underlying uninterpolated field.
     \param uninterpolated field to sum
  */
  void sum_uninterpolated(Uninterpolated const &uninterpolated) {
    _uninterpolated += uninterpolated;
  }

  /**
     \brief Change underlying field data.
     \param field Set field data to this field.
  */
  void set(Field const &field) { _field = field; }

  /**
     \brief Change underlying uninterpolated field data.
     \param uninterpolated Set uninterpolated field data to this field.
  */
  void set_uninterpolated(Uninterpolated const &uninterpolated) {
    _uninterpolated = uninterpolated;
  }

private:
  Field _field;     /**< Scalar field cell data to interpolate.         */
  Locator _locator; /**< Locator to find positions in mesh. */
  Foam::interpolationCellPoint<Scalar> _interpolant{
      _field}; /**< Interpolation object. */
  Uninterpolated
      _uninterpolated; /** Uninterpolated field to add to interpolated field. */
};
template <typename Field, typename Locator, typename CheckOption>
ScalarField_LinearInterpolation_OF(Field &&, Locator &&,
                                   meta::Selector_t<CheckOption>)
    -> ScalarField_LinearInterpolation_OF<Field, Locator, CheckOption,
                                          meta::Empty>;
template <typename Field, typename Locator, typename Uninterpolated,
          typename CheckOption>
ScalarField_LinearInterpolation_OF(Field &&, Locator &&, Uninterpolated &&,
                                   meta::Selector_t<CheckOption>)
    -> ScalarField_LinearInterpolation_OF<Field, Locator, CheckOption,
                                          Uninterpolated>;

/**
   \brief Compute magnitude of average of volumetric field (set of values
   associated with mesh cells). */
template <typename Field, typename Mesh>
auto magnitude_of_average(Field &field, Mesh const &mesh) {
  Foam::scalar mesh_volume = Foam::sum(mesh.cellVolumes());
  auto average_weighted_data = Foam::sum(field * mesh.cellVolumes());
  return Foam::mag(average_weighted_data) / mesh_volume;
}

/**
   \brief Rescale a volumetric field (set of values associated with mesh cells)
   to a given average value. */
template <typename Field, typename Mesh>
void rescale_to_average(Field &field, Mesh const &mesh, double average) {
  field *= average / magnitude_of_average(field, mesh);
}
} // namespace ptof

#endif /* PTOF_FIELD_H */
