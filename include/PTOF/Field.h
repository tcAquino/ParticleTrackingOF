/**
 * @file   Field.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Wed Feb 16 00:00:00 2022
 *
 * @brief Field objects with interpolation and other utilities.
 */

#ifndef PTOF_FIELD_H
#define PTOF_FIELD_H

#include "General/Operation.h"
#include "PTOF/CheckOptions.h"
#include "PTOF/Useful.h"
#include <Vector2D.H>
#include <fieldTypes.H>
#include <interpolationCell.H>
#include <interpolationCellPoint.H>
#include <memory>
#include <point.H>
#include <stdexcept>
#include <type_traits>
#include <vector.H>

namespace ptof {
/** @brief Types to choose interpolation scheme. */
struct InterpolationTypes {
  /** @brief Use local cell value (do not interpolate). */
  struct Cell {};

  /** @brief Linear interpolation. */
  struct Linear {};

  /** @brief Use value at previous time. */
  struct OldTime {};
};

/**
 * @brief Interpolation of vector field cell data based on OpenFOAM routines.
 */
template <typename Field_t, typename Locator_t, typename InterpolationType,
          typename CheckOption>
class VectorField_Interpolation {
  static_assert(
      std::disjunction_v<std::is_same<CheckOption, CheckOptions::NoCheck>,
                         std::is_same<CheckOption, CheckOptions::Check>,
                         std::is_same<CheckOption, CheckOptions::Warn>>,
      "CheckOption must be CheckOptions::NoCheck, CheckOptions::Check, or "
      "CheckOptions::Warn");

  static_assert(
      std::disjunction_v<
          std::is_same<InterpolationType, InterpolationTypes::Cell>,
          std::is_same<InterpolationType, InterpolationTypes::Linear>>,
      "Interpolation type must be InterpolationTypes::Cell or "
      "InterpolationType::Linear");

public:
  using Point = Foam::point;                     /**< 3D point. */
  using Point2D = Foam::Vector2D<Foam::scalar>;  /**< 2D point. */
  using Vector = Foam::vector;                   /**< 3D vector. */
  using Vector2D = Foam::Vector2D<Foam::scalar>; /**< 2D vector. */
  using Scalar = Foam::scalar;                   /**< Scalar (also 1D point). */
  using Index = Foam::label;                     /**< Cell index. */
  using Field = Field_t;                         /**< Underlying field type. */
  using Locator = Locator_t;                     /**< Locator type. */

  /** Whether to check if requested positions are outside of mesh. */
  static constexpr bool check_if_outside =
      !std::is_same_v<CheckOption, CheckOptions::NoCheck>;

  /** Whether to warn when requested positions are outside of mesh. */
  static constexpr bool warn_if_outside =
      std::is_same_v<CheckOption, CheckOptions::Warn>;

  /**
   * @brief Constructor.
   *
   * @param field Vector field cell data to interpolate.
   *
   * @param locator Mesh locator.
   */
  VectorField_Interpolation(Field &&field, Locator &&locator)
      : _field{std::forward<Field>(field)}, _locator{std::forward<Locator>(
                                                locator)} {}

  /**
   * @brief Constructor.
   *
   * @param field Vector field cell data to interpolate.
   *
   * @param locator Mesh locator.
   */
  VectorField_Interpolation(Field &&field, Locator &&locator, InterpolationType,
                            CheckOption)
      : VectorField_Interpolation{std::forward<Field>(field),
                                  std::forward<Locator>(locator)} {}

  /** @brief Deleted copy constructor. */
  VectorField_Interpolation(VectorField_Interpolation const &) = delete;

  /** @brief Move constructor. */
  VectorField_Interpolation(VectorField_Interpolation &&) = default;

  /**
   * @brief Interpolate field.
   *
   * @param position 3D position.
   *
   * @param cell Mesh cell index position is in.
   *
   * @return interpolated field value.
   */
  auto operator()(Point const &position, Index cell) const {
    if constexpr (check_if_outside) {
      if (outside<warn_if_outside>(
              cell, position,
              "Assigning vector field value from nearest cell")) {
        cell = _locator.nearest_cell(position);
        return _field[cell];
      }
    }
    return _interpolator->interpolate(position, cell);
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 2D position.
   *
   * @param cell Mesh cell index position is in.
   *
   * @return interpolated field value.
   */
  auto operator()(Point2D const &position, Index cell) const {
    auto interp = (*this)(make_point(position), cell);
    return Vector2D{interp[0], interp[1]};
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 1D position.
   *
   * @param cell Mesh cell index position is in.
   *
   * @return interpolated field value.
   */
  auto operator()(Scalar position, Index cell) const {
    return (*this)(make_point(position), cell)[0];
  }

  /**
   * @brief Interpolate field.
   *
   * @param state Particle state to interpolate.
   */
  template <typename State> auto operator()(State const &state) const {
    return (*this)(state.position, state.cell);
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 3D position.
   *
   * @return interpolated field value.
   */
  auto operator()(Point const &position) const {
    return (*this)(position, _locator(position));
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 2D position.
   *
   * @return interpolated field value.
   */
  auto operator()(Point2D const &position) const {
    return (*this)(position, _locator(position));
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 1D position.
   *
   * @return interpolated field value.
   */
  auto operator()(Scalar const &position) const {
    return (*this)(position, _locator(position));
  }

  /** @brief Field at cell center. */
  auto operator[](Foam::label cell_id) const { return _field[cell_id]; }

  /** @brief Field at cell center. */
  auto operator()(Foam::label cell_id) const { return _field[cell_id]; }

  /**
   * @brief Locate a state in the mesh.
   *
   * @param state Particle state to locate.
   *
   * @return Mesh cell index.
   */
  template <typename State> auto locate(State const &state) const {
    return _locator(state);
  }

  /** @return Underlying field data. */
  auto const &field() const { return _field; }

  /** @return Locator object to find positions in mesh. */
  auto const &locator() const { return _locator; }

  /** @return Underlying boundary field. */
  auto const &boundaryField() const { return _field.boundaryField(); }

  /** @return Underlying boundary field at patch. */
  auto const &boundaryField(Foam::label patch_id) const {
    return _field.boundaryField()[patch_id];
  }

  /**
   * @return Underlying boundary field value at face of patch.
   *
   * @note Face id in patch is numbered relative to the patch, starting from 0.
   */
  auto const &boundaryField(Foam::label patch_id,
                            Foam::label face_in_patch) const {
    return _field.boundaryField()[patch_id][face_in_patch];
  }

  /** @brief Update interpolator. */
  void update_interpolator() {
    // Note: Not free, but OpenFOAM interpolators usually cash relevant data
    // like weights.
    _interpolator.reset(new Interpolator{_field});
  }

  /**
   * @brief Rescale underlying field.
   *
   * @param factor Scaling factor.
   */
  void rescale(double factor) {
    if (ptof::rescale(_field, factor)) {
      update_interpolator();
    }
  }

  /**
   * @brief Sum to underlying field.
   *
   * @param field Field to sum.
   */
  void sum(Field field) {
    _field += field;
    update_interpolator();
  }

  /**
   * @brief Change underlying field data.
   *
   * @param field Set field data to this field.
   */
  void set(Field field) {
    _field = field;
    update_interpolator();
  }

private:
  Field _field;     /**< Vector field cell data to interpolate. */
  Locator _locator; /**< Locator to find positions in mesh. */

  using Interpolator = std::conditional_t<
      std::is_same_v<InterpolationType, InterpolationTypes::Linear>,
      Foam::interpolationCellPoint<Vector>, Foam::interpolationCell<Vector>>;
  std::unique_ptr<Interpolator> _interpolator{
      make_interpolator()}; /**< Interpolation object. */

  std::unique_ptr<Interpolator> make_interpolator() {
    return std::make_unique<Interpolator>(_field);
  }
};
template <typename Field, typename Locator, typename InterpolationType,
          typename CheckOption>
VectorField_Interpolation(Field &&, Locator &&, InterpolationType, CheckOption)
    -> VectorField_Interpolation<Field, Locator, InterpolationType,
                                 CheckOption>;

/**
 * @brief Interpolation of scalar field cell data based on OpenFOAM routines.
 *
 * @note See ptof::InterpolationTypes class for interpolation options and
 *       ptof::CheckOptions class for bounds checking options.
 */
template <typename Field_t, typename Locator_t, typename InterpolationType,
          typename CheckOption>
class ScalarField_Interpolation {
  static_assert(
      std::disjunction_v<std::is_same<CheckOption, CheckOptions::NoCheck>,
                         std::is_same<CheckOption, CheckOptions::Check>,
                         std::is_same<CheckOption, CheckOptions::Warn>>,
      "CheckOption must be CheckOptions::NoCheck, CheckOptions::Check, or "
      "CheckOptions::Warn");

  static_assert(
      std::disjunction_v<
          std::is_same<InterpolationType, InterpolationTypes::Cell>,
          std::is_same<InterpolationType, InterpolationTypes::Linear>>,
      "Interpolation type must be InterpolationTypes::Cell or "
      "InterpolationType::Linear");

public:
  using Point = Foam::point;                    /**< 3D point. */
  using Point2D = Foam::Vector2D<Foam::scalar>; /**< 2D point. */
  using Scalar = Foam::scalar;                  /**< Scalar (also 1D point). */
  using Index = Foam::label;                    /**< Cell index. */
  using Field = Field_t;                        /**< Underlying field type. */
  using Locator = Locator_t;                    /**< Locator type. */

  /** Whether to check if requested positions are outside of mesh. */
  static constexpr bool check_if_outside =
      !std::is_same_v<CheckOption, CheckOptions::NoCheck>;

  /** Whether to warn when requested positions are outside of mesh. */
  static constexpr bool warn_if_outside =
      std::is_same_v<CheckOption, CheckOptions::Warn>;

  /**
   * @brief Constructor.
   *
   * @param field Scalar field cell data to interpolate.
   *
   * @param locator Mesh locator.
   */
  ScalarField_Interpolation(Field &&field, Locator &&locator)
      : _field{std::forward<Field>(field)}, _locator{std::forward<Locator>(
                                                locator)} {}

  /**
   * @brief Constructor.
   *
   * @param field Scalar field cell data to interpolate.
   *
   * @param locator Mesh locator.
   */
  ScalarField_Interpolation(Field &&field, Locator &&locator, InterpolationType,
                            CheckOption)
      : ScalarField_Interpolation{std::forward<Field>(field),
                                  std::forward<Locator>(locator)} {}

  /** @brief Deleted copy constructor. */
  ScalarField_Interpolation(ScalarField_Interpolation const &) = delete;

  /** @brief Move constructor. */
  ScalarField_Interpolation(ScalarField_Interpolation &&field) = default;

  /**
   * @brief Interpolate field.
   *
   * @param position 3D position.
   *
   * @param cell Mesh cell index position is in.
   *
   * @return interpolated field value.
   */
  auto operator()(Point const &position, Index cell) const {
    if constexpr (check_if_outside) {
      if (outside<warn_if_outside>(
              cell, position,
              "Assigning scalar field value from nearest cell")) {
        cell = _locator.nearest_cell(position);
        return _field[cell];
      }
    }
    return _interpolator->interpolate(position, cell);
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 2D position.
   *
   * @param cell Mesh cell index position is in.
   *
   * @return interpolated field value.
   */
  auto operator()(Point2D const &position, Index cell) const {
    return (*this)(make_point(position), cell);
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 1D position.
   *
   * @param cell Mesh cell index position is in.
   *
   * @return interpolated field value.
   */
  auto operator()(Scalar position, Index cell) const {
    return (*this)(make_point(position), cell);
  }

  /**
   * @brief Interpolate field.
   *
   * @param state Particle state to interpolate.
   */
  template <typename State> auto operator()(State const &state) const {
    return (*this)(state.position, state.cell);
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 3D position.
   *
   * @return interpolated field value.
   */
  auto operator()(Point const &position) const {
    return (*this)(position, _locator(position));
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 2D position.
   *
   * @return interpolated field value.
   */
  auto operator()(Point2D const &position) const {
    return (*this)(position, _locator(position));
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 1D position.
   *
   * @return interpolated field value.
   */
  auto operator()(Scalar const &position) const {
    return (*this)(position, _locator(position));
  }

  /** @brief Field at cell center. */
  auto operator[](Foam::label cell_id) const { return _field[cell_id]; }

  /** @brief Field at cell center. */
  auto operator()(Foam::label cell_id) const { return _field[cell_id]; }

  /**
   * @brief Locate a state in the mesh.
   *
   * @param state Particle state to locate.
   *
   * @return Mesh cell index.
   */
  template <typename State> auto locate(State const &state) const {
    return _locator(state);
  }

  /** @return Underlying field data. */
  auto const &field() const { return _field; }

  /** @return Locator object to find positions in mesh. */
  auto const &locator() const { return _locator; }

  /** @return Underlying boundary field. */
  auto const &boundaryField() const { return _field.boundaryField(); }

  /** @return Underlying boundary field at patch. */
  auto const &boundaryField(Foam::label patch_id) const {
    return _field.boundaryField()[patch_id];
  }

  /**
   * @return Underlying boundary field value at face of patch.
   *
   * @note Face id in patch is numbered relative to the patch, starting from 0.
   */
  auto const &boundaryField(Foam::label patch_id,
                            Foam::label face_in_patch) const {
    return _field.boundaryField()[patch_id][face_in_patch];
  }

  /** @brief Update interpolator. */
  void update_interpolator() {
    // Note: Not free, but OpenFOAM interpolators usually cash relevant data
    // like weights.
    _interpolator.reset(new Interpolator{_field});
  }

  /**
   * @brief Rescale underlying field.
   *
   * @param factor Scaling factor.
   */
  void rescale(double factor) {
    if (ptof::rescale(_field, factor)) {
      update_interpolator();
    }
  }

  /**
   * @brief Sum to underlying field.
   *
   * @param field field to sum
   */
  void sum(Field field) {
    _field += field;
    update_interpolator();
  }

  /**
   * @brief Change underlying field data.
   *
   * @param field Set field data to this field.
   */
  void set(Field field) {
    _field = field;
    update_interpolator();
  }

private:
  Field _field;     /**< Scalar field cell data to interpolate. */
  Locator _locator; /**< Locator to find positions in mesh. */

  using Interpolator = std::conditional_t<
      std::is_same_v<InterpolationType, InterpolationTypes::Linear>,
      Foam::interpolationCellPoint<Scalar>, Foam::interpolationCell<Scalar>>;
  std::unique_ptr<Interpolator> _interpolator{
      make_interpolator()}; /**< Interpolation object. */

  std::unique_ptr<Interpolator> make_interpolator() {
    return std::make_unique<Interpolator>(_field);
  }
};
template <typename Field, typename Locator, typename InterpolationType,
          typename CheckOption>
ScalarField_Interpolation(Field &&, Locator &&, InterpolationType, CheckOption)
    -> ScalarField_Interpolation<Field, Locator, InterpolationType,
                                 CheckOption>;

/**
 * @brief Interpolation of vector field cell data based on OpenFOAM routines.
 */
template <typename Field_t, typename TimeInterpolationType>
class Field_TimeInterpolation {
public:
  using Point = Foam::point;                     /**< 3D point. */
  using Point2D = Foam::Vector2D<Foam::scalar>;  /**< 2D point. */
  using Vector = Foam::vector;                   /**< 3D vector. */
  using Vector2D = Foam::Vector2D<Foam::scalar>; /**< 2D vector. */
  using Scalar = Foam::scalar;                   /**< Scalar (also 1D point). */
  using Index = Foam::label;                     /**< Cell index. */
  using Field = Field_t;                         /**< Underlying field type. */
  using Time = Foam::scalar;                     /**< Time type. */
  using FieldPtr =
      std::unique_ptr<Field>; /**< Pointer to underlying field type. */

  static_assert(
      std::disjunction_v<
          std::is_same<TimeInterpolationType, InterpolationTypes::Linear>,
          std::is_same<TimeInterpolationType, InterpolationTypes::OldTime>>,
      "Interpolation type must be InterpolationTypes::Linear or "
      "InterpolationTypes::OldTime");

  /**
   * @brief Constructor.
   *
   * @param field_new Field at new time.
   *
   * @param field_old Field at old time.
   *
   * @param time_new Time of new field.
   *
   * @param time_old Time of old field.
   */
  Field_TimeInterpolation(FieldPtr field_new, FieldPtr field_old,
                          double time_new, double time_old)
      : _field_new{std::move(field_new)}, _field_old{std::move(field_old)},
        _time_new{time_new}, _time_old{time_old} {}

  /**
   * @brief Constructor.
   *
   * @param field_new Field at new time.
   *
   * @param field_old Field at old time.
   *
   * @param time_new Time of new field.
   *
   * @param time_old Time of old field.
   */
  Field_TimeInterpolation(FieldPtr field_new, FieldPtr field_old,
                          double time_new, double time_old,
                          TimeInterpolationType)
      : _field_new{std::move(field_new)}, _field_old{std::move(field_old)},
        _time_new{time_new}, _time_old{time_old} {}

  /**
   * @brief Interpolate field.
   *
   * @param position Position.
   *
   * @param cell Mesh cell index position is in.
   *
   * @param time Interpolation time.
   *
   * @return interpolated field value.
   */
  template <typename Position>
  auto operator()(Position const &position, Index cell, Time time) const {
    return interpolate((*_field_new)(position, cell),
                       (*_field_old)(position, cell), time);
  }

  /**
   * @brief Interpolate field.
   *
   * @param state Particle state to interpolate.
   */
  template <typename State> auto operator()(State const &state) const {
    return (*this)(state.position, state.cell, state.time);
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 3D position.
   *
   * @param time Interpolation time.
   *
   * @return interpolated field value.
   */
  auto operator()(Point const &position, Time time) const {
    return (*this)(position, locate(position), time);
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 2D position.
   *
   * @param time Interpolation time.
   *
   * @return interpolated field value.
   */
  auto operator()(Point2D const &position, Time time) const {
    return (*this)(position, locate(position), time);
  }

  /**
   * @brief Interpolate field.
   *
   * @param position 1D position.
   *
   * @param time Interpolation time.
   *
   * @return interpolated field value.
   */
  auto operator()(Scalar position, Time time) const {
    return (*this)(position, locate(position), time);
  }

  /** @brief Field at cell center. */
  auto operator()(Foam::label cell_id, Time time) const {
    return interpolate((*_field_new)(cell_id), (*_field_old)(cell_id), time);
  }

  /**
   * @brief Locate a state in the mesh.
   *
   * @param state Particle state to locate.
   *
   * @return Mesh cell index.
   */
  template <typename State> auto locate(State const &state) const {
    return _field_old->locate(state);
  }

  /** @return Locator object to find positions in mesh. */
  auto const &locator() const { return _field_old->locator(); }

  /** @return Full field of underlying field type. */
  auto field(Time time) const {
    return interpolate(_field_new->field(), _field_old->field(), time);
  }

  /** @return Underlying boundary field. */
  auto boundaryField(Time time) const {
    return interpolate(_field_new->boundaryField(), _field_old->boundaryField(),
                       time);
  }

  /** @return Underlying boundary field at patch. */
  auto boundaryField(Foam::label patch_id, Time time) const {
    return interpolate(_field_new->boundaryField(patch_id),
                       _field_old->boundaryField(patch_id), time);
  }

  /**
   * @return Underlying boundary field value at face of patch.
   *
   * @note Face id in patch is numbered relative to the patch, starting from 0.
   */
  auto boundaryField(Foam::label patch_id, Foam::label face_in_patch,
                     Time time) const {
    return interpolate(_field_new->boundaryField(patch_id, face_in_patch),
                       _field_old->boundaryField(patch_id, face_in_patch),
                       time);
  }

  /** @return New field. */
  auto const &field_new() const { return *_field_new; }

  /** @return Old field. */
  auto const &field_old() const { return *_field_old; }

  /** @return Time of new field. */
  auto time_new() const { return _time_new; }

  /** @return Time of old field. */
  auto time_old() const { return _time_old; }

  /**
   * @brief Set new field data.
   */
  template <typename Field> void set_old(Field field) {
    _field_old->set(field);
  }

  /**
   * @brief Set old field data.
   */
  template <typename Field> void set_new(Field field) {
    _field_new->set(field);
  }

  /**
   * @brief Update underlying field by updating the underlying data.
   *
   * @param field Set new underlying field data to this field.
   *
   * @param time_new Time of new field.
   */
  template <typename Field> void update(Field field, Time time_new) {
    _field_new.swap(_field_old);
    set_new(field);
    _time_old = _time_new;
    _time_new = time_new;
  }

  /**
   * @brief Update underlying field.
   *
   * @param field Set new field to this field.
   *
   * @param time_new Time of new field.
   */
  void update(FieldPtr field, Time time_new) {
    _field_old = std::move(_field_new);
    _field_new = std::move(field);
    _time_old = _time_new;
    _time_new = time_new;
  }

  /**
   * @brief Rescale underlying field.
   *
   * @param factor Scaling factor.
   */
  void rescale(double factor) {
    rescale_new(factor);
    rescale_old(factor);
  }

  /**
   * @brief Rescale underlying field.
   *
   * @param factor Scaling factor.
   */
  void rescale_new(double factor) { _field_new->rescale(factor); }

  /**
   * @brief Rescale underlying field.
   *
   * @param factor Scaling factor.
   */
  void rescale_old(double factor) { _field_old->rescale(factor); }

private:
  FieldPtr _field_new; /**< Field at new time. */
  FieldPtr _field_old; /**< Field at old time. */
  Time _time_new;      /**< Time of new field. */
  Time _time_old;      /**< Time of old field. */

  template <typename FieldVal>
  auto interpolate(FieldVal const &field_new, FieldVal const &field_old,
                   Time time) const {
    if constexpr (std::is_same_v<TimeInterpolationType,
                                 InterpolationTypes::Linear>) {
      if (_time_new <= _time_old || time <= _time_old) {
        return field_old;
      }
      if (time >= _time_new) {
        return field_new;
      }
      return static_cast<FieldVal>(field_old + (time - _time_old) /
                                                   (_time_new - _time_old) *
                                                   (field_new - field_old));
    }
    if constexpr (std::is_same_v<TimeInterpolationType,
                                 InterpolationTypes::OldTime>) {
      return field_old;
    }
  }
};
} // namespace ptof

#endif /* PTOF_FIELD_H */
