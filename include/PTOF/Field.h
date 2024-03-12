/**
 \file PTOF/Field.h
 \author Tomás Aquino
 \date 16/02/2022
*/

#ifndef PTOF_FIELD_H
#define PTOF_FIELD_H

#include <type_traits>
#include <fieldTypes.H>
#include <interpolationCellPoint.H>
#include "General/Useful.h"
#include "PTOF/Useful.h"

namespace ptof
{

  /** \class VectorField_LinearInterpolation_OF PTOF/Field.h "PTOF/Field.h"
    *  \brief Linear interpolation of vector field cell data based on OpenFOAM routines. */
  template
  <typename Field, typename Locator,
  bool check_if_outside, bool warn_if_outside,
  typename Uninterpolated = useful::Empty>
  class VectorField_LinearInterpolation_OF
  {
  public:
    using Point = Foam::point;                      /**> 3D point. */
    using Point2D = Foam::Vector2D<Foam::scalar>;   /**> 2D point. */
    using Vector = Foam::vector;                    /**> 3D vector. */
    using Vector2D = Foam::Vector2D<Foam::scalar>;  /**> 2D vector. */
    using Scalar = Foam::scalar;                    /**> Scalar (also 1D point). */
    using Index = Foam::label;                      /**> Cell index.*/
    
    /** Constructor.
    \brief Without checking if position is in bounds.
    \param field Vector field cell data to interpolate.
    \param locator Mesh locator.
    */
    VectorField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, CheckOptions::NoCheck)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    {
      static_assert(check_if_outside == false && warn_if_outside == false,
                    "Bad class template arguments for no bounds checking.");
    }
    
    /** Constructor.
     \brief With checking if position is in bounds but no warning if not.
     \param field Vector field cell data to interpolate.
     \param locator Mesh locator.
    */
    VectorField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, CheckOptions::Check)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    {
      static_assert(check_if_outside == true && warn_if_outside == false,
                    "Bad class template arguments for bounds checking.");
    }
    
    /** Constructor.
     \brief With warning if particle is in bounds.
     \param field Scalar field cell data to interpolate.
     \param locator Mesh locator.
    */
    VectorField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, CheckOptions::Warn)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    {
      static_assert(check_if_outside == true && warn_if_outside == true,
                    "Bad class template arguments for bounds warning.");
    }
    
    /** Constructor.
    \brief Without checking if position is in bounds.
    \param field Vector field cell data to interpolate.
    \param locator Mesh locator.
    \param uninterpolated Uninterpolated field to add to interpolated field.
    */
    VectorField_LinearInterpolation_OF
    (Field&& field, Locator&& locator,
     Uninterpolated&& uninterpolated,
     CheckOptions::NoCheck)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    , _uninterpolated{ uninterpolated }
    {
      static_assert(check_if_outside == false && warn_if_outside == false,
                    "Bad class template arguments for no bounds checking.");
    }
    
    /** Constructor.
     \brief With checking if position is in bounds but no warning if not.
     \param field Vector field cell data to interpolate.
     \param locator Mesh locator.
     \param uninterpolated Uninterpolated field to add to interpolated field.
    */
    VectorField_LinearInterpolation_OF
    (Field&& field, Locator&& locator,
     Uninterpolated&& uninterpolated,
     CheckOptions::Check)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    , _uninterpolated{ uninterpolated }
    {
      static_assert(check_if_outside == true && warn_if_outside == false,
                    "Bad class template arguments for bounds checking.");
    }
    
    /** Constructor.
     \brief With warning if particle is in bounds.
     \param field Scalar field cell data to interpolate.
     \param locator Mesh locator.
     \param uninterpolated Uninterpolated field to add to interpolated field.
    */
    VectorField_LinearInterpolation_OF
    (Field&& field, Locator&& locator,
     Uninterpolated&& uninterpolated,
     CheckOptions::Warn)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    , _uninterpolated{ uninterpolated }
    {
      static_assert(check_if_outside == true && warn_if_outside == true,
                    "Bad class template arguments for bounds warning.");
    }
    
    /**
     \brief Interpolate field.
     \param position 3D position.
     \param cell Hint for mesh cell index position is in.
     \return interpolated field value.
    */
    auto operator()(Point const& position, Index cell = -1) const
    {
      cell = _locator(position, cell);
      if constexpr (check_if_outside)
        if (outside<warn_if_outside>(cell, position, "Assigning zero."))
          return Vector::zero;
      
      if constexpr (!std::is_same_v<Uninterpolated, useful::Empty>)
        return _interpolant.interpolate(position, cell) + _uninterpolated[cell];
      
      return _interpolant.interpolate(position, cell);
    }
    
    /**
     \brief Interpolate field.
     \param position 2D position.
     \param cell Hint for mesh cell index position is in.
     \return interpolated field value.
    */
    auto operator()(Point2D const& position, Index cell = -1) const
    {
      auto interp = (*this)(make_point(position), cell);
      
      return Vector2D{ interp[0], interp[1] };
    }
  
    /**
     \brief Interpolate field.
     \param position 1D position.
     \param cell Hint for mesh cell index position is in.
     \return interpolated field value.
    */
    auto operator()(Scalar position, Index cell = -1) const
    {
      return (*this)(make_point(position), cell)[0];
    }
    
    /**
     \brief Interpolate field.
     \param state Particle state to interpolate.
    */
    template <typename State>
    auto operator()
    (State const& state) const
    {
      return (*this)(state.position, state.cell);
    }
    
    /**
     \brief Locate a state in the mesh.
     \param state Particle state to locate.
     \return Mesh cell index.
    */
    template <typename State>
    auto locate(State const& state) const
    { return _locator(state); }
    
    /** \return Underlying field data. */
    auto const& field() const
    { return _field; }
    
    /** \return Underlying uninterpolated data. */
    auto const& uninterpolated() const
    { return _uninterpolated; }
    
    /** \return Locator object to find positions in mesh. */
    auto const& locator() const
    { return _locator; }
    
    /** \brief Rescale underlying field.
     \param factor Scaling factor.
    */
    void rescale(double factor)
    { _field *= factor; }
    
    /** \brief Rescale underlying uninterpolated field.
     \param factor Scaling factor.
    */
    void rescale_uninterpolated(double factor)
    { _uninterpolated *= factor; }
    
    /**
     \brief Sum to underlying field.
     \param field field to sum
    */
    void sum(Field const& field)
    { _field += field; }
    
    /**
     \brief Sum to underlying uninterpolated field.
     \param uninterpolated field to sum
    */
    void sum_uninterpolated(Uninterpolated const& uninterpolated)
    { _uninterpolated += uninterpolated; }
    
    /**
     \brief Change underlying field data.
     \param field Set field data to this field.
    */
    void set(Field const& field)
    { _field = field; }
    
    /**
     \brief Change underlying uninterpolated field data.
     \param uninterpolated Set uninterpolated field data to this field.
    */
    void set_uninterpolated(Uninterpolated const& uninterpolated)
    { _uninterpolated = uninterpolated; }
    
  private:
    Field _field;                     /**< Vector field cell data to interpolate. */
    Locator _locator;                 /**< Locator to find positions in mesh. */
    Foam::interpolationCellPoint<Vector> _interpolant{ _field };  /**< Interpolation object. */
    Uninterpolated _uninterpolated;   /** Uninterpolated field to add to interpolated field. */
  };
  template
  <typename Field, typename Locator>
  VectorField_LinearInterpolation_OF
  (Field&&, Locator&&,
   CheckOptions::NoCheck) ->
  VectorField_LinearInterpolation_OF
  <Field, Locator, 0, 0, useful::Empty>;
  template
  <typename Field, typename Locator>
  VectorField_LinearInterpolation_OF
  (Field&&, Locator&&,
   CheckOptions::Check) ->
  VectorField_LinearInterpolation_OF
  <Field, Locator, 1, 0, useful::Empty>;
  template
  <typename Field, typename Locator>
  VectorField_LinearInterpolation_OF
  (Field&&, Locator&&,
   CheckOptions::Warn) ->
  VectorField_LinearInterpolation_OF
  <Field, Locator, 1, 1, useful::Empty>;
  template
  <typename Field, typename Locator, typename Uninterpolated>
  VectorField_LinearInterpolation_OF
  (Field&&, Locator&&, Uninterpolated&&,
   CheckOptions::NoCheck) ->
  VectorField_LinearInterpolation_OF
  <Field, Locator, 0, 0, Uninterpolated>;
  template
  <typename Field, typename Locator, typename Uninterpolated>
  VectorField_LinearInterpolation_OF
  (Field&&, Locator&&, Uninterpolated&&,
   CheckOptions::Check) ->
  VectorField_LinearInterpolation_OF
  <Field, Locator, 1, 0, Uninterpolated>;
  template
  <typename Field, typename Locator, typename Uninterpolated>
  VectorField_LinearInterpolation_OF
  (Field&&, Locator&&, Uninterpolated&&,
   CheckOptions::Warn) ->
  VectorField_LinearInterpolation_OF
  <Field, Locator, 1, 1, Uninterpolated>;
  
  /** \class ScalarField_LinearInterpolation_OF PTOF/Field.h "PTOF/Field.h"
   * \brief Linear interpolation of scalar field cell data based on OpenFOAM routines. */
  template
  <typename Field, typename Locator,
  bool check_if_outside, bool warn_if_outside,
  typename Uninterpolated = useful::Empty>
  class ScalarField_LinearInterpolation_OF
  {
  public:
    using Point = Foam::point;                      /**> 3D point. */
    using Point2D = Foam::Vector2D<Foam::scalar>;   /**> 2D point. */
    using Scalar = Foam::scalar;                    /**> Scalar (also 1D point). */
    using Index = Foam::label;                      /**> Cell index.*/
    
    /** Constructor.
     \brief Without checking if position is in bounds.
     \param field Scalar field cell data to interpolate.
     \param locator Mesh locator.
     */
    ScalarField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, CheckOptions::NoCheck)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    {
      static_assert(check_if_outside == false && warn_if_outside == false,
                    "Bad class template arguments for no bounds checking.");
    }
    
    /** Constructor.
     \brief With checking if position is in bounds but no warning if not.
     \param field Scalar field cell data to interpolate.
     \param locator Mesh locator.
    */
    ScalarField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, CheckOptions::Check)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    {
      static_assert(check_if_outside == true && warn_if_outside == false,
                    "Bad class template arguments for bounds checking.");
    }
    
    /** Constructor.
     \brief With warning if particle is in bounds.
     \param field Scalar field cell data to interpolate.
     \param locator Mesh locator.
    */
    ScalarField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, CheckOptions::Warn)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    {
      static_assert(check_if_outside == true && warn_if_outside == true,
                    "Bad class template arguments for bounds warning.");
    }
    
    /** Constructor.
    \brief Without checking if position is in bounds.
    \param field Vector field cell data to interpolate.
    \param locator Mesh locator.
    \param uninterpolated Uninterpolated field to add to interpolated field.
    */
    ScalarField_LinearInterpolation_OF
    (Field&& field, Locator&& locator,
     Uninterpolated&& uninterpolated,
     CheckOptions::NoCheck)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    , _uninterpolated{ uninterpolated }
    {
      static_assert(check_if_outside == false && warn_if_outside == false,
                    "Bad class template arguments for no bounds checking.");
    }
    
    /** Constructor.
     \brief With checking if position is in bounds but no warning if not.
     \param field Vector field cell data to interpolate.
     \param locator Mesh locator.
     \param uninterpolated Uninterpolated field to add to interpolated field.
    */
    ScalarField_LinearInterpolation_OF
    (Field&& field, Locator&& locator,
     Uninterpolated&& uninterpolated,
     CheckOptions::Check)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    , _uninterpolated{ uninterpolated }
    {
      static_assert(check_if_outside == true && warn_if_outside == false,
                    "Bad class template arguments for bounds checking.");
    }
    
    /** Constructor.
     \brief With warning if particle is in bounds.
     \param field Scalar field cell data to interpolate.
     \param locator Mesh locator.
     \param uninterpolated Uninterpolated field to add to interpolated field.
    */
    ScalarField_LinearInterpolation_OF
    (Field&& field, Locator&& locator,
     Uninterpolated&& uninterpolated,
     CheckOptions::Warn)
    : _field{ std::forward<Field>(field) }
    , _locator{ std::forward<Locator>(locator) }
    , _uninterpolated{ uninterpolated }
    {
      static_assert(check_if_outside == true && warn_if_outside == true,
                    "Bad class template arguments for bounds warning.");
    }
    
    /**
     \brief Interpolate field.
     \param position 3D position.
     \param cell Hint for mesh cell index position is in.
     \return interpolated field value.
    */
    auto operator()(Point const& position, Index cell = -1) const
    {
      cell = _locator(position, cell);
      if constexpr (check_if_outside)
        if (outside<warn_if_outside>(cell, position, "Assigning zero."))
          return 0.;
      
      if constexpr (!std::is_same_v<Uninterpolated, useful::Empty>)
        return _interpolant.interpolate(position, cell) + _uninterpolated[cell];
      
      return _interpolant.interpolate(position, cell);
    }
    
    /**
     \brief Interpolate field.
     \param position 2D position.
     \param cell Hint for mesh cell index position is in.
     \return interpolated field value.
    */
    auto operator()(Point2D const& position, Index cell = -1) const
    {
      return (*this)(make_point(position), cell);
    }
    
    /**
     \brief Interpolate field.
     \param position 1D position.
     \param cell Hint for mesh cell index position is in.
     \return interpolated field value.
    */
    auto operator()(Scalar position, Index cell = -1) const
    {
      return (*this)(make_point(position), cell);
    }
    
    /**
     \brief Interpolate field.
     \param state Particle state to interpolate.
    */
    template <typename State>
    auto operator()
    (State const& state) const
    {
      return (*this)(state.position, state.cell);
    }
    
    /**
     \brief Locate a state in the mesh.
     \param state Particle state to locate.
     \return Mesh cell index.
    */
    template <typename State>
    auto locate(State const& state) const
    { return _locator(state); }
    
    /** \return Underlying field data. */
    auto const& field() const
    { return _field; }
    
    /** \return Underlying uninterpolated data. */
    auto const& uninterpolated() const
    { return _uninterpolated; }
    
    /** \return Locator object to find positions in mesh. */
    auto const& locator() const
    { return _locator; }
    
    /** \brief Rescale underlying field.
     \param factor Scaling factor.
    */
    void rescale(double factor)
    { _field *= factor; }
    
    /** \brief Rescale underlying uninterpolated field.
     \param factor Scaling factor.
    */
    void rescale_uninterpolated(double factor)
    { _uninterpolated *= factor; }
    
    /**
     \brief Sum to underlying field.
     \param field field to sum
    */
    void sum(Field const& field)
    { _field += field; }
    
    /**
     \brief Sum to underlying uninterpolated field.
     \param uninterpolated field to sum
    */
    void sum_uninterpolated(Uninterpolated const& uninterpolated)
    { _uninterpolated += uninterpolated; }
    
    /**
     \brief Change underlying field data.
     \param field Set field data to this field.
    */
    void set(Field const& field)
    { _field = field; }
    
    /**
     \brief Change underlying uninterpolated field data.
     \param uninterpolated Set uninterpolated field data to this field.
    */
    void set_uninterpolated(Uninterpolated const& uninterpolated)
    { _uninterpolated = uninterpolated; }
    
  private:
    Field _field;      /**< Scalar field cell data to interpolate.         */
    Locator _locator;  /**< Locator to find positions in mesh. */
    Foam::interpolationCellPoint<Scalar> _interpolant{ _field }; /**< Interpolation object. */
    Uninterpolated _uninterpolated;   /** Uninterpolated field to add to interpolated field. */
  };
  template
  <typename Field, typename Locator>
  ScalarField_LinearInterpolation_OF
  (Field&&, Locator&&,
   CheckOptions::NoCheck) ->
  ScalarField_LinearInterpolation_OF
  <Field, Locator, 0, 0, useful::Empty>;
  template
  <typename Field, typename Locator>
  ScalarField_LinearInterpolation_OF
  (Field&&, Locator&&,
   CheckOptions::Check) ->
  ScalarField_LinearInterpolation_OF
  <Field, Locator, 1, 0, useful::Empty>;
  template
  <typename Field, typename Locator>
  ScalarField_LinearInterpolation_OF
  (Field&&, Locator&&,
   CheckOptions::Warn) ->
  ScalarField_LinearInterpolation_OF
  <Field, Locator, 1, 1, useful::Empty>;
  template
  <typename Field, typename Locator, typename Uninterpolated>
  ScalarField_LinearInterpolation_OF
  (Field&&, Locator&&, Uninterpolated&&,
   CheckOptions::NoCheck) ->
  ScalarField_LinearInterpolation_OF
  <Field, Locator, 0, 0, Uninterpolated>;
  template
  <typename Field, typename Locator, typename Uninterpolated>
  ScalarField_LinearInterpolation_OF
  (Field&&, Locator&&, Uninterpolated&&,
   CheckOptions::Check) ->
  ScalarField_LinearInterpolation_OF
  <Field, Locator, 1, 0, Uninterpolated>;
  template
  <typename Field, typename Locator, typename Uninterpolated>
  ScalarField_LinearInterpolation_OF
  (Field&&, Locator&&, Uninterpolated&&,
   CheckOptions::Warn) ->
  ScalarField_LinearInterpolation_OF
  <Field, Locator, 1, 1, Uninterpolated>;
  
  /** \brief Compute magnitude of average of volumetric field (set of values associated with mesh cells). */
  template <typename Field, typename Mesh>
  auto magnitude_of_average(Field& field, Mesh const& mesh)
  {
    Foam::scalar mesh_volume
      = Foam::sum(mesh.cellVolumes());
    auto average_weighted_data
      = Foam::sum(field*mesh.cellVolumes());
    return Foam::mag(average_weighted_data)/mesh_volume;
  }
  
  /** \brief Rescale a volumetric field (set of values associated with mesh cells)  to a given average value. */
  template <typename Field, typename Mesh>
  void rescale_to_average
  (Field& field, Mesh const& mesh, double average)
  {
    field *= average/magnitude_of_average(field, mesh);
  }
}


#endif /* PTOF_FIELD_H */
