//
//  Field.h
//
//  Created by Tomás Aquino on 16/02/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Field_OF_h
#define Field_OF_h

#include <fieldTypes.H>
#include <interpolationCellPoint.H>
#include <iostream>
#include "general/useful.h"
#include "ParticleTrackingOF/useful.h"


namespace ptof
{
  // Options for interpolation constructors
  namespace FieldOptions
  {
    using NoCheck = useful::Selector<int, 0>; // No bounds checking
    using Check = useful::Selector<int, 1>; // Bounds checking, no warning
    using Warn = useful::Selector<int, 2>; // Bounds checking, warning
  }
  
  // Linear interpolation of vector field data
  // based on OpenFOAM routines
  template
  <typename Field, typename Locator,
  bool check_if_outside, bool warn_if_outside>
  class VectorField_LinearInterpolation_OF
  {
  public:
    // Typedefs for position, vector, and mesh cell index types
    using Point = Foam::point;
    using Point2D = Foam::Vector2D<Foam::scalar>;
    using Vector = Foam::vector;
    using Vector2D = Foam::Vector2D<Foam::scalar>;
    using Scalar = Foam::scalar;
    using Index = Foam::label;
    
    // Construct without bounds checking and without bounds warning
    VectorField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, FieldOptions::NoCheck)
    : field{ std::forward<Field>(field) }
    , locator{ std::forward<Locator>(locator) }
    {}
    
    // Construct with bounds checking and without bounds warning
    VectorField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, FieldOptions::Check)
    : field{ std::forward<Field>(field) }
    , locator{ std::forward<Locator>(locator) }
    {}
    
    // Construct with bounds checking and with bounds warning
    VectorField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, FieldOptions::Warn)
    : field{ std::forward<Field>(field) }
    , locator{ std::forward<Locator>(locator) }
    {}
    
    // Interpolate field according to 3D position and mesh cell
    auto operator()(Point const& position, Index cell) const
    {
      if constexpr (check_if_outside)
        if (check_bounds(position, cell))
          return Vector::zero;
      
      return interpolant.interpolate(position, cell);
    }
    
    // Interpolate field according to 3D position
    auto operator()(Point const& position) const
    {
      return this->operator()(position, locator(position));
    }
    
    // Interpolate field according to 2D position and mesh cell
    auto operator()(Point2D const& position, Index cell) const
    {
      auto interp = this->operator()(make_point(position), cell);
      
      return Vector2D{ interp[0], interp[1] };
    }
    
    // Interpolate field according to 2D position
    auto operator()(Point2D const& position) const
    {
      return this->operator()(position, locator(position));
    }
    
    // Interpolate field according to 1D position and mesh cell
    auto operator()(Scalar position, Index cell) const
    {
      return this->operator()(make_point(position), cell)[0];
    }
    
    // Interpolate field according to scalar position
    auto operator()(Scalar position) const
    {
      return this->operator()(position, locator(position));
    }
    
    // Interpolate field according to state
    template <typename State>
    auto operator()
    (State const& state) const
    {
      return this->operator()(state.position, locator(state));
    }
    
    // Locate a state in the mesh
    template <typename State>
    auto locate(State const& state) const
    { return locator(state); }
    
    // Get underlying field
    auto const& get_field() const
    { return field; }
    
    // Get underlying locator object
    auto const& get_locator() const
    { return locator; }
    
    // Rescale underlying field
    auto rescale(double factor)
    { field *= factor; }
    
  private:
    Field field;      // Field data to interpolate
    Locator locator;  // Locator to find positions in mesh
    
    // Interpolation object
    Foam::interpolationCellPoint<Vector> interpolant{ field };
    
    //Bounds checking
    bool check_bounds
    (Point const& position, Index cell) const
    {
      if (cell == -1)
      {
        if constexpr (warn_if_outside)
          std::cerr
            << "Warning: Requested field value at position "
            << "("
            << position[0] << ", "
            << position[1] << ", "
            << position[2] << ")"
            << " outside mesh; assigning zero\n";
        return 1;
      }
      return 0;
    }
  };
  template
  <typename Field, typename Locator>
  VectorField_LinearInterpolation_OF
  (Field&&, Locator&&,
   FieldOptions::NoCheck) ->
  VectorField_LinearInterpolation_OF
  <Field, Locator, 0, 0>;
  template
  <typename Field, typename Locator>
  VectorField_LinearInterpolation_OF
  (Field&&, Locator&&,
   FieldOptions::Check) ->
  VectorField_LinearInterpolation_OF
  <Field, Locator, 1, 0>;
  template
  <typename Field, typename Locator>
  VectorField_LinearInterpolation_OF
  (Field&&, Locator&&,
   FieldOptions::Warn) ->
  VectorField_LinearInterpolation_OF
  <Field, Locator, 1, 1>;
  
  // Linear interpolation of scalar field data
  // based on OpenFOAM routines
  template
  <typename Field, typename Locator,
  bool check_if_outside, bool warn_if_outside>
  class ScalarField_LinearInterpolation_OF
  {
  public:
    // Typedefs for position, vector, and mesh cell index types
    using Point = Foam::point;
    using Point2D = Foam::Vector2D<Foam::scalar>;
    using Scalar = Foam::scalar;
    using Index = Foam::label;
    
    // Construct without bounds checking and without bounds warning
    ScalarField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, FieldOptions::NoCheck)
    : field{ std::forward<Field>(field) }
    , locator{ std::forward<Locator>(locator) }
    {}
    
    // Construct with bounds checking and without bounds warning
    ScalarField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, FieldOptions::Check)
    : field{ std::forward<Field>(field) }
    , locator{ std::forward<Locator>(locator) }
    {}
    
    // Construct with bounds checking and with bounds warning
    ScalarField_LinearInterpolation_OF
    (Field&& field, Locator&& locator, FieldOptions::Warn)
    : field{ std::forward<Field>(field) }
    , locator{ std::forward<Locator>(locator) }
    {}
    
    // Interpolate field according to 3D position and mesh cell
    auto operator()(Point const& position, Index cell) const
    {
      if constexpr (check_if_outside)
        if (check_bounds(position, cell))
          return 0.;
      
      return interpolant.interpolate(position, cell);
    }
    
    // Interpolate field according to 3D position
    auto operator()(Point const& position) const
    {
      return this->operator()(position, locator(position));
    }
    
    // Interpolate field according to 2D position and mesh cell
    auto operator()(Point2D const& position, Index cell) const
    {
      return this->operator()(make_point(position), cell);
    }
    
    // Interpolate field according to 2D position
    auto operator()(Point2D const& position) const
    {
      return this->operator()(position, locator(position));
    }
    
    // Interpolate field according to 1D position and mesh cell
    auto operator()(Scalar position, Index cell) const
    {
      return this->operator()(make_point(position), cell);
    }
    
    // Interpolate field according to scalar position
    auto operator()(Scalar position) const
    {
      return this->operator()(position, locator(position));
    }
    
    // Interpolate field according to state
    template <typename State>
    auto operator()
    (State const& state) const
    {
      return this->operator()(state.position, locator(state));
    }
    
    // Locate a state in the mesh
    template <typename State>
    auto locate(State const& state) const
    { return locator(state); }
    
    // Get underlying field
    auto const& get_field() const
    { return field; }
    
    // Get underlying locator object
    auto const& get_locator() const
    { return locator; }
    
    // Rescale underlying field
    auto rescale(double factor)
    { field *= factor; }
    
  private:
    Field field;      // Field data to interpolate
    Locator locator;  // Locator to find positions in mesh
    
    // Interpolation object
    Foam::interpolationCellPoint<Scalar> interpolant{ field };
    
    //Bounds checking
    bool check_bounds
    (Point const& position, Index cell) const
    {
      if (cell == -1)
      {
        if constexpr (warn_if_outside)
          std::cerr
            << "Warning: Requested field value at position "
            << "("
            << position[0] << ", "
            << position[1] << ", "
            << position[2] << ")"
            << " outside mesh; assigning zero\n";
        return 1;
      }
      return 0;
    }
  };
  template
  <typename Field, typename Locator>
  ScalarField_LinearInterpolation_OF
  (Field&&, Locator&&,
   FieldOptions::NoCheck) ->
  ScalarField_LinearInterpolation_OF
  <Field, Locator, 0, 0>;
  template
  <typename Field, typename Locator>
  ScalarField_LinearInterpolation_OF
  (Field&&, Locator&&,
   FieldOptions::Check) ->
  ScalarField_LinearInterpolation_OF
  <Field, Locator, 1, 0>;
  template
  <typename Field, typename Locator>
  ScalarField_LinearInterpolation_OF
  (Field&&, Locator&&,
   FieldOptions::Warn) ->
  ScalarField_LinearInterpolation_OF
  <Field, Locator, 1, 1>;
  
  // Compute magnitude of average of volumetric field
  // (set of values associated with mesh cells)
  template <typename Field, typename Mesh>
  auto magnitude_of_average(Field& field, Mesh const& mesh)
  {
    Foam::scalar mesh_volume
      = Foam::sum(mesh.cellVolumes());
    auto average_weighted_data
      = Foam::sum(field*mesh.cellVolumes());
    return Foam::mag(average_weighted_data)/mesh_volume;
  }
  
  // Rescale a volumetric field
  // (set of values associated with mesh cells)
  // to a given average value
  template <typename Field, typename Mesh>
  void rescale_to_average
  (Field& field, Mesh const& mesh, double average)
  {
    field *= average/magnitude_of_average(field, mesh);
  }
}


#endif /* Field_OF_h */
