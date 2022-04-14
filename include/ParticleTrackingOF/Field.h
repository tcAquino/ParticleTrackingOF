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
#include "general/useful.h"
#include "ParticleTrackingOF/useful.h"

namespace ptof
{
  // Linear interpolation of field data
  // based on OpenFOAM routines
  template
  <typename Field, typename Locator>
  class Field_LinearInterpolation_OF
  {
  public:
    // Typedefs for position, vector, and mesh cell index types
    using Point = Foam::point;
    using Point2D = Foam::Vector2D<Foam::scalar>;
    using Vector = Foam::vector;
    using Vector2D = Foam::Vector2D<Foam::scalar>;
    using Scalar = Foam::scalar;
    using Index = Foam::label;
    
    // Construct
    Field_LinearInterpolation_OF
    (Field&& field, Locator&& locator)
    : field{ std::forward<Field>(field) }
    , locator{ std::forward<Locator>(locator) }
    {}
    
    // Interpolate field according to 3D position and mesh cell
    auto operator()(Point const& position, Index cell) const
    {
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
      auto interp = interpolant.interpolate(make_point(position), cell);
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
      return interpolant.interpolate(make_point(position), cell)[0];
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
    
  private:
    Field field;      // Field data to interpolate
    Locator locator;  // Locator to find positions in mesh
    
    // Interpolation object
    Foam::interpolationCellPoint<Foam::vector> interpolant{ field };
  };
  template
  <typename Field, typename Locator>
  Field_LinearInterpolation_OF
  (Field&&, Locator&&) ->
  Field_LinearInterpolation_OF
  <Field, Locator>;
}


#endif /* Field_OF_h */
