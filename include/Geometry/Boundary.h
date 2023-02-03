/**
* \file Geometry/Boundary.h
* \author Tomás Aquino
* \date 08/02/2019
*/

#ifndef GEOMETRY_BOUNDARY_H
#define GEOMETRY_BOUNDARY_H

#include <cmath>
#include <numeric>
#include <utility>
#include <vector>
#include "General/Constants.h"
#include "General/Operations.h"
#include "Geometry/SymmetryPlanes.h"

// A boundary class should implement the following basic functionality:
//
// class Example_Boundary
// {
// public:
//
//   // Check if position is out of bounds
//   template <typename Position>
//   bool outOfBounds(Position const& position) const
//   {
//     // Return true if boundary condition is to be applied, false otherwise
//   }
//
//   // Enforce boundary condition
//   template <typename State>
//   bool operator()(State& state, State const& state_old) const
//   {
//     // Apply the boundary condition given current and previous state
//     // Return true if anything was done, false otherwise
//   }
// };

namespace geometry
{
  /** Return displacement to be added to 1d position to place it within
   periodic boundaries. */
	double boundary_periodic
  (double position, std::pair<double, double> const& boundaries)
	{
		double box_size = boundaries.second - boundaries.first;
		double pos = position - boundaries.first;
		return -std::floor(pos/box_size)*box_size;
	}
  
  /** Return pair of displacement to be added to 1d position
   * and the associated signed number of domain displacements out
   * to place it within periodic boundaries. */
  std::pair<double, int> boundary_periodic_with_outside_info
  (double position, std::pair<double, double> const& boundaries)
  {
    double box_size = boundaries.second - boundaries.first;
    double pos = position - boundaries.first;
    int outside = int(std::floor(pos/box_size));
    return { -outside*box_size, outside };
  }

  /** Return displacement to be added to 1d position
   * to reflect it off boundaries. */
	double boundary_reflecting
  (double position, std::pair<double, double> const& boundaries)
	{
		double box_size = boundaries.second - boundaries.first;
		double pos = position - boundaries.first;
		double nr_boxes_jumped = std::trunc(pos / box_size);
		if (int(nr_boxes_jumped) % 2 == 0)
			return -pos + std::abs(pos - nr_boxes_jumped * box_size);
		else
			return -pos + box_size - std::abs(pos - nr_boxes_jumped * box_size);
	}

  /** Return whether 1d position is outside boundaries. */
  bool outOfBounds_box(double position, std::pair<double,double> const& boundaries)
  {
    return position < boundaries.first || position > boundaries.second;
  }

  /** Return whether position is outside boundaries along each dimension. */
  template <typename Position = std::vector<double>,
  typename Boundaries = std::vector<std::pair<double, double>>>
  bool outOfBounds_box(Position const& position, Boundaries const& boundaries)
  {
    position[0];
    for (std::size_t dd=0; dd<position.size(); ++dd)
      if (outOfBounds_box(position[dd], boundaries[dd]))
        return true;

    return false;
  }

  /** \class Boundary_DoNothing Geometry/Boundary.h "Geometry/Boundary.h"
   *  \brief Boundary with no effect. */
	class Boundary_DoNothing
	{
	public:
    /** Check if position is out of bounds (it never is). */
		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return false; }

    /** Enforce boundary condition (do nothing). */
		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    { return false; }
	};

  /** \class Boundary_Periodic_1d Geometry/Boundary.h "Geometry/Boundary.h"
   * \brief Periodic boundary in 1d./
   * 
   * State must define: 
   * - position (scalar type) */
	class Boundary_Periodic_1d
	{
	public:
		const std::pair<double, double> boundaries;  /**< Lower and upper boundaries. */

    /** Construct given lower and upper boundaries. */
		Boundary_Periodic_1d(std::pair<double, double> boundaries)
		: boundaries(boundaries)
		{}

    /** Check if position is out of bounds.*/
		template <typename Position>
		bool outOfBounds(Position position) const
		{ return OutOfBounds_Box(position, boundaries); }

    /** Enforce boundary condition. */
		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
		{
      if (!outOfBounds(state.position))
        return false;
      state.position += boundary_periodic(state.position, boundaries);
      return true;
    }
    
    template <typename Position, typename Projection = double>
    void translate(Position& position, Projection projection) const
    {
      operation::plus_InPlace(position,
                      (boundaries.second-boundaries.first)*projection);
    }
	};

  /** \class Boundary_Periodic Geometry/Boundary.h "Geometry/Boundary.h"
   * \brief  Periodic boundaries along each dimension.
   * 
   * State must define: 
   * - position (vector type) */
	class Boundary_Periodic
	{
	public:
    /** Lower and upper boundaries along each dimension. */
		const std::vector<std::pair<double, double>> boundaries;
    const std::vector<double> domain_dimensions;

    /** Construct given lower and upper boundaries along each dimension. */
		Boundary_Periodic(std::vector<std::pair<double, double>> boundaries)
    : boundaries{ boundaries }
    , domain_dimensions{ make_dimensions(boundaries) }
		{}

    /** Check if position is out of bounds. */
		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return outOfBounds_box(position, boundaries); }

    /** Enforce boundary condition. */
		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
		{
      if (!outOfBounds(state.position))
        return false;
      for (std::size_t dd = 0; dd < boundaries.size(); ++dd)
        state.position[dd] += boundary_periodic(state.position[dd], boundaries[dd]);
      return true;
		}
    
    /** Translate position according to symmetry planes. */
    template <typename Position, typename Projections = std::vector<double>>
    void translate
    (Position& position, Projections const& projections) const
    {
      operation::plus_InPlace(position,
                      operation::times(domain_dimensions, projections));
    }
    
  private:
    std::vector<double> make_dimensions
    (std::vector<std::pair<double, double>> const& boundaries) const
    {
      std::vector<double> dimensions;
      for (auto const& val : boundaries)
        dimensions.push_back(val.second - val.first);
        
      return dimensions;
    }
	};
  
  /** \class Boundary_Periodic_Dim Geometry/Boundary.h "Geometry/Boundary.h"
   * \brief Periodic boundary along dimension dd.
   * 
   * State must define: 
   * - position (vector type) */
  template <std::size_t dd>
  class Boundary_Periodic_Dim
  {
  public:
    
    const std::pair<double, double> boundaries;  /**< Lower and upper boundaries along dimension dd. */
    static const std::size_t dim = dd;           /**< Dimension for outside visibility. */

    /** Construct with boundaries at same distance from origin given domain halfwidth. */
    Boundary_Periodic_Dim(double half_width)
    : boundaries{ -half_width, half_width }
    {}

    /** Construct given lower and upper boundaries along dimension dim. */
    Boundary_Periodic_Dim(std::pair<double, double> boundaries)
    : boundaries{ boundaries }
    {}

    /** Check if position is out of bounds. */
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return outOfBounds_box(position[dd], boundaries); }

    /** Enforce boundary condition. */
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      if (!outOfBounds(state.position))
        return false;
      state.position[dd] +=
        boundary_periodic(state.position[dd], boundaries);
      return true;
    }
  };
  
  /** \class Boundary_Periodic_WithOutsideInfo Geometry/Boundary.h "Geometry/Boundary.h"
   * \brief Periodic boundaries along each dimension
   * with information about where position would be outside domain.
   * 
   * State must define: 
   * - position (vector type)
   * - periodicity (std::vector<int>) counting how many domains have been traveled along each dimension */
  class Boundary_Periodic_WithOutsideInfo
  {
  public:
    const std::vector<std::pair<double, double>>boundaries;   /**< Lower and upper boundaries along each dimension. */
    const std::vector<double> domain_dimensions;

    /** Construct given lower and upper boundaries. */
    Boundary_Periodic_WithOutsideInfo
    (std::vector<std::pair<double, double>> boundaries)
    : boundaries(boundaries)
    , domain_dimensions{ make_dimensions(boundaries) }
    {}

    /** Check if position is out of bounds. */
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return outOfBounds_box(position, boundaries); }

    /** Enforce boundary condition. */
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      if (!outOfBounds(state.position))
        return false;
      for (std::size_t dd = 0; dd < boundaries.size(); ++dd)
      {
        auto change_outside =
          boundary_periodic_with_outside_info(state.position[dd], boundaries[dd]);
        state.position[dd] += change_outside.first;
        state.periodicity[dd] += change_outside.second;
      }
      return true;
    }
    
    /** Translate position according to symmetry planes. */
    template <typename Position, typename Projections = std::vector<double>>
    void translate
    (Position& position, Projections const& projections) const
    {
      operation::plus_InPlace(position,
                              operation::times(domain_dimensions, projections));
    }
    
  private:
    std::vector<double> make_dimensions
    (std::vector<std::pair<double, double>> const& boundaries) const
    {
      std::vector<double> dimensions;
      for (auto const& val : boundaries)
        dimensions.push_back(val.second - val.first);
        
      return dimensions;
    }
  };

  /** \class Boundary_Reflecting_1d Geometry/Boundary.h "Geometry/Boundary.h"
   * \brief Reflecting boundary in one dimension.
   * 
   * State must define: 
   * - position (scalar type) */
	class Boundary_Reflecting_1d
	{
	public:
		const std::pair<double, double> boundaries;  /**< Lower and upper boundaries. */

    /** Construct given lower and upper boundaries. */
		Boundary_Reflecting_1d(std::pair<double, double> boundaries)
		: boundaries(boundaries)
		{}

    /** Check if position is out of bounds. */
		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return outOfBounds_box(position, boundaries); }

    /** Enforce boundary condition. */
		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
		{
      if (!outOfBounds(state.position))
        return false;
      state.position += boundary_reflecting(state.position, boundaries);
      return true;
    }
	};

  /** \class Boundary_Reflecting Geometry/Boundary.h "Geometry/Boundary.h"
   * \brief Reflecting boundaries along along each dimension. 
   * 
   * State must define: 
   * - position (vector type) */
	class Boundary_Reflecting
	{
	public:
		const std::vector<std::pair<double, double>>boundaries; /**< Lower and upper boundaries along each dimension. */

    /** Construct given lower and upper boundaries along each dimension. */
		Boundary_Reflecting(std::vector<std::pair<double, double>> boundaries)
		: boundaries(boundaries)
		{}

    /** Check if position is out of bounds. */
		template <typename Position>
		bool outOfBounds(Position const& position) const
		{ return outOfBounds_box(position, boundaries); }

    /** Enforce boundary condition. */
		template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
		{
      if (!outOfBounds(state.position))
        return false;
			for (std::size_t dd = 0; dd < boundaries.size(); ++dd)
				state.position[dd] +=
          boundary_reflecting(state.position[dd], boundaries[dd]);
      return true;
		}
	};

  /** \class Boundary_Reflecting_Dim Geometry/Boundary.h "Geometry/Boundary.h"
   * \brief Reflecting boundary along dimension dd.
   * 
   * State must define: 
   * - position (vector type) */
  template <std::size_t dd>
  class Boundary_Reflecting_Dim
  {
  public:
    
    const std::pair<double, double> boundaries;  /**< Lower and upper boundaries along dimension dd. */
    static const std::size_t dim = dd;           /**< Dimension for outside visibility.*/    

    /** Construct with boundaries at same distance from origin given domain halfwidth. */
    Boundary_Reflecting_Dim(double half_width)
    : boundaries{ -half_width, half_width }
    {}

    /** Construct given lower and upper boundaries along dimension dim. */
    Boundary_Reflecting_Dim(std::pair<double, double> boundaries)
    : boundaries{ boundaries }
    {}

    /** Check if position is out of bounds. */
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return outOfBounds_box(position[dd], boundaries); }

    /** Enforce boundary condition. */
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      if (!outOfBounds(state.position))
        return false;
      state.position[dd] +=
        boundary_reflecting(state.position[dd], boundaries);
      return true;
    }
  };
  
  /** \class Boundary_RadialReflecting_2d Geometry/Boundary.h "Geometry/Boundary.h"
   * \brief Reflecting boundaries on the inside of a circle. 
   * 
   * State must define: 
   * - position (vector type) */
  class Boundary_RadialReflecting_2d
  {
  public:
    const double radius;

    /** Construct given domain radius along reflecting dimensions. */
    Boundary_RadialReflecting_2d(double radius)
    : radius{ radius }
    {}

    /** Check if position is out of bounds. */
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      return operation::abs_sq(position) > radius_sq;
    }

    /** Enforce boundary condition. */
    template <typename State>
    bool operator()(State& state, State const& state_old) const
    {
      if (!outOfBounds(state.position))
        return false;
      
      auto position_old = state_old.position;
      std::vector<double> jump = operation::minus(state.position, position_old);
      std::vector<double> tangent_at_contact(2);
      
      // Iterate reflection procedure until position is inside boundary
      while (1)
      {
        // Fraction alpha of jump up to the boundary
        double pos_dot_jump = operation::dot(position_old, jump);
        double norm_sq_jump = operation::abs_sq(jump);
        double pos_sq_old = operation::abs_sq(position_old);
        double aux = pos_dot_jump/norm_sq_jump;
        double alpha = std::sqrt(aux*aux
                                 +(radius_sq-pos_sq_old)/norm_sq_jump)-aux;
        
        // Avoid numerical issues
        if (alpha <= 0. || useful::isnan(alpha))
        {
          operation::times_scalar_InPlace(radius/operation::abs(state.position),
                                          state.position);
          break;
        }
        
        // Place old position at contact
        for (std::size_t dd = 0; dd < 2; ++dd)
          position_old[dd] += alpha*jump[dd];

        // Tangent to boundary at contact
        tangent_at_contact[0] = position_old[1]/radius;
        tangent_at_contact[1] = -position_old[0]/radius;

        // Angle between jump and tangent to boundary at contact
        double tangent_dot_jump = operation::dot(tangent_at_contact, jump);
        double cos_contact_angle = tangent_dot_jump/std::sqrt(norm_sq_jump);
        double contact_angle = std::acos(cos_contact_angle);

        // Outgoing jump at collision, rotate outer jump and flip the direction
        operation::rotate(jump, -2.*contact_angle);
        operation::times_scalar_InPlace(1.-alpha, jump);
        
        // Place particle at reflected position
        operation::plus(position_old, jump,
                        state.position);
        
        // If inside boundary, done
        if(!outOfBounds(state.position))
          break;
      }
      
      return true;
    }

  private:
    const double radius_sq{ radius*radius };
  };

  /** \class Boundary_Open_RadialReflecting_3d "Geometry/Boundary.h"
  * \brief Reflecting boundaries on the inside of an infinite cylinder in 3d.
  * 
  * Cylinder longitudinal axis is dd_open, which may be 0 or 2. \n
  * State must define: 
  * - position (vector type) */
  template <std::size_t dd_open = 0>
  class Boundary_Open_RadialReflecting_3d
  {
    static_assert(dd_open == 0 || dd_open == 2,
                  "Cylinder axis dimension must be 0 or 2");
  public:
    const double radius;

    /** Construct given domain radius along reflecting dimension. */
    Boundary_Open_RadialReflecting_3d(double radius)
    : radius{ radius }
    {}

    /** Check if position is out of bounds. */
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      double radial_pos_sq = std::inner_product(
        position.cbegin()+begin_transverse, position.cbegin()+begin_transverse+2,
        position.cbegin()+begin_transverse, 0.);
      return radial_pos_sq > radius_sq;
    }

    /** Enforce boundary condition/ */
    template <typename State>
    bool operator()(State& state, State const& state_old) const
    {
      if (!outOfBounds(state.position))
        return false;
      
      std::vector<double> position_radial_old(2);
      std::vector<double> jump_radial(2);
      std::vector<double> tangent_at_contact(2);
      
      for (std::size_t dd = 0; dd < 2; ++dd)
      {
        jump_radial[dd] = state.position[dd+begin_transverse]
          -state_old.position[dd+begin_transverse];
        position_radial_old[dd] = state_old.position[dd+begin_transverse];
      }
      
      // Iterate reflection procedure until position is inside boundary
      while (1)
      {
        // Fraction alpha of jump up to the boundary
        double pos_dot_jump = std::inner_product(
          position_radial_old.cbegin(), position_radial_old.cend(),
          jump_radial.cbegin(), 0.);
        double norm_sq_jump = std::inner_product(
          jump_radial.cbegin(), jump_radial.cend(),
          jump_radial.cbegin(), 0.);
        double radial_pos_sq_old = operation::abs_sq(position_radial_old);
        double aux = pos_dot_jump/norm_sq_jump;
        double alpha = std::sqrt(aux*aux
                                 +(radius_sq-radial_pos_sq_old)/norm_sq_jump)-aux;
        
        // Avoid numerical issues
        if (alpha <= 0. || useful::isnan(alpha))
        {
           double radial_pos_sq = std::inner_product(
             state.position.cbegin()+begin_transverse,
             state.position.cbegin()+begin_transverse+2,
             state.position.cbegin()+begin_transverse, 0.);
           for (std::size_t dd = 0; dd < 2; ++dd)
             state.position[dd+begin_transverse] *= radius/std::sqrt(radial_pos_sq);
           break;
        }

        // Place old position at contact
        for (std::size_t dd = 0; dd < 2; ++dd)
          position_radial_old[dd] += alpha*jump_radial[dd];

        // Tangent to boundary at contact
        tangent_at_contact[0] = position_radial_old[1]/radius;
        tangent_at_contact[1] = -position_radial_old[0]/radius;

        // Angle between jump and tangent to boundary at contact
        double tangent_dot_jump = std::inner_product(
          tangent_at_contact.cbegin(), tangent_at_contact.cend(),
          jump_radial.cbegin(), 0.);
        double cos_contact_angle = tangent_dot_jump/std::sqrt(norm_sq_jump);
        double contact_angle = std::acos(cos_contact_angle);

        // Outgoing jump at collision, rotate outer jump and flip the direction
        operation::rotate(jump_radial, -2.*contact_angle);
        operation::times_scalar_InPlace(1.-alpha, jump_radial);
        
        // Place particle at reflected position
        for (std::size_t dd = 0; dd < 2; ++dd)
          state.position[dd+begin_transverse] =
            position_radial_old[dd]+jump_radial[dd];
        
        // If inside boundary, done
        if(!outOfBounds(state.position))
          break;
      }
      
      return true;
    }

  private:
    const double radius_sq{ radius*radius };
    const std::size_t begin_transverse{ dd_open == 0 ? 1 : 0 };
  };
  
  /** \class Boundary_Periodic_SymmetryPlanes Geometry/Boundary.h "Geometry/Boundary.h"
   * \brief Periodic boundaries along each symmetry plane. 
   * 
   * State must define:
   * - position (vector type) */
  template <typename SymmetryPlanes>
  class Boundary_Periodic_SymmetryPlanes
  {
  public:
    const SymmetryPlanes symmetry_planes;
    const double scale;
    const std::vector<double> origin;
    
    /** Construct given symmetry plane object,
     * overall scale factor, and coordinate origin. */
    Boundary_Periodic_SymmetryPlanes
    (SymmetryPlanes symmetry_planes, double scale = 1.,
     std::vector<double> origin = {})
    : symmetry_planes{ symmetry_planes }
    , scale{ scale }
    , origin{ origin.size() == 0.
      ? std::vector<double>(symmetry_planes.dim, 0.)
      : origin }
    {}
    
    /** Construct given overall scale factor and coordinate origin. */
    Boundary_Periodic_SymmetryPlanes
    (double scale = 1., std::vector<double> origin = {})
    : symmetry_planes{}
    , scale{ scale }
    , origin{ origin.size() == 0.
      ? std::vector<double>(symmetry_planes.dim, 0.)
      : origin }
    {}

    /** Check if position is out of bounds. */
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
        if (std::floor(geometry::project(position,
                                         symmetry_planes, dd,
                                         scale, origin)) != 0.)
          return 1;
      return 0;
    }

    /** Enforce boundary condition. */
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      auto projections = geometry::place_in_unit_cell(state.position,
                                                      symmetry_planes,
                                                      scale, origin);
      
      for (auto const& val : projections)
        if (val != 0)
          return 1;
      
      return 0;
    }
    
    /** Translate position according to symmetry planes. */
    template <typename Position, typename Projections = std::vector<double>>
    void translate
    (Position& position, Projections const& projections) const
    {
      geometry::translate(position, symmetry_planes, projections, scale);
    }
  };
  
  /** \class Boundary_Periodic_SymmetryPlanes_WithOutsideInfo Geometry/Boundary.h "Geometry/Boundary.h"
   * \brief Periodic boundaries along each symmetry plane
   * with information about where position would be outside domain.
   * 
   * State must define:
   * - position (vector type)
   * - periodicity (std::vector<int>) counting how many domains have been traveled along each dimension*/
  template <typename SymmetryPlanes>
  class Boundary_Periodic_SymmetryPlanes_WithOutsideInfo
  {
  public:
    const SymmetryPlanes symmetry_planes;
    const double scale;
    const std::vector<double> origin;
    
    /** Construct given symmetry plane object,
     * overall scale factor, and coordinate origin. */
    Boundary_Periodic_SymmetryPlanes_WithOutsideInfo
    (SymmetryPlanes symmetry_planes, double scale = 1.,
     std::vector<double> origin = {})
    : symmetry_planes{ symmetry_planes }
    , scale{ scale }
    , origin{ origin.size() == 0.
      ? std::vector<double>(symmetry_planes.dim, 0.)
      : origin }
    {}
    
    /** Construct given overall scale factor and coordinate origin. */
    Boundary_Periodic_SymmetryPlanes_WithOutsideInfo
    (double scale = 1., std::vector<double> origin = {})
    : symmetry_planes{}
    , scale{ scale }
    , origin{ origin.size() == 0.
      ? std::vector<double>(symmetry_planes.dim, 0.)
      : origin }
    {}
    
    /** Check if position is out of bounds. */
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
        if (std::floor(geometry::project(position,
                                         symmetry_planes, dd,
                                         scale, origin)) != 0.)
          return 1;
      return 0;
    }

    /** Enforce boundary condition. */
    template <typename State>
    bool operator()(State& state, State const& state_old = {}) const
    {
      auto projections =
        geometry::place_in_unit_cell(state.position,
                                     symmetry_planes,
                                     scale, origin);
      operation::plus_InPlace(state.periodicity, projections);
      
      for (auto const& val : projections)
        if (val != 0)
          return 1;
      
      return 0;
    }
    
    /** Translate position according to symmetry planes. */
    template <typename Position, typename Projections = std::vector<double>>
    void translate
    (Position& position, Projections const& projections) const
    {
      geometry::translate(position, symmetry_planes, projections, scale);
    }
  };
}

#endif /* GEOMETRY_BOUNDARY_H */
