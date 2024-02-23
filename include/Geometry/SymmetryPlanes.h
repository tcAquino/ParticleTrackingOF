/**
 \file Geometry/SymmetryPlanes.h
 \author Tomás Aquino
 \date 13/05/2020
 
 \brief Symmetry planes to describe general symmetries.
 
 A SymmetryPlanes object must define:
 - dim : The spatial dimension.
 - length : The length of the unit cell along the normal vector associated with each plane.
 - normal : A vector of unit normal vectors associated with each plane.
 - translation : A vector of translation vectors associated with each plane.
*/

#ifndef SymmetryPlanes_h
#define SymmetryPlanes_h

#include <cmath>
#include "General/Operations.h"

namespace geometry
{
  /** \struct SymmetryPlanes_Bcc Geometry/SymmetryPlanes.h "Geometry/SymmetryPlanes.h"
   \brief Symmetry plane for body centered cubic sphere packing. */
  struct SymmetryPlanes_Bcc
  {
    /** Constructor.
    \brief Symmetry planes for unit sphere radius. */
    SymmetryPlanes_Bcc()
    : normal{
        operation::times_scalar(1./std::sqrt(2.),
                                std::vector<double>{ 0., 1., 1. }),
        operation::times_scalar(1./std::sqrt(2.),
                                std::vector<double>{ 1., 0., 1. }),
        operation::times_scalar(1./std::sqrt(2.),
                                std::vector<double>{ 1., 1., 0. })
    }
    , translation{
        operation::times_scalar(2./std::sqrt(3.),
                                std::vector<double>{ -1., 1., 1. }),
        operation::times_scalar(2./std::sqrt(3.),
                                std::vector<double>{ 1., -1., 1. }),
        operation::times_scalar(2./std::sqrt(3.),
                                std::vector<double>{ 1., 1., -1. })
    }
    {}
    
    static constexpr std::size_t dim{ 3 };              /**< Spatial dimension. */
    const std::vector<double> length{
      std::vector<double>(dim, 4./std::sqrt(6.)) };     /**< Primitive cell side. */
    const std::vector<std::vector<double>> normal;      /**< Unit normal vectors associated with each plane. */
    const std::vector<std::vector<double>> translation; /**< Translation vectors associated with each plane.*/
  };
  
  /** \struct SymmetryPlanes_Sc Geometry/SymmetryPlanes.h "Geometry/SymmetryPlanes.h"
  \brief Symmetry plane for simple cubic sphere packing. */
  struct SymmetryPlanes_Sc
  {
    /** Constructor.
    \brief Symmetry planes for unit sphere radius. */
    SymmetryPlanes_Sc()
    : normal{
        { 1., 0., 0. },
        { 0., 1., 0. },
        { 0., 0., 1. }
    }
    , translation{
        std::vector<double>{ 2., 0., 0. },
        std::vector<double>{ 0., 2., 0. },
        std::vector<double>{ 0., 0., 2. }
    }
    {}
    
    static constexpr std::size_t dim{ 3 };                /**< Spatial dimension. */
    const std::vector<double> length{
      std::vector<double>(dim, 1.) };                     /**< Primitive cell side. */
    const std::vector<std::vector<double>> normal;        /**< Unit normal vectors associated with each plane. */
    const std::vector<std::vector<double>> translation;   /**< Translation vectors associated with each plane.*/
  };
  
  /** \brief Compute component of projection along basis vectors of symmetry planes.
   \param position Position to project.
   \param symmetry_planes Symmetry planes defining the periodic unit cell.
   \param dd Component of project to compute.
   \param scale Scale factor.
   \return Projection component.
   */
  template <typename Position, typename SymmetryPlanes>
  auto project
  (Position const& position, SymmetryPlanes const& symmetry_planes,
   std::size_t dd, double scale)
  {
    return operation::dot(position, symmetry_planes.normal[dd])
      /(scale*symmetry_planes.length[dd]);
  }
  
  /** \brief Compute projection along basis vectors of symmetry planes
  \param position Position to project.
  \param symmetry_planes Symmetry planes defining the periodic unit cell.
  \param scale Scale factor.
  \return Projection vector.
  */
  template <typename Position, typename SymmetryPlanes>
  auto project
  (Position const& position, SymmetryPlanes const& symmetry_planes,
   double scale)
  {
    Position projections;
    projections.reserve(symmetry_planes.dim);
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
      projections.push_back(project(position, symmetry_planes,
                                    dd, scale));
    return projections;
  }
  
  /** \brief Compute component of projection along basis vectors of symmetry planes.
   \param position Position to project.
   \param symmetry_planes Symmetry planes defining the periodic unit cell.
   \param dd Component of project to compute.
   \param scale Scale factor.
   \param origin Coordinate system origin.
   \return Projection component.
  */
  template <typename Position, typename SymmetryPlanes, typename Origin>
  auto project
  (Position const& position, SymmetryPlanes const& symmetry_planes,
   std::size_t dd, double scale, Origin const& origin)
  {
    return project(operation::minus(position, origin),
                   symmetry_planes, dd, scale);
  }
  
  /** \brief Compute projection along basis vectors of symmetry planes.
   \param position Position to project.
   \param symmetry_planes Symmetry planes defining the periodic unit cell.
   \param scale Scale factor.
   \param origin Coordinate system origin.
   \return Projection.
  */
  template <typename Position, typename SymmetryPlanes, typename Origin>
  auto project
  (Position const& position, SymmetryPlanes const& symmetry_planes,
   double scale, Origin const& origin)
  {
    return project(operation::minus(position, origin),
                   symmetry_planes, scale);
  }
  
  /** \brief Translate position along basis vectors of symmetry planes.
   \param position Position to translate.
   \param symmetry_planes Symmetry planes defining the periodic unit cell.
   \param projections Pre-computed projections along basis vectors.
   \param scale Scale factor.
  */
  template <typename Position, typename SymmetryPlanes,
  typename Projections = std::vector<double>>
  void translate
  (Position& position, SymmetryPlanes const& symmetry_planes,
   Projections const& projections,
   double scale)
  {
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
      operation::plus_InPlace(position,
        operation::times_scalar(scale*projections[dd],
                                symmetry_planes.translation[dd]));
  }
  
  /** \brief Translate position backwards along basis vectors of symmetry planes.
   \param position Position to translate.
   \param symmetry_planes Symmetry planes defining the periodic unit cell.
   \param projections Pre-computed projections along basis vectors.
   \param scale Scale factor.
  */
  template <typename Position, typename SymmetryPlanes,
  typename Projections = std::vector<double>>
  void translate_back
  (Position& position, SymmetryPlanes const& symmetry_planes,
   Projections const& projections,
   double scale)
  {
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
      operation::minus_InPlace(position,
        operation::times_scalar(scale*projections[dd],
                                symmetry_planes.translation[dd]));
  }
  
  /** \brief Place particle in periodic unit cell according to symmetry planes.
   \param position Position to place in unit cell.
   \param symmetry_planes Symmetry planes defining the periodic unit cell.
   \param scale Scale factor.
   \return Projection along basis vectors of symmetry planes.
  */
  template <typename Position, typename SymmetryPlanes>
  auto place_in_unit_cell
  (Position& position, SymmetryPlanes const& symmetry_planes,
   double scale)
  {
    std::vector<int> projections;
    projections.reserve(symmetry_planes.dim);
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
      projections.push_back(std::floor(project(position, symmetry_planes,
                                               dd, scale)));
    translate_back(position, symmetry_planes, projections, scale);

    return projections;
  }
  
  /** \brief Place particle in periodic unit cell according to symmetry planes.
   \param position Position to place in unit cell.
   \param symmetry_planes Symmetry planes defining the periodic unit cell.
   \param scale Scale factor.
   \param origin Coordinate system origin.
   \return Projection along basis vectors of symmetry planes.
  */
  template <typename Position, typename SymmetryPlanes, typename Origin>
  auto place_in_unit_cell
  (Position& position, SymmetryPlanes const& symmetry_planes,
   double scale, Origin const& origin)
  {
    std::vector<int> projections;
    projections.reserve(symmetry_planes.dim);
    for (std::size_t dd = 0; dd < symmetry_planes.dim; ++dd)
      projections.push_back(std::floor(project(position, symmetry_planes,
                                               dd, scale, origin)));
    translate_back(position, symmetry_planes, projections, scale);

    return projections;
  }
  
  /** \brief Add periodic velocity point images outside domain.
   \details If the projection of a point is within a given fraction of the boundary along each symmetry plane, add an image along that plane.
   \param points Point positions, to add new points.
   \param velocities Velocity values, to add new velocities.
   \param fraction Fraction of original points to add.
   \param boundary_periodic Object to enforce periodic boundary conditions.
   \note \c Boundary_Periodic must define:
   - symmetry_planes
   - scale
   - origin
   - translate(Positions&, Projections const&)
  */
  template
  <typename Points = std::vector<std::vector<double>>,
   typename Velocities = std::vector<std::vector<double>>,
   typename Boundary_Periodic>
  void add_periodic_images
  (Points& points, Velocities& velocities,
   double fraction, Boundary_Periodic const& boundary_periodic)
  {
    std::size_t nr_points = points.size();
    if (nr_points == 0)
      return;
    
    std::size_t dim = points.back().size();
    std::vector<std::vector<int>>
      projections_forward(dim, std::vector<int>(dim));
    std::vector<std::vector<int>>
      projections_backward(dim, std::vector<int>(dim));
    for (std::size_t dd = 0; dd < dim; ++dd)
    {
      projections_forward[dd][dd] = 1.;
      projections_backward[dd][dd] = -1.;
    }
    
    for (std::size_t pp = 0; pp < nr_points; ++pp)
    {
      auto projections =
        geometry::project(points[pp],
                          boundary_periodic.symmetry_planes,
                          boundary_periodic.scale,
                          boundary_periodic.origin);
      for (std::size_t dd = 0; dd < dim; ++dd)
      {
        // If close to lower bound, make an image forward
        if (projections[dd] < fraction)
        {
          points.push_back(points[pp]);
          velocities.push_back(velocities[pp]);
          boundary_periodic.translate(points.back(),
                                      projections_forward[dd]);
        }
        // If close to upper bound, make an image backward
        if (projections[dd] > 1.-fraction)
        {
          points.push_back(points[pp]);
          velocities.push_back(velocities[pp]);
          boundary_periodic.translate(points.back(),
                                      projections_backward[dd]);
        }
      }
    }
  }
}


#endif /* SymmetryPlanes_h */
