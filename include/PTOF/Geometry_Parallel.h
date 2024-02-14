/**
* \file PTOF/Geometry_Parallel.h
* \author Tomás Aquino
* \date 13/02/2024
*/
#ifndef PTOF_GEOMETRY_PARALLEL_H
#define PTOF_GEOMETRY_PARALLEL_H

#include <algorithm>
#include <cstddef>
#include <string>
#include <iterator>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <fvMesh.H>
#include <IOobject.H>
#include <meshSearch.H>
#include <unordered_map>
#include "Geometry/Boundary.h"
#include "PTOF/Boundary.h"
#include "PTOF/Directories.h"
#include "PTOF/Geometry.h"

namespace ptof
{
  /** \struct Geometry_Parallel PTOF/Geometry_Parallel.h "PTOF/Geometry_Parallel.h"
   * Basic geometry class */
  template
  <std::size_t dim_val,
  BoundaryConditionSet::Type dynamics>
  struct Geometry_Parallel
  {
    static constexpr std::size_t dim{ dim_val };
    using Mesh = Foam::fvMesh;
    using MeshSearch = Foam::meshSearch;
    using Locator = Locator_Cell;
    
    /** Store only absorptions for transport simulations,
     * store absorptions and reinjections for fpt simulations. */
    using BoundaryInfo
      = std::conditional_t<
          dynamics
            == BoundaryConditionSet::Type::transport,
          Store_Absorbed,
          Store_Absorbed_Reinjections>;
    
    /** Construct with mesh from specified OpenFOAM case time. */
    template <typename Parameters>
    Geometry_Parallel
    (DirectoriesOF const& directories_of,
     Directories const& directories,
     std::size_t num_threads,
     Parameters const& parameters = useful::Empty{})
    : mesh{
        Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ } }
    {
      mesh_copies.reserve(num_threads-1);
      for (std::size_t thread = 0; thread < num_threads-1; ++thread)
        mesh_copies.emplace_back(std::make_unique<Mesh>(Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ }));
      mesh_search.reserve(num_threads);
      mesh_search.emplace_back(std::make_unique<MeshSearch>(mesh));
      for (std::size_t thread = 1; thread < num_threads; ++thread)
        mesh_search.emplace_back(std::make_unique<MeshSearch>(*mesh_copies[thread-1]));
      locator.reserve(num_threads);
      for (std::size_t thread = 0; thread < num_threads; ++thread)
        locator.emplace_back(*mesh_search[thread]);
    }
    
    /** Make a Boundary object for the preset type of dynamics:
     * - 'transport': No parameters, custom boundary does nothing;
     * - 'firstpassage': Parameters should hold InitialCondition object;
     * - custom boundary reinjects according to initial condition. */
    template
    <typename TransportParameters,
    typename ReactionParameters,
    typename SolverParameters,
    typename SurfaceReaction,
    typename Parameters = useful::Empty>
    auto makeBoundary
    (Directories const& directories,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     SurfaceReaction&& reaction,
     std::size_t thread,
     Parameters&& parameters = {}) const
    {
      if constexpr (dynamics == BoundaryConditionSet::Type::transport)
        return ptof::Boundary_Cases{
          get_boundary_conditions<dynamics>(directories),
          *mesh_search[thread],
          locator[thread],
          BoundaryInfo{},
          Boundary_DoNothing{},
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(reaction) };
      if constexpr (dynamics == BoundaryConditionSet::Type::firstpassage)
        return ptof::Boundary_Cases{
          get_boundary_conditions<dynamics>(directories),
          *mesh_search[thread],
          locator[thread],
          BoundaryInfo{},
          Boundary_DoNothing{},
          Boundary_Reinject{ parameters } };
      throw std::runtime_error{
        std::string{ "Boundary conditions for " }
          + BoundaryConditionSet::name(dynamics)
          + " not supported" };
    }
    
    /** Output generic information about object. */
    template <typename OStream>
    static void info(OStream& output)
    {
      output <<
        "--------------------------------------------------\n"
        "Geometry\n"
        "--------------------------------------------------\n"
        "Spatial dimension: " + std::to_string(dim) + "\n"
        "Boundary conditions for: " + BoundaryConditionSet::name(dynamics) + "\n"
        "Boundary condition types:\n"
        "\treflecting: Reflecting\n"
        "\tperiodic: No effect\n"
        "\tabsorbing: Absorbing\n"
        "\tinfo: Information upon crossing\n"
        "\tcustom: ";
      output <<
        (BoundaryConditionSet::name(dynamics) == "transport"
         ? "No effect\n"
         : "Reinject according to initial condition\n");
      output <<
        "\tempty: No effect (default for unspecified patches in mesh)\n"
        "--------------------------------------------------\n";
    }
    
    /** Output information about current object. */
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Mesh\n"
        "--------------------------------------------------\n";
      output << "cells: " << mesh.nCells() << "\n"
             << "faces: " << mesh.nFaces() << "\n"
             << "edges: " << mesh.nEdges() << "\n";
      output <<
        "--------------------------------------------------\n";
    }
    
    Mesh mesh;                                               /**< Mesh for fields and boundaries. */
    std::vector<std::unique_ptr<Mesh>> mesh_copies;          /**< Mesh copies needed for parallel searches */
    std::vector<std::unique_ptr<MeshSearch>> mesh_search;         /**< To find positions in mesh, etc. */
    std::vector<Locator> locator;                            /**< Wrapper on mesh search.         */
  };
  
  /** \struct Geometry_Periodic_Cartesian_Parallel PTOF/Geometry_Parallel.h "PTOF/Geometry_Parallel.h"
   *  \brief Geometry class for geometry
   *  with Cartesian periodic BCs along some dimensions.*/
  template
  <std::size_t dim_val,
  BoundaryConditionSet::Type dynamics>
  struct Geometry_Periodic_Cartesian_Parallel
  {
    static constexpr std::size_t dim{ dim_val };
    using Mesh = Foam::fvMesh;
    using MeshSearch = Foam::meshSearch;
    using Locator = Locator_Cell;
    using Boundary_Periodic = geometry::Boundary_Periodic_WithOutsideInfo;
    
    /** Store only absorptions. */
    using BoundaryInfo = Store_Absorbed;
    
    /** Construct with mesh from specified OpenFOAM case time. */
    template <typename Parameters>
    Geometry_Periodic_Cartesian_Parallel
    (DirectoriesOF const& directories_of,
     Directories const& directories,
     std::size_t num_threads,
     Parameters const& params = useful::Empty{})
    : mesh{
        Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ } }
    , boundary_periodic{
        make_boundary_periodic(directories) }
    {
      mesh_copies.reserve(num_threads-1);
      for (std::size_t thread = 0; thread < num_threads-1; ++thread)
        mesh_copies.emplace_back(std::make_unique<Mesh>(Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ }));
      mesh_search.reserve(num_threads);
      mesh_search.emplace_back(std::make_unique<MeshSearch>(mesh));
      for (std::size_t thread = 1; thread < num_threads; ++thread)
        mesh_search.emplace_back(std::make_unique<MeshSearch>(*mesh_copies[thread-1]));
      locator.reserve(num_threads);
      for (std::size_t thread = 0; thread < num_threads; ++thread)
        locator.emplace_back(*mesh_search[thread]);
    }
    
    /** Make a Boundary object for Cartesian periodicity. */
    template
    <typename TransportParameters,
    typename ReactionParameters,
    typename SolverParameters,
    typename SurfaceReaction,
    typename Parameters = useful::Empty>
    auto makeBoundary
    (Directories const& directories,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     SurfaceReaction&& reaction,
     std::size_t thread,
     Parameters&& parameters = {}) const
    {
      if constexpr (dynamics == BoundaryConditionSet::Type::transport)
        return ptof::Boundary_Cases{
          get_boundary_conditions<dynamics>(directories),
          *mesh_search[thread],
          locator[thread],
          BoundaryInfo{},
          Boundary_Periodic_OF{
            boundary_periodic,
            *mesh_search[thread] },
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(reaction) };
      if constexpr (dynamics == BoundaryConditionSet::Type::firstpassage)
        return ptof::Boundary_Cases{
          get_boundary_conditions<dynamics>(directories),
          *mesh_search[thread],
          locator[thread],
          BoundaryInfo{},
          Boundary_Periodic_OF{
            boundary_periodic,
            *mesh_search[thread] },
          Boundary_Reinject{ parameters } };
      throw std::runtime_error{
        std::string{ "Boundary conditions for " }
          + BoundaryConditionSet::name(dynamics)
          + " not supported" };
    }
    
    /** Output generic information about object. */
    template <typename OStream>
    static void info(OStream& output)
    {
      output <<
        "--------------------------------------------------\n"
        "Geometry\n"
        "--------------------------------------------------\n"
        "Spatial dimension: " + std::to_string(dim) + "\n"
        "Boundary conditions for: " + BoundaryConditionSet::name(dynamics) + "\n"
        "Boundary condition types:\n"
        "\treflecting: Reflecting\n"
        "\tperiodic: Cartesian periodic, position extracted from mesh\n";
      output <<
        "\tabsorbing: Absorbing\n"
        "\tinfo: Information upon crossing\n"
        "\tcustom: ";
      output <<
        (BoundaryConditionSet::name(dynamics) == "transport"
         ? "No effect\n"
         : "Reinject according to initial condition\n");
      output <<
        "\tempty: No effect (default for unspecified patches in mesh)\n"
      "--------------------------------------------------\n";
    }
    
    /** Output information about current object. */
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Mesh\n"
        "--------------------------------------------------\n";
      output << "cells: " << mesh.nCells() << "\n"
             << "faces: " << mesh.nFaces() << "\n"
             << "edges: " << mesh.nEdges() << "\n"
        "--------------------------------------------------\n";
      output <<
        "--------------------------------------------------\n"
        "Periodic boundaries\n"
        "--------------------------------------------------\n";
      output << "Automatically extracted Cartesian periodic boundaries:\n";
      if (boundary_periodic.boundaries.size() == 0)
        output << "None";
      for (std::size_t dd = 0;
           dd < boundary_periodic.boundaries.size(); ++dd)
        if (boundary_periodic.boundaries[dd].second
            < std::numeric_limits<double>::infinity())
          output << "Dimension " << dd
                 << " at " << boundary_periodic.boundaries[dd].first
                 << " and " << boundary_periodic.boundaries[dd].second
                 << "\n";
      output <<
        "--------------------------------------------------\n";
    }
    
    Mesh mesh;                                             /**< Mesh for fields and boundaries.        */
    Boundary_Periodic boundary_periodic;                   /**< Boundary object to enforce periodicity.*/
    std::vector<std::unique_ptr<Mesh>> mesh_copies;          /**< Mesh copies needed for parallel searches */
    std::vector<std::unique_ptr<MeshSearch>> mesh_search;  /**< To find positions in mesh, etc.        */
    std::vector<Locator> locator;                          /**< Wrapper on mesh search.                */
    
  private:
    Boundary_Periodic make_boundary_periodic
    (Directories const& directories)
    {
      auto boundary_conditions =
        get_boundary_conditions<dynamics>(directories);
      
      std::vector<std::vector<double>> boundaries_dim(dim);
      for (auto const& bc : boundary_conditions)
        if (bc.second == "periodic")
        {
          /** Find average and variance of patch face centers. */
          auto const& patch_id =
            mesh.boundaryMesh().findPatchID(bc.first, false);
          auto const& patch = mesh.boundaryMesh()[patch_id];
          auto face_start = patch.start();
          Foam::point average_face_center = Foam::point::zero;
          Foam::point variance_face_center = Foam::point::zero;
          Foam::scalar total_area = 0.;
          for (auto face = face_start;
               face < face_start + patch.size(); ++face)
          {
            auto area = face_area(face, mesh);
            auto weighted_center = area*face_center(face, mesh);
            total_area += area;
            average_face_center += weighted_center;
            for (std::size_t dd = 0; dd < variance_face_center.size(); ++dd)
              variance_face_center[dd] +=
                weighted_center[dd]*weighted_center[dd];
          }
          average_face_center /= total_area;
          variance_face_center /= total_area*total_area;
          for (std::size_t dd = 0; dd < variance_face_center.size(); ++dd)
            variance_face_center[dd] -=
              average_face_center[dd]*average_face_center[dd];
          
          /** Choose dimension with least variance. */
          auto it_min = std::min_element(variance_face_center.begin(),
                                         variance_face_center.end());
          std::size_t dd = std::distance(variance_face_center.begin(),
                                         it_min);
          boundaries_dim[dd].push_back(average_face_center[dd]);
        }
      
      /** Check which dimensions have periodic BCs
       * If a dimension has periodic BCs they must be exactly two.*/
      std::vector<bool> dim_has_periodic_bcs(dim, 0);
      for (std::size_t dd = 0; dd < dim; ++dd)
        if (boundaries_dim[dd].size() != 0)
        {
          if (boundaries_dim[dd].size() != 2)
            throw std::runtime_error{
              "Extracted "
              + std::to_string(boundaries_dim[dd].size())
              + " periodic boundaries "
              "along cartesian dimension "
              + std::to_string(dd) };
          
          dim_has_periodic_bcs[dd] = 1;
        }
      
      /** Periodic conditions are assumed to be checked in order
       * across Cartesian dimensions.
       * If there are periodic BCs in all dimensions
       * starting from 0 up to some dimension, add only those boundaries
       * If there are gap dimensions between periodic dimensions,
       * add {-inf, inf} boundaries.*/
      bool found_nonperiodic = 0;
      bool nonconsecutive_dims = 0;
      for (bool periodic : dim_has_periodic_bcs)
      {
        if (found_nonperiodic && periodic)
        {
          nonconsecutive_dims = 1;
          break;
        }
        if (!periodic)
          found_nonperiodic = 1;
      }
      std::vector<std::pair<double, double>> boundaries;
      for (std::size_t dd = 0; dd < dim; ++dd)
      {
        if (dim_has_periodic_bcs[dd])
        {
          boundaries.push_back({
            std::min(boundaries_dim[dd][0],
                     boundaries_dim[dd][1]),
            std::max(boundaries_dim[dd][0],
                     boundaries_dim[dd][1])
          });
        }
        else if (nonconsecutive_dims)
        {
          boundaries.push_back({
            -std::numeric_limits<double>::infinity(),
            std::numeric_limits<double>::infinity()
          });
        }
      }
      return { boundaries };
    };
  };
  
  /** \struct Geometry_Bcc_Parallel PTOF/Geometry_Parallel.h "PTOF/Geometry_Parallel.h"
   * \brief Geometry class for periodic geometry in a
   * body centered cubic beadpack
   * Holds spatial dimension info,
   * mesh, mesh search objects,
   * and periodic boundary to enforce periodicity. */
  template
  <std::size_t dim_val,
  BoundaryConditionSet::Type dynamics,
  BoundaryConditionSet::Type periodicity_type>
  struct Geometry_Bcc_Parallel
  {
    static constexpr std::size_t dim{ dim_val };
    using Mesh = Foam::fvMesh;
    using MeshSearch = Foam::meshSearch;
    using Locator = Locator_Cell;
    using Boundary_Periodic = std::conditional_t<
      periodicity_type
      == BoundaryConditionSet::Type::cartesian,
    geometry::Boundary_Periodic_WithOutsideInfo,
    geometry::Boundary_Periodic_SymmetryPlanes_WithOutsideInfo<
      geometry::SymmetryPlanes_Bcc>>;
    
    /** Store only absorptions.*/
    using BoundaryInfo = Store_Absorbed;
    
    /** Construct with mesh from specified OpenFOAM case time. */
    template <typename Parameters>
    Geometry_Bcc_Parallel
    (DirectoriesOF const& directories_of,
     Directories const& directories,
     std::size_t num_threads,
     Parameters const& params)
    : mesh{
        Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ } }
    , boundary_periodic{
        make_boundary_periodic(useful::Selector<
                                BoundaryConditionSet::Type,
                                periodicity_type>{},
                               params) }
    {
      mesh_copies.reserve(num_threads-1);
      for (std::size_t thread = 0; thread < num_threads-1; ++thread)
        mesh_copies.emplace_back(std::make_unique<Mesh>(Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ }));
      mesh_search.reserve(num_threads);
      mesh_search.emplace_back(std::make_unique<MeshSearch>(mesh));
      for (std::size_t thread = 1; thread < num_threads; ++thread)
        mesh_search.emplace_back(std::make_unique<MeshSearch>(*mesh_copies[thread-1]));
      locator.reserve(num_threads);
      for (std::size_t thread = 0; thread < num_threads; ++thread)
        locator.emplace_back(*mesh_search[thread]);
    }
    
    /** Make a Boundary object for the preset type of periodicity
     * Parameters holds periodic boundary
     * - 'cartesian': Cartesian periodicity
     * - 'symmetryplanes': Periodicity according to symmetry planes */
    template
    <typename TransportParameters,
    typename ReactionParameters,
    typename SolverParameters,
    typename SurfaceReaction,
    typename Parameters = useful::Empty>
    auto makeBoundary
    (Directories const& directories,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     SurfaceReaction&& reaction,
     std::size_t thread,
     Parameters&& parameters = {}) const
    {
      if constexpr (dynamics == BoundaryConditionSet::Type::transport)
        return ptof::Boundary_Cases{
          get_boundary_conditions<dynamics>(directories),
          *mesh_search[thread],
          locator[thread],
          BoundaryInfo{},
          Boundary_Periodic_OF{
            boundary_periodic,
            *mesh_search[thread] },
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(reaction) };
      if constexpr (dynamics == BoundaryConditionSet::Type::firstpassage)
        return ptof::Boundary_Cases{
          get_boundary_conditions<dynamics>(directories),
          *mesh_search[thread],
          locator[thread],
          BoundaryInfo{},
          Boundary_Periodic_OF{
            boundary_periodic,
            *mesh_search[thread] },
          Boundary_Reinject{ parameters } };
      throw std::runtime_error{
        std::string{ "Boundary conditions for " }
          + BoundaryConditionSet::name(dynamics)
          + " not supported" };
    }
    
    /** Output generic information about object.*/
    template <typename OStream>
    static void info(OStream& output)
    {
      output <<
        "--------------------------------------------------\n"
        "Geometry\n"
        "--------------------------------------------------\n"
        "Spatial dimension: " + std::to_string(dim) + "\n"
        "Boundary conditions for: " + BoundaryConditionSet::name(dynamics) + "\n"
        "Boundary condition types:\n"
        "\treflecting: Reflecting\n"
        "\tperiodic: ";
      output <<
        (BoundaryConditionSet::name(periodicity_type) == "cartesian"
         ? "Periodic in the primitive unit cell\n"
         : "Periodic in the minimal unit cell\n");
      output <<
        "\tabsorbing: Absorbing\n"
        "\tinfo: Information upon crossing\n"
        "\tcustom: ";
      output <<
        (BoundaryConditionSet::name(dynamics) == "transport"
         ? "No effect\n"
         : "Reinject according to initial condition\n");
      output <<
        "\tempty: No effect (default for unspecified patches in mesh)\n"
      "--------------------------------------------------\n";
    }
    
    /** Output information about current object.*/
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Mesh\n"
        "--------------------------------------------------\n";
      output << "cells: " << mesh.nCells() << "\n"
             << "faces: " << mesh.nFaces() << "\n"
             << "edges: " << mesh.nEdges() << "\n";
      output <<
        "--------------------------------------------------\n";
    }
    
    Mesh mesh;                                             /**< Mesh for fields and boundaries.        */
    Boundary_Periodic boundary_periodic;                   /**< Boundary object to enforce periodicity.*/
    std::vector<std::unique_ptr<Mesh>> mesh_copies;          /**< Mesh copies needed for parallel searches */
    std::vector<std::unique_ptr<MeshSearch>> mesh_search;  /**< To find positions in mesh, etc.        */
    std::vector<Locator> locator;                          /**< Wrapper on mesh search.                */
  
  private:
    template <typename Parameters>
    Boundary_Periodic make_boundary_periodic
    (useful::Selector<BoundaryConditionSet::Type,
      BoundaryConditionSet::Type::cartesian>,
      Parameters const& params)
    {
      return { params.primitive_cell_boundaries };
    };
  
    template <typename Parameters>
    Boundary_Periodic make_boundary_periodic
    (useful::Selector<BoundaryConditionSet::Type,
      BoundaryConditionSet::Type::symmetryplanes>,
     Parameters const& params)
    {
      std::vector<double> midpoint;
      for (auto const& boundary
           : params.primitive_cell_boundaries)
        midpoint.push_back( (boundary.first + boundary.second)/2. );
      return { params.radius, midpoint };
    }
  };
}

#endif /* PTOF_GEOMETRY_PARALLEL_H */
