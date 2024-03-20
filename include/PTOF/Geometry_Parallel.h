/**
 \file PTOF/Geometry_Parallel.h
 \author Tomás Aquino
 \date 13/02/2024
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
#include "PTOF/Locator_Parallel.h"
#include "PTOF/Store.h"

namespace ptof
{
  /** \struct Geometry_Parallel PTOF/Geometry_Parallel.h "PTOF/Geometry_Parallel.h"
   * \brief Basic geometry class. */
  template
  <std::size_t dim_val,
  BoundaryConditionSet::Type dynamics,
  typename SearchOption = SearchOptions::SecondNeighborPrecheck>
  struct Geometry_Parallel
  {
    static constexpr std::size_t dim{ dim_val };         /**< Spatial dimension. */
    using Mesh = Foam::fvMesh;                           /**< Mesh type. */
    using Locator = Locator_Cell_Parallel<SearchOption>; /**< Mesh locator type. */
    using BoundaryInfo
      = std::conditional_t<
          dynamics
            == BoundaryConditionSet::Type::transport,
          Store_Absorbed,
          Store_Absorbed_Reinjections>;                  /**< What to store upon hitting boundary. */
    
    /** Constructor.
     \brief Use mesh from specified OpenFOAM case time.
     \param directories_of OpenFOAM case directory information.
     \param directories Current case directory information.
     \param num_threads Number of parallel threads.
     */
    Geometry_Parallel
    (DirectoriesOF const& directories_of,
     Directories const& directories,
     std::size_t num_threads)
    : _meshes(make_meshes(directories_of, num_threads))
    , _mesh_searches(make_mesh_searches(directories_of, num_threads))
    , locator{ _mesh_searches }
    {}
    
    /**
      \brief Make a Boundary object for the preset type of dynamics.
      \details
      - 'transport': No parameters, custom boundary does nothing.
      - 'firstpassage': Parameters should hold InitialCondition object.
      - 'custom': Boundary reinjects according to initial condition.
      \param directories Current case directory information.
      \param params_transport Transport parameters.
      \param params_reaction Reaction parameters.
      \param params_solvers Solver parameters.
      \param surface_reaction Surface reaction.
      \param initial_condition Initial condition object for reinjection (for FPT dynamics).
      \return Boundary object.
     */
    template
    <typename TransportParameters,
    typename ReactionParameters,
    typename SolverParameters,
    typename SurfaceReaction,
    typename InitialCondition = useful::Empty>
    auto makeBoundary
    (Directories const& directories,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     SurfaceReaction&& surface_reaction,
     InitialCondition&& initial_condition = {}) const
    {
      auto boundary_conditions = get_boundary_conditions<dynamics>(directories);
      for (auto const& bc : boundary_conditions)
        if (bc.second == "periodic")
          throw std::runtime_error{
            "Boundary conditions of type periodic not supported by Geometry type" };
      if constexpr (dynamics != BoundaryConditionSet::Type::firstpassage)
        for (auto const& bc : boundary_conditions)
          if (bc.second == "custom")
            throw std::runtime_error{
              std::string{ "Boundary conditions of type custom not supported for " }
              + BoundaryConditionSet::name(dynamics)
              + " dynamics" };
      
      if constexpr (dynamics == BoundaryConditionSet::Type::transport)
        return ptof::Boundary_Cases{
          boundary_conditions,
          locator,
          BoundaryInfo{},
          Boundary_DoNothing{},
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(surface_reaction) };
      if constexpr (dynamics == BoundaryConditionSet::Type::firstpassage)
        return ptof::Boundary_Cases{
          boundary_conditions,
          locator,
          BoundaryInfo{},
          Boundary_DoNothing{},
          Boundary_Reinject{ std::forward<InitialCondition>(initial_condition) } };
      throw std::runtime_error{
        std::string{ "Boundary conditions for " }
          + BoundaryConditionSet::name(dynamics)
          + " not supported" };
    }
    
    /**
     \brief Output generic information about object.
     \param output Output stream.
     */
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
        "\treacting_reflecting: Reflection and surface reaction\n"
        "\tabsorbing: Absorbing\n"
        "\tinfo: Information upon crossing\n"
        "\tcustom: ";
      if constexpr (dynamics != BoundaryConditionSet::Type::firstpassage)
        output <<
          "\tcustom: Reinject according to initial condition\n";
      output <<
        "\tempty: No effect (default for unspecified patches in mesh)\n"
        "--------------------------------------------------\n";
    }
    
    /**
     \brief Output information about current object.
     */
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Mesh\n"
        "--------------------------------------------------\n";
      output << "cells: " << mesh().nCells() << "\n"
             << "faces: " << mesh().nFaces() << "\n"
             << "edges: " << mesh().nEdges() << "\n";
      output <<
        "--------------------------------------------------\n";
    }
    
    /** \return Mesh for the current thread. */
    auto const& mesh() const
    { return *_meshes[omp_get_thread_num()]; }
    
  private:
    using MeshSearch = Foam::meshSearch;            /**< Mesh searching tools type. */
    
    std::vector<std::unique_ptr<Mesh>> _meshes;               /**< Mesh copies needed for parallel searches. */
    std::vector<std::unique_ptr<MeshSearch>> _mesh_searches;  /**< Mesh searching tools. */
    
  public:
    Locator locator;                                          /**< Object to locate positions in mesh. */
    
  private:
    /**
     \param directories_of OpenFOAM case directory information.
     \param num_threads Number of parallel threads.
     \return Mesh copies for parallel searches
    */
    std::vector<std::unique_ptr<Mesh>> make_meshes
    (DirectoriesOF const& directories_of, std::size_t num_threads)
    {
      std::vector<std::unique_ptr<Mesh>> meshes;
      meshes.reserve(num_threads);
      for (std::size_t thread = 0; thread < num_threads; ++thread)
        meshes.emplace_back(std::make_unique<Mesh>(Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ }));
      
      return meshes;
    }
    
    /**
     \param directories_of OpenFOAM case directory information.
     \param num_threads Number of parallel threads.
     \return Mesh search tools for parallel searches
    */
    std::vector<std::unique_ptr<MeshSearch>> make_mesh_searches
    (DirectoriesOF const& directories_of, std::size_t num_threads)
    {
      std::vector<std::unique_ptr<MeshSearch>> mesh_searches;
      mesh_searches.reserve(num_threads);
      for (std::size_t thread = 0; thread < num_threads; ++thread)
        mesh_searches.emplace_back(std::make_unique<MeshSearch>(*_meshes[thread]));

      return mesh_searches;
    }
  };
  
  /** \struct Geometry_Periodic_Cartesian_Parallel PTOF/Geometry_Parallel.h "PTOF/Geometry_Parallel.h"
   *  \brief Geometry class for geometry with Cartesian periodic BCs along some dimensions.*/
  template
  <std::size_t dim_val,
  BoundaryConditionSet::Type dynamics,
  typename SearchOption = SearchOptions::SecondNeighborPrecheck>
  struct Geometry_Periodic_Cartesian_Parallel
  {
    static constexpr std::size_t dim{ dim_val };         /**< Spatial dimension. */
    using Mesh = Foam::fvMesh;                           /**< Mesh type. */
    using Locator = Locator_Cell_Parallel<SearchOption>; /**< Mesh locator type. */
    using Boundary_Periodic =
      geometry::Boundary_Periodic_WithOutsideInfo;       /**< Periodic boundary type. */
    using BoundaryInfo = Store_Absorbed;                 /**< What to store upon hitting boundary. */
    
    /** Constructor.
     \brief Use mesh from specified OpenFOAM case time.
     \param directories_of OpenFOAM case directory information.
     \param directories Current case directory information.
     \param num_threads Number of parallel threads.
    */
    Geometry_Periodic_Cartesian_Parallel
    (DirectoriesOF const& directories_of,
     Directories const& directories,
     std::size_t num_threads)
    : _meshes{ make_meshes(directories_of, num_threads) }
    , _mesh_searches{ make_mesh_searches(directories_of, num_threads) }
    , locator{ _mesh_searches }
    , boundary_periodic{ makeboundary_periodic(directories) }
    {}
    
    /**
     \brief Make a Boundary object for Cartesian periodicity.
     \param directories Current case directory information.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
     \param surface_reaction Surface reaction.
     \param initial_condition Initial condition object for reinjection (for FPT dynamics).
     \return Boundary object. */
    template
    <typename TransportParameters,
    typename ReactionParameters,
    typename SolverParameters,
    typename SurfaceReaction,
    typename InitialCondition = useful::Empty>
    auto makeBoundary
    (Directories const& directories,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     SurfaceReaction&& surface_reaction,
     InitialCondition&& initial_condition = {}) const
    {
      auto boundary_conditions = get_boundary_conditions<dynamics>(directories);
      if constexpr (dynamics != BoundaryConditionSet::Type::firstpassage)
        for (auto const& bc : boundary_conditions)
          if (bc.second == "custom")
            throw std::runtime_error{
              std::string{ "Boundary conditions of type custom not supported for " }
              + BoundaryConditionSet::name(dynamics)
              + " dynamics" };
      
      if constexpr (dynamics == BoundaryConditionSet::Type::transport)
        return ptof::Boundary_Cases{
          boundary_conditions,
          locator,
          BoundaryInfo{},
          Boundary_Periodic_OF{
            boundary_periodic,
            locator },
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(surface_reaction) };
      if constexpr (dynamics == BoundaryConditionSet::Type::firstpassage)
        return ptof::Boundary_Cases{
          boundary_conditions,
          locator,
          BoundaryInfo{},
          Boundary_Periodic_OF{
            boundary_periodic,
            locator },
          Boundary_Reinject{ std::forward<InitialCondition>(initial_condition) } };
      throw std::runtime_error{
        std::string{ "Boundary conditions for " }
          + BoundaryConditionSet::name(dynamics)
          + " not supported" };
    }
    
    /** \brief Output generic information about object. */
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
        "\treacting_reflecting: Reflection and surface reaction\n"
        "\tperiodic: Cartesian periodic, position extracted from mesh\n";
      output <<
        "\tabsorbing: Absorbing\n"
        "\tinfo: Information upon crossing\n";
      if constexpr (dynamics != BoundaryConditionSet::Type::firstpassage)
        output <<
          "\tcustom: Reinject according to initial condition\n";
      output <<
        "\tempty: No effect (default for unspecified patches in mesh)\n"
      "--------------------------------------------------\n";
    }
    
    /** \brief Output information about current object. */
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Mesh\n"
        "--------------------------------------------------\n";
      output << "cells: " << mesh().nCells() << "\n"
             << "faces: " << mesh().nFaces() << "\n"
             << "edges: " << mesh().nEdges() << "\n"
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
    
    /** \return Mesh for the current thread. */
    auto const& mesh() const
    { return *_meshes[omp_get_thread_num()]; }
    
  private:
    using MeshSearch = Foam::meshSearch;                      /**< Mesh searching tools type. */
    
    std::vector<std::unique_ptr<Mesh>> _meshes;               /**< Mesh copies needed for parallel searches. */
    std::vector<std::unique_ptr<MeshSearch>> _mesh_searches;  /**< Mesh searching tools. */
    
  public:
    Locator locator;                                          /**< Object to locate positions in mesh. */
    Boundary_Periodic boundary_periodic;                      /**< Boundary object to enforce periodicity. */
    
  private:
    /**
    \brief Make a boundary object to enforce periodicity based on cartesian planes.
    \param directories Current case directory information.
    \return Periodic boundary object.
    */
    Boundary_Periodic makeboundary_periodic
    (Directories const& directories)
    {
      auto boundary_conditions = get_boundary_conditions<dynamics>(directories);
      
      std::vector<std::vector<double>> boundaries_dim(dim);
      for (auto const& bc : boundary_conditions)
        if (bc.second == "periodic")
        {
          // Find average and variance of patch face centers.
          auto const& patch_id =
            mesh().boundaryMesh().findPatchID(bc.first, false);
          auto const& patch = mesh().boundaryMesh()[patch_id];
          auto face_start = patch.start();
          Foam::point average_face_center = Foam::point::zero;
          Foam::point variance_face_center = Foam::point::zero;
          Foam::scalar total_area = 0.;
          for (auto face = face_start;
               face < face_start + patch.size(); ++face)
          {
            auto area = face_area(face, mesh());
            auto weighted_center = area*face_center(face, mesh());
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
          
          // Choose dimension with least variance.
          auto it_min = std::min_element(variance_face_center.begin(),
                                         variance_face_center.end());
          std::size_t dd = std::distance(variance_face_center.begin(),
                                         it_min);
          boundaries_dim[dd].push_back(average_face_center[dd]);
        }
      
      // Check which dimensions have periodic BCs.
      // If a dimension has periodic BCs they must be exactly two.
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
      
      // Periodic conditions are assumed to be checked in order
      // across Cartesian dimensions.
      // If there are periodic BCs in all dimensions
      // starting from 0 up to some dimension, add only those boundaries
      // If there are gap dimensions between periodic dimensions,
      // add {-inf, inf} boundaries.
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
    }
    
    std::vector<std::unique_ptr<Mesh>> make_meshes
    (DirectoriesOF const& directories_of, std::size_t num_threads)
    {
      std::vector<std::unique_ptr<Mesh>> meshes;
      meshes.reserve(num_threads);
      for (std::size_t thread = 0; thread < num_threads; ++thread)
        meshes.emplace_back(std::make_unique<Mesh>(Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ }));
      
      return meshes;
    }
    
    /**
     \param directories_of OpenFOAM case directory information.
     \param num_threads Number of parallel threads.
     \return Mesh search tools for parallel searches
    */
    std::vector<std::unique_ptr<MeshSearch>> make_mesh_searches
    (DirectoriesOF const& directories_of, std::size_t num_threads)
    {
      std::vector<std::unique_ptr<MeshSearch>> mesh_searches;
      mesh_searches.reserve(num_threads);
      for (std::size_t thread = 0; thread < num_threads; ++thread)
        mesh_searches.emplace_back(std::make_unique<MeshSearch>(*_meshes[thread]));

      return mesh_searches;
    }
  };
  
  /** \struct Geometry_Bcc_Parallel PTOF/Geometry_Parallel.h "PTOF/Geometry_Parallel.h"
   \brief Geometry class for periodic representation of a body centered cubic beadpack.
   \details Holds spatial dimension info, mesh, mesh search objects, and periodic boundary to enforce periodicity. */
  template
  <std::size_t dim_val,
  BoundaryConditionSet::Type dynamics,
  BoundaryConditionSet::Type periodicity_type,
  typename SearchOption = SearchOptions::SecondNeighborPrecheck>
  struct Geometry_Bcc_Parallel
  {
    static constexpr std::size_t dim{ dim_val };          /**< Spatial dimension. */
    using Mesh = Foam::fvMesh;                            /**< Mesh type. */
    using Locator = Locator_Cell_Parallel<SearchOption>;  /**< Mesh locator type. */
    using Boundary_Periodic = std::conditional_t<
      periodicity_type
      == BoundaryConditionSet::Type::cartesian,
    geometry::Boundary_Periodic_WithOutsideInfo,
    geometry::Boundary_Periodic_SymmetryPlanes_WithOutsideInfo<
      geometry::SymmetryPlanes_Bcc>>;                     /**< Periodic boundary type. */
    using BoundaryInfo = Store_Absorbed;                  /**< What to store upon hitting boundary. */
    
    /** Constructor.
     \brief Use mesh from specified OpenFOAM case time.
     \brief Use mesh from specified OpenFOAM case time.
     \param directories_of OpenFOAM case directory information.
     \param directories Current case directory information.
     \param num_threads Number of parallel threads.
     */
    Geometry_Bcc_Parallel
    (DirectoriesOF const& directories_of,
     Directories const& directories,
     std::size_t num_threads)
    : _meshes(make_meshes(directories_of, num_threads))
    , _mesh_searches(make_mesh_searches(directories_of, num_threads))
    , locator{ _mesh_searches }
    , radius{
      std::is_same_v<periodicity_type, BoundaryConditionSet::Type::cartesian>
      ? std::pow(sum(mesh().V())/(1.-std::sqrt(3.)*constants::pi/8.), 1./3.)*std::sqrt(3.)/4.
      : 1. }
    , boundary_periodic{ makeboundary_periodic(useful::Selector<
                                               BoundaryConditionSet::Type,
                                               periodicity_type>{},
                                               directories) }
    {}
    
    /**
      \brief Make a Boundary object for the preset type of dynamics and periodicity.
      \details
      - 'cartesian': Cartesian periodicity.
      - 'symmetryplanes': Periodicity according to symmetry planes.
      \param directories Current case directory information.
      \param params_transport Transport parameters.
      \param params_reaction Reaction parameters.
      \param params_solvers Solver parameters.
      \param surface_reaction Surface reaction.
      \param initial_condition Initial condition object for reinjection (for FPT dynamics). */
    template
    <typename TransportParameters,
    typename ReactionParameters,
    typename SolverParameters,
    typename SurfaceReaction,
    typename InitialCondition = useful::Empty>
    auto makeBoundary
    (Directories const& directories,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     SurfaceReaction&& surface_reaction,
     InitialCondition&& initial_condition = {}) const
    {
      auto boundary_conditions = get_boundary_conditions<dynamics>(directories);
      if constexpr (dynamics != BoundaryConditionSet::Type::firstpassage)
        for (auto const& bc : boundary_conditions)
          if (bc.second == "custom")
            throw std::runtime_error{
              std::string{ "Boundary conditions of type custom not supported for " }
              + BoundaryConditionSet::name(dynamics)
              + " dynamics" };
      
      if constexpr (dynamics == BoundaryConditionSet::Type::transport)
        return ptof::Boundary_Cases{
          boundary_conditions,
          locator,
          BoundaryInfo{},
          Boundary_Periodic_OF{
            boundary_periodic,
            locator },
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(surface_reaction) };
      if constexpr (dynamics == BoundaryConditionSet::Type::firstpassage)
        return ptof::Boundary_Cases{
          boundary_conditions,
          locator,
          BoundaryInfo{},
          Boundary_Periodic_OF{
            boundary_periodic,
            locator },
          Boundary_Reinject{ std::forward<InitialCondition>(initial_condition) } };
      throw std::runtime_error{
        std::string{ "Boundary conditions for " }
          + BoundaryConditionSet::name(dynamics)
          + " not supported" };
    }
    
    /**
     \brief Output generic information about object.
     */
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
        "\treacting_reflecting: Reflection and surface reaction\n"
        "\tperiodic: ";
      output <<
        (BoundaryConditionSet::name(periodicity_type) == "cartesian"
         ? "Periodic in the primitive unit cell\n"
         : "Periodic in the minimal unit cell\n");
      output <<
        "\tabsorbing: Absorbing\n"
        "\tinfo: Information upon crossing\n";
      if constexpr (dynamics == BoundaryConditionSet::Type::firstpassage)
        output <<
          "\tcustom: Reinject according to initial condition\n";
      output <<
        "\tempty: No effect (default for unspecified patches in mesh)\n"
      "--------------------------------------------------\n";
    }
    
    /**
     \brief Output information about current object.
     */
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Mesh\n"
        "--------------------------------------------------\n";
      output << "cells: " << mesh().nCells() << "\n"
             << "faces: " << mesh().nFaces() << "\n"
             << "edges: " << mesh().nEdges() << "\n";
      output <<
        "--------------------------------------------------\n";
    }
    
    /** \return Mesh for the current thread. */
    auto const& mesh() const
    { return *_meshes[omp_get_thread_num()]; }
    
  private:
    using MeshSearch = Foam::meshSearch;                      /**< Mesh searching tools type. */
    
    std::vector<std::unique_ptr<Mesh>> _meshes;               /**< Mesh copies needed for parallel searches. */
    std::vector<std::unique_ptr<MeshSearch>> _mesh_searches;  /**< Mesh searching tools. */
    
  public:
    const double radius;                                      /**< Bead radius. */
    Locator locator;                                          /**< Object to locate positions in mesh. */
    Boundary_Periodic boundary_periodic;                      /**< Boundary object to enforce periodicity.*/
  
  private:
    /**
     \brief Make a boundary object to enforce cartesian periodicity.
     \return Periodic boundary object.
     */
    Boundary_Periodic makeboundary_periodic
    (useful::Selector<BoundaryConditionSet::Type,
      BoundaryConditionSet::Type::cartesian>,
      Directories const& directories)
    {
      return {
        extract_cartesian_periodic_boundaries(get_boundary_conditions<dynamics>(directories),
                                              mesh(),
                                              dim) };
    };
  
    /**
    \brief Make a boundary object to enforce periodicity based on symmetry planes.
    \return Periodic boundary object.
    */
    Boundary_Periodic makeboundary_periodic
    (useful::Selector<BoundaryConditionSet::Type,
      BoundaryConditionSet::Type::symmetryplanes>,
     Directories const& directories)
    { return { radius }; }
    
    std::vector<std::unique_ptr<Mesh>> make_meshes
    (DirectoriesOF const& directories_of, std::size_t num_threads)
    {
      std::vector<std::unique_ptr<Mesh>> meshes;
      meshes.reserve(num_threads);
      for (std::size_t thread = 0; thread < num_threads; ++thread)
        meshes.emplace_back(std::make_unique<Mesh>(Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ }));
      
      return meshes;
    }
    
    /**
     \param directories_of OpenFOAM case directory information.
     \param num_threads Number of parallel threads.
     \return Mesh search tools for parallel searches
    */
    std::vector<std::unique_ptr<MeshSearch>> make_mesh_searches
    (DirectoriesOF const& directories_of, std::size_t num_threads)
    {
      std::vector<std::unique_ptr<MeshSearch>> mesh_searches;
      mesh_searches.reserve(num_threads);
      for (std::size_t thread = 0; thread < num_threads; ++thread)
        mesh_searches.emplace_back(std::make_unique<MeshSearch>(*_meshes[thread]));

      return mesh_searches;
    }
  };
}

#endif /* PTOF_GEOMETRY_PARALLEL_H */
