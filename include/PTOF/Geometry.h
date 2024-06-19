/**
 \file PTOF/Geometry.h
 \author Tomás Aquino
 \date 09/03/2022
*/
#ifndef PTOF_GEOMETRY_H
#define PTOF_GEOMETRY_H

#include "General/Constants.h"
#include "General/Meta.h"
#include "General/Serial.h"
#include "General/Useful.h"
#include "Geometry/Boundary.h"
#include "PTOF/Boundary.h"
#include "PTOF/Directories.h"
#include "PTOF/Locator.h"
#include "PTOF/Store.h"
#include "PTOF/Useful.h"
#include <IOobject.H>
#include <cmath>
#include <cstddef>
#include <fvMesh.H>
#include <limits>
#include <meshSearch.H>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ptof {
/** \class Geometry PTOF/Geometry.h "PTOF/Geometry.h"
 * \brief Basic geometry class.
 \note Does not handle periodic boundaries. */
template <std::size_t dim_val, typename ParallelOption_t,
          Dynamics::Type dynamics = Dynamics::Type::transport,
          typename SearchOption = SearchOptions::SecondNeighborPrecheck>
struct Geometry {
  static constexpr std::size_t dim{dim_val}; /**< Spatial dimension. */
  using ParallelOption =
      ParallelOption_t;      /**< Choose serial or parallel implementations. */
  using Mesh = Foam::fvMesh; /**< Mesh. */
  using MeshSearch = Foam::meshSearch; /**< Mesh searching tools. */
  using Locator =
      Locator_Cell<Geometry, SearchOption>; /**< Locate positions in mesh.*/
  using BoundaryInfo =
      std::conditional_t<dynamics == Dynamics::Type::transport, Store_Absorbed,
                         Store_Absorbed_Reinjections>; /**< What to store upon
                                                          hitting a boundary. */

  /** Constructor.
  \brief Use mesh from specified OpenFOAM case time.
  \param directories_of OpenFOAM case directory information.
  \param directories Current case directory information.
  */
  Geometry(DirectoriesOF const &directories_of, Directories const &directories)
      : _meshes(make_meshes(directories_of)),
        _mesh_searches(make_mesh_searches(directories_of)) {}

  /**
   \brief Make a Boundary object for the preset type of dynamics.

   - 'transport': No parameters, custom boundary does nothing.
   - 'firstpassage': Parameters should hold InitialCondition object.
   - 'custom': Boundary reinjects according to initial condition.
   \param directories Current case directory information.
   \param params_transport Transport parameters.
   \param params_reaction Reaction parameters.
   \param params_solvers Solver parameters.
   \param surface_reaction Surface reaction.
   \param initial_condition Initial condition object for reinjection (for FPT
   dynamics). \return Boundary object.
  */
  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters, typename SurfaceReaction,
            typename InitialCondition = useful::Empty>
  auto makeBoundary(Directories const &directories,
                    TransportParameters const &params_transport,
                    ReactionParameters const &params_reaction,
                    SolverParameters const &params_solvers,
                    SurfaceReaction &&surface_reaction,
                    InitialCondition &&initial_condition = {}) const {
    auto boundary_conditions = get_boundary_conditions<dynamics>(directories);
    for (auto const &bc : boundary_conditions)
      if (bc.second == "periodic")
        throw std::runtime_error{"Boundary conditions of type periodic not "
                                 "supported by Geometry type"};
    if constexpr (dynamics != Dynamics::Type::firstpassage)
      for (auto const &bc : boundary_conditions)
        if (bc.second == "custom")
          throw std::runtime_error{
              std::string{
                  "Boundary conditions of type custom not supported for "} +
              Dynamics::name(dynamics) + " dynamics"};

    if constexpr (dynamics == Dynamics::Type::transport)
      return ptof::Boundary_Cases{
          boundary_conditions,
          locator,
          BoundaryInfo{},
          Boundary_DoNothing{},
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(surface_reaction)};
    if constexpr (dynamics == Dynamics::Type::firstpassage)
      return ptof::Boundary_Cases{
          boundary_conditions, locator, BoundaryInfo{}, Boundary_DoNothing{},
          Boundary_Reinject{std::forward<InitialCondition>(initial_condition)}};
    throw std::runtime_error{std::string{"Boundary conditions for "} +
                             Dynamics::name(dynamics) + " not supported"};
  }

  /**
  \brief Output generic information about object.
  \param output Output stream
  */
  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------------------\n"
              "Geometry\n"
              "--------------------------------------------------------------\n"
              "Spatial dimension: " +
                  std::to_string(dim) +
                  "\n"
                  "Boundary conditions for: " +
                  Dynamics::name(dynamics) +
                  "\n"
                  "Boundary condition types:\n"
                  "\treflecting: Reflecting\n"
                  "\treacting_reflecting: Reflection and surface reaction\n"
                  "\tabsorbing: Absorbing\n"
                  "\tinfo: Information upon crossing\n";
    if constexpr (dynamics != Dynamics::Type::firstpassage)
      output << "\tcustom: Reinject according to initial condition\n";
    output
        << "\tempty: No effect (default for unspecified patches in mesh)\n"
           "--------------------------------------------------------------\n";
  }

  /**
     \brief Output information about current object.
   */
  template <typename OStream> void info_runtime(OStream &output) const {
    output
        << "--------------------------------------------------------------\n"
           "Mesh\n"
           "--------------------------------------------------------------\n";
    output << "cells: " << mesh().nCells() << "\n"
           << "faces: " << mesh().nFaces() << "\n"
           << "edges: " << mesh().nEdges() << "\n";
    output
        << "--------------------------------------------------------------\n";
  }

  /** \return Mesh. */
  auto const &mesh() const {
    return *_meshes[useful::get_thread_num(ParallelOption{})];
  }

  /** \return Mesh search tools. */
  auto const &mesh_search() const {
    return *_mesh_searches[useful::get_thread_num(ParallelOption{})];
  }

private:
  std::vector<std::unique_ptr<Mesh>>
      _meshes; /**< Mesh or mesh copies needed for parallel searches. */
  std::vector<std::unique_ptr<MeshSearch>>
      _mesh_searches; /**< Mesh searching tools. */

  /**
     \param directories_of OpenFOAM case directory information.
     \return Mesh copies for parallel searches
    */
  std::vector<std::unique_ptr<Mesh>>
  make_meshes(DirectoriesOF const &directories_of) {
    std::size_t num_threads = useful::get_num_threads(ParallelOption{});
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.reserve(num_threads);
    for (std::size_t thread = 0; thread < num_threads; ++thread)
      meshes.emplace_back(std::make_unique<Mesh>(Foam::IOobject{
          Foam::fvMesh::defaultRegion, directories_of.time.timeName(),
          directories_of.time, Foam::IOobject::MUST_READ}));
    return meshes;
  }

  /**
   \param directories_of OpenFOAM case directory information.
   \return Mesh search tools for parallel searches
  */
  std::vector<std::unique_ptr<MeshSearch>>
  make_mesh_searches(DirectoriesOF const &directories_of) {
    std::size_t num_threads = useful::get_num_threads(ParallelOption{});
    std::vector<std::unique_ptr<MeshSearch>> mesh_searches;
    mesh_searches.reserve(num_threads);
    for (std::size_t thread = 0; thread < num_threads; ++thread)
      mesh_searches.emplace_back(
          std::make_unique<MeshSearch>(*_meshes[thread]));

    return mesh_searches;
  }

public:
  Locator locator{*this}; /**< To locate positions in mesh.*/
};

/** \class Geometry_Periodic_Cartesian PTOF/Geometry.h "PTOF/Geometry.h"
 *  \brief Geometry class for geometry to handle periodic BCs along any number
 * of Cartesian dimensions.*/
template <std::size_t dim_val, typename ParallelOption_t,
          Dynamics::Type dynamics = Dynamics::Type::transport,
          typename SearchOption = SearchOptions::SecondNeighborPrecheck>
struct Geometry_Periodic_Cartesian {
  static constexpr std::size_t dim{dim_val}; /**< Spatial dimension. */
  using ParallelOption =
      ParallelOption_t;      /**< Choose serial or parallel implementations. */
  using Mesh = Foam::fvMesh; /**< Mesh. */
  using MeshSearch = Foam::meshSearch; /**< Mesh searching tools. */
  using Locator = Locator_Cell<Geometry_Periodic_Cartesian,
                               SearchOption>; /**< Locate positions in mesh.*/
  using Boundary_Periodic =
      geom::Boundary_Periodic_WithOutsideInfo; /**< Periodic boundary. */
  using BoundaryInfo =
      Store_Absorbed; /**< What to store upon hitting a boundary. */

  /** Constructor.
   \brief Use mesh from specified OpenFOAM case time.
   \param directories_of OpenFOAM case directory information.
   \param directories Current case directory information.
  */
  Geometry_Periodic_Cartesian(DirectoriesOF const &directories_of,
                              Directories const &directories)
      : _meshes(make_meshes(directories_of)),
        _mesh_searches(make_mesh_searches(directories_of)),
        boundary_periodic{extract_cartesian_periodic_boundaries(
            get_boundary_conditions<dynamics>(directories), mesh(), dim)} {}

  /**
   \brief Make a Boundary object for Cartesian periodicity.
   \param directories Current case directory information.
   \param params_transport Transport parameters.
   \param params_reaction Reaction parameters.
   \param params_solvers Solver parameters.
   \param surface_reaction Surface reaction.
   \param initial_condition Initial condition object for reinjection (for FPT
   dynamics). \return Boundary object.
  */
  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters, typename SurfaceReaction,
            typename InitialCondition = useful::Empty>
  auto makeBoundary(Directories const &directories,
                    TransportParameters const &params_transport,
                    ReactionParameters const &params_reaction,
                    SolverParameters const &params_solvers,
                    SurfaceReaction &&surface_reaction,
                    InitialCondition &&initial_condition = {}) const {
    auto boundary_conditions = get_boundary_conditions<dynamics>(directories);
    if constexpr (dynamics != Dynamics::Type::firstpassage)
      for (auto const &bc : boundary_conditions)
        if (bc.second == "custom")
          throw std::runtime_error{
              std::string{
                  "Boundary conditions of type custom not supported for "} +
              Dynamics::name(dynamics) + " dynamics"};

    if constexpr (dynamics == Dynamics::Type::transport)
      return ptof::Boundary_Cases{
          boundary_conditions,
          locator,
          BoundaryInfo{},
          Boundary_Periodic_OF{boundary_periodic, locator},
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(surface_reaction)};
    if constexpr (dynamics == Dynamics::Type::firstpassage)
      return ptof::Boundary_Cases{
          boundary_conditions, locator, BoundaryInfo{},
          Boundary_Periodic_OF{boundary_periodic, locator},
          Boundary_Reinject{std::forward<InitialCondition>(initial_condition)}};
    throw std::runtime_error{std::string{"Boundary conditions for "} +
                             Dynamics::name(dynamics) + " not supported"};
  }

  /** \brief Output generic information about object. */
  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------------------\n"
              "Geometry\n"
              "--------------------------------------------------------------\n"
              "Spatial dimension: " +
                  std::to_string(dim) +
                  "\n"
                  "Boundary conditions for: " +
                  Dynamics::name(dynamics) +
                  "\n"
                  "Boundary condition types:\n"
                  "\treflecting: Reflecting\n"
                  "\treacting_reflecting: Reflection and surface reaction\n"
                  "\tperiodic: Cartesian periodic, position extracted from "
                  "mesh\n";
    output << "\tabsorbing: Absorbing\n"
              "\tinfo: Information upon crossing\n";
    if constexpr (dynamics != Dynamics::Type::firstpassage)
      output << "\tcustom: Reinject according to initial condition\n";
    output
        << "\tempty: No effect (default for unspecified patches in mesh)\n"
           "--------------------------------------------------------------\n";
  }

  /** \brief Output information about current object. */
  template <typename OStream> void info_runtime(OStream &output) const {
    output << "--------------------------------------------------------------\n"
              "Mesh\n"
              "--------------------------------------------------------------\n"
           << "cells: " << mesh().nCells() << "\n"
           << "faces: " << mesh().nFaces() << "\n"
           << "edges: " << mesh().nEdges()
           << "\n"
              "--------------------------------------------------------------\n"
           << "--------------------------------------------------------------\n"
              "Periodic boundaries\n"
              "--------------------------------------------------------------\n"
           << "Automatically extracted Cartesian periodic boundaries:\n";
    if (boundary_periodic.boundaries.size() == 0)
      output << "None";
    for (std::size_t dd = 0; dd < boundary_periodic.boundaries.size(); ++dd)
      if (boundary_periodic.boundaries[dd].second <
          std::numeric_limits<double>::infinity())
        output << "Dimension " << dd << " at "
               << boundary_periodic.boundaries[dd].first << " and "
               << boundary_periodic.boundaries[dd].second << "\n";
    output
        << "--------------------------------------------------------------\n";
  }

  /** \return Mesh. */
  auto const &mesh() const {
    return *_meshes[useful::get_thread_num(ParallelOption{})];
  }

  /** \return Mesh search tools. */
  auto const &mesh_search() const {
    return *_mesh_searches[useful::get_thread_num(ParallelOption{})];
  }

private:
  std::vector<std::unique_ptr<Mesh>>
      _meshes; /**< Mesh or mesh copies needed for parallel searches. */
  std::vector<std::unique_ptr<MeshSearch>>
      _mesh_searches; /**< Mesh searching tools. */

  /**
   \param directories_of OpenFOAM case directory information.
   \return Mesh copies for parallel searches
  */
  std::vector<std::unique_ptr<Mesh>>
  make_meshes(DirectoriesOF const &directories_of) {
    std::size_t num_threads = useful::get_num_threads(ParallelOption{});
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.reserve(num_threads);
    for (std::size_t thread = 0; thread < num_threads; ++thread)
      meshes.emplace_back(std::make_unique<Mesh>(Foam::IOobject{
          Foam::fvMesh::defaultRegion, directories_of.time.timeName(),
          directories_of.time, Foam::IOobject::MUST_READ}));
    return meshes;
  }

  /**
   \param directories_of OpenFOAM case directory information.
   \return Mesh search tools for parallel searches
  */
  std::vector<std::unique_ptr<MeshSearch>>
  make_mesh_searches(DirectoriesOF const &directories_of) {
    std::size_t num_threads = useful::get_num_threads(ParallelOption{});
    std::vector<std::unique_ptr<MeshSearch>> mesh_searches;
    mesh_searches.reserve(num_threads);
    for (std::size_t thread = 0; thread < num_threads; ++thread)
      mesh_searches.emplace_back(
          std::make_unique<MeshSearch>(*_meshes[thread]));

    return mesh_searches;
  }

public:
  Locator locator{*this}; /**< To locate positions in mesh.*/
  Boundary_Periodic
      boundary_periodic; /**< Boundary object to enforce periodicity.*/
};

/** \class Geometry_Bcc PTOF/Geometry.h "PTOF/Geometry.h"
 \brief Geometry class for periodic representation of a body centered cubic
 beadpack. \details Holds spatial dimension info, mesh, mesh search objects,
 and periodic boundary to enforce periodicity. */
template <typename ParallelOption_t,
          Periodicity::Type periodicity = Periodicity::Type::cartesian,
          Dynamics::Type dynamics = Dynamics::Type::transport,
          typename SearchOption = SearchOptions::SecondNeighborPrecheck>
struct Geometry_Bcc {
  static constexpr std::size_t dim{3}; /**< Spatial dimension. */
  using ParallelOption =
      ParallelOption_t;      /**< Choose serial or parallel implementations. */
  using Mesh = Foam::fvMesh; /**< Mesh. */
  using MeshSearch = Foam::meshSearch; /**< Mesh searching tools. */
  using Locator =
      Locator_Cell<Geometry_Bcc, SearchOption>; /**< Locate positions in mesh.*/
  using Boundary_Periodic =
      std::conditional_t<periodicity == Periodicity::Type::cartesian,
                         geom::Boundary_Periodic_WithOutsideInfo,
                         geom::Boundary_Periodic_SymmetryPlanes_WithOutsideInfo<
                             geom::SymmetryPlanes_Bcc>>; /**< Periodic boundary
                                                            type. */
  using BoundaryInfo =
      Store_Absorbed; /**< What to store upon hitting a boundary. */

  /** Constructor.
  \brief Use mesh from specified OpenFOAM case time.
  \param directories_of OpenFOAM case directory information.
  \param directories Current case directory information.
  */
  Geometry_Bcc(DirectoriesOF const &directories_of,
               Directories const &directories)
      : _meshes(make_meshes(directories_of)),
        _mesh_searches(make_mesh_searches(directories_of)),
        radius{periodicity == Periodicity::Type::cartesian
                   ? std::pow(sum(mesh().cellVolumes()) /
                                  (1. - std::sqrt(3.) * constants::pi / 8.),
                              1. / 3.) *
                         std::sqrt(3.) / 4.
                   : 1.},
        boundary_periodic{make_boundary_periodic(
            meta::Selector<Periodicity::Type, periodicity>{}, directories)} {}

  /**
   \brief Make a Boundary object for the preset type of dynamics and
   periodicity \details
   - 'cartesian': Cartesian periodicity
   - 'symmetryplanes': Periodicity according to symmetry planes
   \param directories Current case directory information.
   \param params_transport Transport parameters.
   \param params_reaction Reaction parameters.
   \param params_solvers Solver parameters.
   \param surface_reaction Surface reaction.
   \param initial_condition Initial condition object for reinjection (for FPT
   dynamics). */
  template <typename TransportParameters, typename ReactionParameters,
            typename SolverParameters, typename SurfaceReaction,
            typename InitialCondition = useful::Empty>
  auto makeBoundary(Directories const &directories,
                    TransportParameters const &params_transport,
                    ReactionParameters const &params_reaction,
                    SolverParameters const &params_solvers,
                    SurfaceReaction &&surface_reaction,
                    InitialCondition &&initial_condition = {}) const {
    auto boundary_conditions = get_boundary_conditions<dynamics>(directories);
    if constexpr (dynamics != Dynamics::Type::firstpassage)
      for (auto const &bc : boundary_conditions)
        if (bc.second == "custom")
          throw std::runtime_error{
              std::string{
                  "Boundary conditions of type custom not supported for "} +
              Dynamics::name(dynamics) + " dynamics"};

    if constexpr (dynamics == Dynamics::Type::transport)
      return ptof::Boundary_Cases{
          boundary_conditions,
          locator,
          BoundaryInfo{},
          Boundary_Periodic_OF{boundary_periodic, locator},
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(surface_reaction)};
    if constexpr (dynamics == Dynamics::Type::firstpassage)
      return ptof::Boundary_Cases{
          boundary_conditions, locator, BoundaryInfo{},
          Boundary_Periodic_OF{boundary_periodic, locator},
          Boundary_Reinject{std::forward<InitialCondition>(initial_condition)}};
    throw std::runtime_error{std::string{"Boundary conditions for "} +
                             Dynamics::name(dynamics) + " not supported"};
  }

  /**
  \brief Output generic information about object.
  */
  template <typename OStream> static void info(OStream &output) {
    output << "--------------------------------------------------------------\n"
              "Geometry\n"
              "--------------------------------------------------------------\n"
              "Spatial dimension: " +
                  std::to_string(dim) +
                  "\n"
                  "Boundary conditions for: " +
                  Dynamics::name(dynamics) +
                  "\n"
                  "Boundary condition types:\n"
                  "\treflecting: Reflecting\n"
                  "\treacting_reflecting: Reflection and surface reaction\n"
                  "\tperiodic: "
           << (Periodicity::name(periodicity) == "cartesian"
                   ? "Periodic in the primitive unit cell\n"
                   : "Periodic in the minimal unit cell\n")
           << "\tabsorbing: Absorbing\n"
              "\tinfo: Information upon crossing\n";
    if constexpr (dynamics == Dynamics::Type::firstpassage)
      output << "\tcustom: Reinject according to initial condition\n";
    output
        << "\tempty: No effect (default for unspecified patches in mesh)\n"
           "--------------------------------------------------------------\n";
  }

  /**
  \brief Output information about current object.
  */
  template <typename OStream> void info_runtime(OStream &output) const {
    output
        << "--------------------------------------------------------------\n"
           "Mesh\n"
           "--------------------------------------------------------------\n"
        << "cells: " << mesh().nCells() << "\n"
        << "faces: " << mesh().nFaces() << "\n"
        << "edges: " << mesh().nEdges() << "\n"
        << "--------------------------------------------------------------\n";
  }

  /** \return Mesh. */
  auto const &mesh() const {
    return *_meshes[useful::get_thread_num(ParallelOption{})];
  }

  /** \return Mesh search tools. */
  auto const &mesh_search() const {
    return *_mesh_searches[useful::get_thread_num(ParallelOption{})];
  }

private:
  std::vector<std::unique_ptr<Mesh>>
      _meshes; /**< Mesh or mesh copies needed for parallel searches. */
  std::vector<std::unique_ptr<MeshSearch>>
      _mesh_searches; /**< Mesh searching tools. */

  /**
   \param directories_of OpenFOAM case directory information.
   \return Mesh copies for parallel searches
  */
  std::vector<std::unique_ptr<Mesh>>
  make_meshes(DirectoriesOF const &directories_of) {
    std::size_t num_threads = useful::get_num_threads(ParallelOption{});
    std::vector<std::unique_ptr<Mesh>> meshes;
    meshes.reserve(num_threads);
    for (std::size_t thread = 0; thread < num_threads; ++thread)
      meshes.emplace_back(std::make_unique<Mesh>(Foam::IOobject{
          Foam::fvMesh::defaultRegion, directories_of.time.timeName(),
          directories_of.time, Foam::IOobject::MUST_READ}));
    return meshes;
  }

  /**
   \param directories_of OpenFOAM case directory information.
   \return Mesh search tools for parallel searches
  */
  std::vector<std::unique_ptr<MeshSearch>>
  make_mesh_searches(DirectoriesOF const &directories_of) {
    std::size_t num_threads = useful::get_num_threads(ParallelOption{});
    std::vector<std::unique_ptr<MeshSearch>> mesh_searches;
    mesh_searches.reserve(num_threads);
    for (std::size_t thread = 0; thread < num_threads; ++thread)
      mesh_searches.emplace_back(
          std::make_unique<MeshSearch>(*_meshes[thread]));

    return mesh_searches;
  }

public:
  Locator locator{*this}; /**< To locate positions in mesh.*/
  const double radius;    /**< Bead radius. */
  Boundary_Periodic
      boundary_periodic; /**< Boundary object to enforce periodicity. */

private:
  /**
  \brief Make a boundary object to enforce cartesian periodicity.
  \return Periodic boundary object.
  */
  Boundary_Periodic make_boundary_periodic(
      meta::Selector<Periodicity::Type, Periodicity::Type::cartesian>,
      Directories const &directories) {
    return {extract_cartesian_periodic_boundaries(
        get_boundary_conditions<dynamics>(directories), mesh(), dim)};
  };

  /**
  \brief Make a boundary object to enforce periodicity based on symmetry
  planes. \return Periodic boundary object.
  */
  Boundary_Periodic make_boundary_periodic(
      meta::Selector<Periodicity::Type, Periodicity::Type::symmetryplanes>,
      Directories const &directories) {
    return {radius};
  }
};
} // namespace ptof

#endif /* PTOF_GEOMETRY_H */
