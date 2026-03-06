/**
   \file PTOF/Geometry.h
   \author Tomas Aquino
   \date 09/03/2022
   \brief Objects to handle domain geometry and boundaries.
*/
#ifndef PTOF_GEOMETRY_H
#define PTOF_GEOMETRY_H

#include "General/Constant.h"
#include "General/Meta.h"
#include "General/Parallel.h"
#include "Geometry/Boundary.h"
#include "PTOF/Boundary_Cases.h"
#include "PTOF/Directories.h"
#include "PTOF/DynamicsList.h"
#include "PTOF/Locator.h"
#include "PTOF/PeriodicityList.h"
#include "PTOF/SearchOptions.h"
#include "PTOF/Useful.h"
#include <IOobject.H>
#include <cmath>
#include <cstddef>
#include <fvMesh.H>
#include <iostream>
#include <limits>
#include <meshSearch.H>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ptof {
/**
   \class Geometry_Generic PTOF/Geometry.h "PTOF/Geometry.h"
   \brief Generic geometry class.
   \note Does not handle periodic boundaries.
*/
template <std::size_t dim_val,
          typename ParallelOption_t = par::ParallelOptions::Serial,
          DynamicsList::Type dynamics = DynamicsList::Type::transport,
          typename SearchOption = SearchOptions::SecondNeighborPrecheck>
struct Geometry_Generic {
  static constexpr std::size_t dim{dim_val}; /**< Spatial dimension. */
  using ParallelOption =
      ParallelOption_t;      /**< Choose serial or parallel implementations. */
  using Mesh = Foam::fvMesh; /**< Mesh. */
  using MeshSearch = Foam::meshSearch; /**< Mesh searching tools. */
  using Locator = Locator_Cell<Geometry_Generic,
                               SearchOption>; /**< Locate positions in mesh.*/

  /**
     \brief Constructor.
     \details Use mesh from specified OpenFOAM case time.
     \param directories_of OpenFOAM case directory information.
     \param directories Current case directory information.
     \param logger Logger object for optional construction information.
  */
  Geometry_Generic(DirectoriesOF const &directories_of,
                   Directories const &directories,
                   io::Logger &&logger = io::NullLogger{})
      : _mesh(make_mesh(directories_of)), _mesh_search{_mesh} {
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Parallel>) {
      logger("Precomputing demand-driven mesh data...");
      compute_demand_driven_mesh_data(_mesh);
      logger("Done!");
      logger("Precomputing demand-driven mesh search data...");
      compute_demand_driven_meshSearch_data(_mesh_search);
      logger("Done!");
    }
  }

  /**
     \brief Make a Boundary object for the \tparam dynamics type.
     \details
     - transport: No parameters, custom boundary is not supported.
     - firstpassage: Parameters should hold InitialCondition object, custom
     boundary reinjects particles.
     \note Periodic boundaries are not supported.
     \param directories Current case directory information.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
     \param velocity_field Velocity field as a function of state.
     \param surface_reaction Surface reaction.
     \param initial_condition Initial condition object for reinjection (for FPT
     dynamics).
     \return Boundary object.
  */
  template <typename State, typename TransportParameters,
            typename ReactionParameters, typename SolverParameters,
            typename VelocityField, typename SurfaceReaction,
            typename InitialCondition = meta::Empty>
  auto makeBoundary(Directories const &directories,
                    TransportParameters const &params_transport,
                    ReactionParameters const &params_reaction,
                    SolverParameters const &params_solvers,
                    VelocityField &&velocity_field,
                    SurfaceReaction &&surface_reaction,
                    InitialCondition &&initial_condition = {}) const {
    auto boundary_conditions = get_boundary_conditions<dynamics>(directories);
    for (auto const &bc : boundary_conditions) {
      if (bc.second == "periodic") {
        throw std::runtime_error{"Boundary conditions of type periodic not "
                                 "supported by Geometry type"};
      }
    }

    if constexpr (dynamics != DynamicsList::Type::firstpassage) {
      for (auto const &bc : boundary_conditions) {
        if (bc.second == "custom") {
          throw std::runtime_error{std::string{"Boundary conditions for "} +
                                   DynamicsList::name(dynamics) + " : " +
                                   "Boundary condition type custom" + " : " +
                                   "Not supported"};
        }
      }
    }

    if constexpr (dynamics == DynamicsList::Type::transport) {
      return ptof::Boundary_Cases{
          meta::Selector_t<State>{},
          boundary_conditions,
          locator,
          std::forward<VelocityField>(velocity_field),
          Boundary_DoNothing{},
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(surface_reaction)};
    }
    if constexpr (dynamics == DynamicsList::Type::firstpassage) {
      return ptof::Boundary_Cases{
          meta::Selector_t<State>{},
          boundary_conditions,
          locator,
          std::forward<VelocityField>(velocity_field),
          Boundary_DoNothing{},
          Boundary_Reinject{std::forward<InitialCondition>(initial_condition)}};
    }
    throw std::runtime_error{std::string{"Boundary conditions for "} +
                             DynamicsList::name(dynamics) + " not supported"};
  }

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  static std::ostream &info(std::ostream &output) {
    output << io::line() << "Geometry\n"
           << io::line()
           << "Spatial dimension: " + std::to_string(dim) +
                  R"(
Boundary conditions for: )" +
                  DynamicsList::name(dynamics) +
                  R"(
Boundary condition types:
  - reflecting
    - Reflecting.
  - reacting
    - Reflection and surface reaction.
  - absorbing
    - Absorbing.
  - inlet
    - Absorbing for outward flow velocities, reflecting otherwise.)";
    if constexpr (dynamics != DynamicsList::Type::firstpassage) {
      output << R"(
  - custom
    - Reinject according to initial condition.)";
    }
    output << R"(
  - empty
    - No effect (default for unspecified patches in mesh).
)" << io::line();
    return output;
  }

  /** \brief Output information about current object. */
  std::ostream &info_runtime(std::ostream &output) const {
    output << io::line() << "Mesh\n"
           << io::line() << "Cells: " << mesh().nCells() << "\n"
           << "Faces: " << mesh().nFaces() << "\n"
           << "Edges: " << mesh().nEdges() << "\n"
           << io::line();
    return output;
  }

  /** \return Mesh. */
  auto const &mesh() const { return _mesh; }

  /** \return Mesh search tools. */
  auto const &mesh_search() const { return _mesh_search; }

private:
  Mesh _mesh;              /**< Mesh object. */
  MeshSearch _mesh_search; /**< Mesh searching tools. */

  /**
     \param directories_of OpenFOAM case directory information.
     \return Mesh.
  */
  Mesh make_mesh(DirectoriesOF const &directories_of) const {
    return Mesh{Foam::IOobject{Foam::fvMesh::defaultRegion,
                               directories_of.time.timeName(),
                               directories_of.time, Foam::IOobject::MUST_READ}};
  }

public:
  Locator locator{*this}; /**< To locate positions in mesh.*/
};

/**
   \class Geometry_Periodic_Cartesian PTOF/Geometry.h "PTOF/Geometry.h"
   \brief Geometry class for geometry to handle periodic BCs along any number
   of Cartesian dimensions.
*/
template <std::size_t dim_val, typename ParallelOption_t,
          DynamicsList::Type dynamics = DynamicsList::Type::transport,
          typename SearchOption = SearchOptions::SecondNeighborPrecheck>
struct Geometry_Periodic_Cartesian {
  static constexpr std::size_t dim{dim_val}; /**< Spatial dimension. */
  using ParallelOption =
      ParallelOption_t;      /**< Choose serial or parallel implementations. */
  using Mesh = Foam::fvMesh; /**< Mesh. */
  using MeshSearch = Foam::meshSearch; /**< Mesh searching tools. */
  using Locator = Locator_Cell<Geometry_Periodic_Cartesian,
                               SearchOption>; /**< Locate positions in mesh.*/
  using BoundaryPeriodic =
      geom::Boundary_Periodic_WithOutsideInfo; /**< Periodic boundary. */

  /**
     \brief Constructor.
     \details Use mesh from specified OpenFOAM case time.
     \param directories_of OpenFOAM case directory information.
     \param directories Current case directory information.
     \param logger Logger object for optional construction information.
  */
  Geometry_Periodic_Cartesian(DirectoriesOF const &directories_of,
                              Directories const &directories,
                              io::Logger &&logger = io::NullLogger{})
      : _mesh(make_mesh(directories_of)),
        _mesh_search{_mesh},
        boundary_periodic{extract_cartesian_periodic_boundaries(
            get_boundary_conditions<dynamics>(directories), mesh(), dim)} {
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Parallel>) {
      logger("Precomputing demand-driven mesh data...");
      compute_demand_driven_mesh_data(_mesh);
      logger("Done!");
      logger("Precomputing demand-driven mesh search data...");
      compute_demand_driven_meshSearch_data(_mesh_search);
      logger("Done!");
    }
  }

  /**
     \brief Make a Boundary object for the \tparam dynamics type with Cartesian
     periodicity along periodic boundary pairs.
     \details
     - transport: No parameters, custom boundary is not supported.
     - firstpassage: Parameters should hold InitialCondition object, custom
     boundary reinjects particles.
     \param directories Current case directory information.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
     \param velocity_field Velocity field as a function of state.
     \param surface_reaction Surface reaction.
     \param initial_condition Initial condition object for reinjection (for FPT
     dynamics).
     \return Boundary object.
  */
  template <typename State, typename TransportParameters,
            typename ReactionParameters, typename SolverParameters,
            typename VelocityField, typename SurfaceReaction,
            typename InitialCondition = meta::Empty>
  auto makeBoundary(Directories const &directories,
                    TransportParameters const &params_transport,
                    ReactionParameters const &params_reaction,
                    SolverParameters const &params_solvers,
                    VelocityField &&velocity_field,
                    SurfaceReaction &&surface_reaction,
                    InitialCondition &&initial_condition = {}) const {
    auto boundary_conditions = get_boundary_conditions<dynamics>(directories);
    if constexpr (dynamics != DynamicsList::Type::firstpassage) {
      for (auto const &bc : boundary_conditions) {
        if (bc.second == "custom") {
          throw std::runtime_error{std::string{"Boundary conditions for "} +
                                   DynamicsList::name(dynamics) + " : " +
                                   "Boundary condition type custom" + " : " +
                                   "Not supported"};
        }
      }
    }

    if constexpr (dynamics == DynamicsList::Type::transport) {
      return ptof::Boundary_Cases{
          meta::Selector_t<State>{},
          boundary_conditions,
          locator,
          std::forward<VelocityField>(velocity_field),
          Boundary_Periodic{boundary_periodic, locator},
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(surface_reaction)};
    }
    if constexpr (dynamics == DynamicsList::Type::firstpassage) {
      return ptof::Boundary_Cases{
          meta::Selector_t<State>{},
          boundary_conditions,
          locator,
          std::forward<VelocityField>(velocity_field),
          Boundary_Periodic{boundary_periodic, locator},
          Boundary_Reinject{std::forward<InitialCondition>(initial_condition)}};
    }
    throw std::runtime_error{std::string{"Boundary conditions for "} +
                             DynamicsList::name(dynamics) + " not supported"};
  }

  /** \brief Output generic information about object. */
  static std::ostream &info(std::ostream &output) {
    output << io::line() << "Geometry\n"
           << io::line()
           << "Spatial dimension: " + std::to_string(dim) +
                  R"(
Boundary conditions for:)" +
                  DynamicsList::name(dynamics) +
                  R"(
Boundary condition types:
  - reflecting
    - Reflecting.
  - reacting
    - Reflection and surface reaction.
  - periodic
    - Cartesian periodic (pairs extracted from mesh geometry).
  - absorbing
    - Absorbing.)";
    if constexpr (dynamics == DynamicsList::Type::firstpassage) {
      output << R"(
  - custom
    - Reinject according to initial condition.
)";
    }
    output << R"(
  - empty
    - No effect (default for unspecified patches in mesh).
)" << io::line();
    return output;
  }

  /** \brief Output information about current object. */
  std::ostream &info_runtime(std::ostream &output) const {
    output << io::line() << "Mesh\n"
           << io::line() << "Cells: " << mesh().nCells() << "\n"
           << "Faces: " << mesh().nFaces() << "\n"
           << "Edges: " << mesh().nEdges() << "\n"
           << io::line();

    output << io::line()
           << "Automatically extracted Cartesian periodic boundaries\n"
           << io::line();
    if (boundary_periodic.boundaries.size() == 0) {
      output << "None";
    }
    for (std::size_t dd = 0; dd < boundary_periodic.boundaries.size(); ++dd) {
      if (boundary_periodic.boundaries[dd].second <
          std::numeric_limits<double>::infinity()) {
        output << "Dimension " << dd << " at "
               << boundary_periodic.boundaries[dd].first << " and "
               << boundary_periodic.boundaries[dd].second << "\n";
      }
    }
    output << io::line();

    return output;
  }

  /** \return Mesh. */
  auto const &mesh() const { return _mesh; }

  /** \return Mesh search tools. */
  auto const &mesh_search() const { return _mesh_search; }

private:
  Mesh _mesh;              /**< Mesh object. */
  MeshSearch _mesh_search; /**< Mesh searching tools. */

  /**
     \param directories_of OpenFOAM case directory information.
     \return Mesh.
  */
  Mesh make_mesh(DirectoriesOF const &directories_of) {
    return Mesh{Foam::IOobject{Foam::fvMesh::defaultRegion,
                               directories_of.time.timeName(),
                               directories_of.time, Foam::IOobject::MUST_READ}};
  }

public:
  Locator locator{*this}; /**< To locate positions in mesh.*/
  BoundaryPeriodic
      boundary_periodic; /**< Boundary object to enforce periodicity.*/
};

/**
   \class Geometry_Bcc PTOF/Geometry.h "PTOF/Geometry.h"
   \brief Geometry class for periodic representation of a body centered cubic
   beadpack.
   \details Holds spatial dimension info, mesh, mesh search objects, and
   periodic boundary to enforce periodicity.
*/
template <typename ParallelOption_t,
          PeriodicityList::Type periodicity = PeriodicityList::Type::cartesian,
          DynamicsList::Type dynamics = DynamicsList::Type::transport,
          typename SearchOption = SearchOptions::SecondNeighborPrecheck>
struct Geometry_Bcc {
  static constexpr std::size_t dim{3}; /**< Spatial dimension. */
  using ParallelOption =
      ParallelOption_t;      /**< Choose serial or parallel implementations. */
  using Mesh = Foam::fvMesh; /**< Mesh. */
  using MeshSearch = Foam::meshSearch; /**< Mesh searching tools. */
  using Locator =
      Locator_Cell<Geometry_Bcc, SearchOption>; /**< Locate positions in mesh.*/
  using BoundaryPeriodic =
      std::conditional_t<periodicity == PeriodicityList::Type::cartesian,
                         geom::Boundary_Periodic_WithOutsideInfo,
                         geom::Boundary_Periodic_SymmetryPlanes_WithOutsideInfo<
                             geom::SymmetryPlanes_Bcc>>; /**< Periodic boundary
                                                            type. */

  /**
     \brief Constructor.
     \details Use mesh from specified OpenFOAM case time.
     \param directories_of OpenFOAM case directory information.
     \param directories Current case directory information.
     \param logger Logger object for optional construction information.
  */
  Geometry_Bcc(DirectoriesOF const &directories_of,
               Directories const &directories,
               io::Logger &&logger = io::NullLogger{})
      : _mesh(make_mesh(directories_of)), _mesh_search{_mesh},
        radius{periodicity == PeriodicityList::Type::cartesian
                   ? std::pow(sum(mesh().cellVolumes()) /
                                  (1. - std::sqrt(3.) * cnst::pi / 8.),
                              1. / 3.) *
                         std::sqrt(3.) / 4.
                   : 1.},
        boundary_periodic{make_boundary_periodic(
            meta::Selector<PeriodicityList::Type, periodicity>{},
            directories)} {
    if constexpr (std::is_same_v<ParallelOption,
                                 par::ParallelOptions::Parallel>) {
      logger("Precomputing demand-driven mesh data...");
      compute_demand_driven_mesh_data(_mesh);
      logger("Done!");
      logger("Precomputing demand-driven mesh search data...");
      compute_demand_driven_meshSearch_data(_mesh_search);
      logger("Done!");
    }
  }

  /**
     \brief Make a Boundary object for the \tparam dynamics and \tparam
     periodicity types.
     \details
     - Dynamics types:
         - transport: No parameters, custom boundary does nothing.
         - firstpassage: Parameters should hold InitialCondition object,
         custom boundary reinjects particles.
     - Periodicity types:
         - cartesian: Cartesian periodicity.
         - symmetryplanes: Periodicity according to symmetry planes.
     \param directories Current case directory information.
     \param params_transport Transport parameters.
     \param params_reaction Reaction parameters.
     \param params_solvers Solver parameters.
     \param velocity_field Velocity field as a function of state.
     \param surface_reaction Surface reaction.
     \param initial_condition Initial condition object for reinjection (for
     FPT dynamics).
  */
  template <typename State, typename TransportParameters,
            typename ReactionParameters, typename SolverParameters,
            typename VelocityField, typename SurfaceReaction,
            typename InitialCondition = meta::Empty>
  auto makeBoundary(Directories const &directories,
                    TransportParameters const &params_transport,
                    ReactionParameters const &params_reaction,
                    SolverParameters const &params_solvers,
                    VelocityField &&velocity_field,
                    SurfaceReaction &&surface_reaction,
                    InitialCondition &&initial_condition = {}) const {
    auto boundary_conditions = get_boundary_conditions<dynamics>(directories);
    if constexpr (dynamics != DynamicsList::Type::firstpassage) {
      for (auto const &bc : boundary_conditions) {
        if (bc.second == "custom") {
          throw std::runtime_error{
              std::string{"Boundary conditions for dynamics type "} +
              DynamicsList::name(dynamics) + " : " +
              "Boundary condition type custom" + " : " + "Not supported"};
        }
      }
    }

    if constexpr (dynamics == DynamicsList::Type::transport) {
      return ptof::Boundary_Cases{
          meta::Selector_t<State>{},
          boundary_conditions,
          locator,
          std::forward<VelocityField>(velocity_field),
          Boundary_Periodic{boundary_periodic, locator},
          Boundary_DoNothing{},
          std::forward<SurfaceReaction>(surface_reaction)};
    }
    if constexpr (dynamics == DynamicsList::Type::firstpassage) {
      return ptof::Boundary_Cases{
          meta::Selector_t<State>{},
          boundary_conditions,
          locator,
          std::forward<VelocityField>(velocity_field),
          Boundary_Periodic{boundary_periodic, locator},
          Boundary_Reinject{std::forward<InitialCondition>(initial_condition)}};
    }
    throw std::runtime_error{std::string{"Boundary conditions for "} +
                             DynamicsList::name(dynamics) + " not supported"};
  }

  /** \brief Output generic information about object. */
  static std::ostream &info(std::ostream &output) {
    output << io::line() << "Geometry\n"
           << io::line()
           << "Spatial dimension: " + std::to_string(dim) +
                  R"(
Boundary conditions for: )" +
                  DynamicsList::name(dynamics) +
                  R"(
Boundary condition types:
  - reflecting:
    - Reflecting.
  - reacting
    - Reflection and surface reaction.
  - periodic)"
           << (PeriodicityList::name(periodicity) == "cartesian" ? R"(
    - Periodic in the primitive unit cell.)"
                                                                 : R"(
    - Periodic in the minimal unit cell.)")
           << R"(
  - absorbing
    - Absorbing
  - inlet
    - Absorbing for outward flow velocities, reflecting otherwise.)";
    if constexpr (dynamics == DynamicsList::Type::firstpassage) {
      output << R"(
  - custom
    - Reinject according to initial condition.)";
    }
    output << R"(
  - empty
    - No effect (default for unspecified patches in mesh).
)" << io::line();
    return output;
  }

  /** \brief Output information about current object. */
  std::ostream &info_runtime(std::ostream &output) const {
    output << io::line() << "Mesh\n"
           << io::line() << "Cells: " << mesh().nCells() << "\n"
           << "Faces: " << mesh().nFaces() << "\n"
           << "Edges: " << mesh().nEdges() << "\n"
           << io::line();
    return output;
  }

  /** \return Mesh. */
  auto const &mesh() const { return _mesh; }

  /** \return Mesh search tools. */
  auto const &mesh_search() const { return _mesh_search; }

private:
  Mesh _mesh;              /**< Mesh object. */
  MeshSearch _mesh_search; /**< Mesh searching tools. */

  /**
     \param directories_of OpenFOAM case directory information.
     \return Mesh.
  */
  Mesh make_mesh(DirectoriesOF const &directories_of) const {
    return Mesh{Foam::IOobject{Foam::fvMesh::defaultRegion,
                               directories_of.time.timeName(),
                               directories_of.time, Foam::IOobject::MUST_READ}};
  }

public:
  Locator locator{*this}; /**< To locate positions in mesh.*/
  const double radius;    /**< Bead radius. */
  BoundaryPeriodic
      boundary_periodic; /**< Boundary object to enforce periodicity. */

private:
  /**
     \brief Make a boundary object to enforce cartesian periodicity.
     \return Periodic boundary object.
  */
  BoundaryPeriodic make_boundary_periodic(
      meta::Selector<PeriodicityList::Type, PeriodicityList::Type::cartesian>,
      Directories const &directories) {
    return {extract_cartesian_periodic_boundaries(
        get_boundary_conditions<dynamics>(directories), mesh(), dim)};
  };

  /**
     \brief Make a boundary object to enforce periodicity based on symmetry
     planes.
     \return Periodic boundary object.
  */
  BoundaryPeriodic
  make_boundary_periodic(meta::Selector<PeriodicityList::Type,
                                        PeriodicityList::Type::symmetryplanes>,
                         Directories const &directories) {
    return {radius};
  }
};
} // namespace ptof

#endif /* PTOF_GEOMETRY_H */
