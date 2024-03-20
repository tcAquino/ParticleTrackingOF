/**
 \file PTOF/Geometry.h
 \author Tomás Aquino
 \date 09/03/2022
*/
#ifndef PTOF_GEOMETRY_H
#define PTOF_GEOMETRY_H

#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#include <fvMesh.H>
#include <IOobject.H>
#include <meshSearch.H>
#include "General/Constants.h"
#include "Geometry/Boundary.h"
#include "PTOF/Boundary.h"
#include "PTOF/Directories.h"
#include "PTOF/Locator.h"
#include "PTOF/Store.h"

namespace ptof
{
  /** \struct BoundaryConditionSet PTOF/Geometry.h "PTOF/Geometry.h"
     \brief Keep track of names of types of boundary condition sets  for different types of simulations. */
  struct BoundaryConditionSet
  {
    /** \enum Type
     *  \brief Implemented types. */  
    enum class Type
    {
      transport,       /**< No parameters, custom boundary does nothing.   */
      firstpassage,    /**< Parameters should hold InitialCondition object.*/
      cartesian,       /**< Cartesian periodicity.                         */
      symmetryplanes   /**< Periodicity according to symmetry planes.      */
    };
    
    /** \brief Get type from name. */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** \brief Get name from type. */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** \brief Check if name exists. */
    static bool contains(std::string const& name)
    { return name_to_type.count(name); }
    
    /** \brief Map names to types. */
    inline static const
    std::unordered_map<std::string, Type> name_to_type
    {
      { "transport", Type::transport },
      { "firstpassage", Type::firstpassage },
      { "cartesian", Type::cartesian },
      { "symmetryplanes", Type::symmetryplanes }
    };
    
    /** \brief Map types to names. */
    inline static const
    std::unordered_map<Type, std::string> type_to_name
    {
      { Type::transport, "transport" },
      { Type::firstpassage, "firstpassage" },
      { Type::cartesian, "cartesian" },
      { Type::symmetryplanes, "symmetryplanes" }
    };
  };
  
  /**
   \brief Get boundary conditions for dynamics type.
   \param directories Current case directory information.
   */
  template <BoundaryConditionSet::Type dynamics>
  auto get_boundary_conditions(Directories const& directories)
  {
    if constexpr (dynamics == BoundaryConditionSet::Type::transport)
      return
        ptof::get_boundary_conditions(directories.dir_boundaryconditions
                                      + "/boundary_conditions_"
                                      + BoundaryConditionSet::name(dynamics)
                                      + ".dat");
    if constexpr (dynamics == BoundaryConditionSet::Type::firstpassage)
      return
        ptof::get_boundary_conditions(directories.dir_boundaryconditions
                                      + "/boundary_conditions_"
                                      + BoundaryConditionSet::name(dynamics)
                                      + ".dat");
  }
  
  /** \class Geometry PTOF/Geometry.h "PTOF/Geometry.h"
   * \brief Basic geometry class.
   \note Does not handle periodic boundaries. */
  template
  <std::size_t dim_val,
  BoundaryConditionSet::Type dynamics = BoundaryConditionSet::Type::transport,
  typename SearchOption = SearchOptions::SecondNeighborPrecheck>
  struct Geometry
  {
    static constexpr std::size_t dim{ dim_val };    /**< Spatial dimension. */
    using Mesh = Foam::fvMesh;                      /**< Mesh. */
    using MeshSearch = Foam::meshSearch;            /**< Mesh searching tools. */
    using Locator = Locator_Cell<SearchOption>;     /**< Locate positions in mesh.*/
    using BoundaryInfo
      = std::conditional_t<
          dynamics
            == BoundaryConditionSet::Type::transport,
          Store_Absorbed,
          Store_Absorbed_Reinjections>;             /**< What to store upon hitting a boundary. */
    
    /** Constructor.
    \brief Use mesh from specified OpenFOAM case time.
    \param directories_of OpenFOAM case directory information.
    \param directories Current case directory information.
    */
    Geometry
    (DirectoriesOF const& directories_of,
     Directories const& directories)
    : _mesh{
        Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ } }
    {}
    
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
    \param output Output stream
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
    
    /** \return Mesh. */
    auto const& mesh() const
    { return _mesh; }
    
  private:
    Mesh _mesh;                          /**< Mesh for fields and boundaries. */
    MeshSearch _mesh_search{ _mesh };    /**< Mesh searching tools.*/
    
  public:
    Locator locator{ _mesh_search };    /**< To locate positions in mesh.*/
  };
  
  /** \class Geometry_Periodic_Cartesian PTOF/Geometry.h "PTOF/Geometry.h"
   *  \brief Geometry class for geometry to handle periodic BCs along any number of Cartesian dimensions.*/
  template
  <std::size_t dim_val,
  BoundaryConditionSet::Type dynamics,
  typename SearchOption = SearchOptions::SecondNeighborPrecheck>
  struct Geometry_Periodic_Cartesian
  {
    static constexpr std::size_t dim{ dim_val };    /**< Spatial dimension. */
    using Mesh = Foam::fvMesh;                      /**< Mesh. */
    using MeshSearch = Foam::meshSearch;            /**< Mesh searching tools. */
    using Locator = Locator_Cell<SearchOption>;     /**< Locate positions in mesh.*/
    using Boundary_Periodic =
      geometry::Boundary_Periodic_WithOutsideInfo;  /**< Periodic boundary. */
    using BoundaryInfo = Store_Absorbed;            /**< What to store upon hitting a boundary. */
    
    /** Constructor.
     \brief Use mesh from specified OpenFOAM case time.
     \param directories_of OpenFOAM case directory information.
     \param directories Current case directory information.
    */
    Geometry_Periodic_Cartesian
    (DirectoriesOF const& directories_of,
     Directories const& directories)
    : _mesh{
        Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ } }
    , boundary_periodic{
        extract_cartesian_periodic_boundaries(get_boundary_conditions<dynamics>(directories),
                                              mesh(),
                                              dim) }
    {}
    
    /**
     \brief Make a Boundary object for Cartesian periodicity.
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
    
    /** \return Mesh. */
    auto const& mesh() const
    { return _mesh; }
    
  private:
    Mesh _mesh;                            /**< Mesh for fields and boundaries.        */
    MeshSearch _mesh_search{ _mesh };       /**< Mesh searching tools. */
    
  public:
    Locator locator{ _mesh_search };       /**< To locate positions in mesh.*/
    Boundary_Periodic boundary_periodic;  /**< Boundary object to enforce periodicity.*/
    
  private:
    /**
    \brief Make a boundary object to enforce periodicity based on cartesian planes.
    \param directories Current case directory information.
    \return Periodic boundary object.
    */
  };
  
  /** \class Geometry_Bcc PTOF/Geometry.h "PTOF/Geometry.h"
   \brief Geometry class for periodic representation of a body centered cubic beadpack.
   \details Holds spatial dimension info, mesh, mesh search objects, and periodic boundary to enforce periodicity. */
  template
  <std::size_t dim_val,
  BoundaryConditionSet::Type dynamics,
  BoundaryConditionSet::Type periodicity_type,
  typename SearchOption = SearchOptions::SecondNeighborPrecheck>
  struct Geometry_Bcc
  {
    static constexpr std::size_t dim{ dim_val };    /**< Spatial dimension. */
    using Mesh = Foam::fvMesh;                      /**< Mesh. */
    using MeshSearch = Foam::meshSearch;            /**< Mesh searching tools. */
    using Locator = Locator_Cell<SearchOption>;     /**< Locate positions in mesh.*/
    using Boundary_Periodic = std::conditional_t<
      periodicity_type
      == BoundaryConditionSet::Type::cartesian,
    geometry::Boundary_Periodic_WithOutsideInfo,
    geometry::Boundary_Periodic_SymmetryPlanes_WithOutsideInfo<
      geometry::SymmetryPlanes_Bcc>>;              /**< Periodic boundary type. */
    using BoundaryInfo = Store_Absorbed;           /**< What to store upon hitting a boundary. */
    
    /** Constructor.
    \brief Use mesh from specified OpenFOAM case time.
    \param directories_of OpenFOAM case directory information.
    \param directories Current case directory information.
    */
    Geometry_Bcc
    (DirectoriesOF const& directories_of,
     Directories const& directories)
    : _mesh{
        Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories_of.time.timeName(),
          directories_of.time,
          Foam::IOobject::MUST_READ } }
    , radius{
      std::is_same_v<periodicity_type, BoundaryConditionSet::Type::cartesian>
      ? std::pow(sum(mesh().V())/(1.-std::sqrt(3.)*constants::pi/8.), 1./3.)*std::sqrt(3.)/4.
      : 1. }
    , boundary_periodic{
        make_boundary_periodic(useful::Selector<
                                BoundaryConditionSet::Type,
                                periodicity_type>{},
                               directories) }
    {}
    
    /**
     \brief Make a Boundary object for the preset type of dynamics and periodicity
     \details
     - 'cartesian': Cartesian periodicity
     - 'symmetryplanes': Periodicity according to symmetry planes
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
    
    /** \return Mesh. */
    auto const& mesh() const
    { return _mesh; }
    
  private:
    Mesh _mesh;                             /**< Mesh for fields and boundaries. */
    MeshSearch _mesh_search{ _mesh };       /**< Mesh searching tools. */
    
  public:
    const double radius;                    /**< Bead radius. */
    Locator locator{ _mesh_search };        /**< To locate positions in mesh.*/
    Boundary_Periodic boundary_periodic;    /**< Boundary object to enforce periodicity. */
  
  private:
    /**
    \brief Make a boundary object to enforce cartesian periodicity.
    \return Periodic boundary object.
    */
    Boundary_Periodic make_boundary_periodic
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
    Boundary_Periodic make_boundary_periodic
    (useful::Selector<BoundaryConditionSet::Type,
     BoundaryConditionSet::Type::symmetryplanes>,
     Directories const& directories)
    { return { radius }; }
  };
}

#endif /* PTOF_GEOMETRY_H */
