//
//  Geometry.h
//
//  Created by Tomás Aquino on 09/03/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Geometry_OF_h
#define Geometry_OF_h

#include <cstddef>
#include <string>
#include <type_traits>
#include <fvMesh.H>
#include <IOobject.H>
#include <meshSearch.H>
#include <unordered_map>
#include "ParticleTrackingOF/Boundary.h"
#include "ParticleTrackingOF/Directories.h"
#include "ParticleTrackingOF/Locator.h"
#include "ParticleTrackingOF/Store.h"
#include "Stochastic/CTRW/Boundary.h"

namespace ptof
{
  // Keep track of names of
  // types of boundary condition sets
  // for different types of simulations
  struct BoundaryConditionSet
  {
    // Implemented types
    enum class Type
    {
      transport,
      firstpassage,
      cartesian,
      symmetryplanes
    };
    
    // Type from name
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    // Name from type
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    // Check if name exists
    static bool contains(std::string const& name)
    { return name_to_type.count(name); }
    
    // Map names to types
    inline static const
    std::unordered_map<std::string, Type> name_to_type
    {
      { "transport", Type::transport },
      { "firstpassage", Type::firstpassage },
      { "cartesian", Type::cartesian },
      { "symmetryplanes", Type::symmetryplanes }
    };
    
    // Map types to names
    inline static const
    std::unordered_map<Type, std::string> type_to_name
    {
      { Type::transport, "transport" },
      { Type::firstpassage, "firstpassage" },
      { Type::cartesian, "cartesian" },
      { Type::symmetryplanes, "symmetryplanes"}
    };
  };
  
  // Geometry class
  // Holds spatial dimension info,
  // mesh, and mesh search objects,
  // and handles boundary construction
  // for different types of dynamics
  template
  <std::size_t dim_val,
  BoundaryConditionSet::Type dynamics>
  struct Geometry
  {
    static constexpr std::size_t dim{ dim_val };
    using Mesh = Foam::fvMesh;
    using MeshSearch = Foam::meshSearch;
    using Locator = Locator_Cell;
    
    // Store only absorptions for transport simulations,
    // store absorptions and reinjections for fpt simulations
    using BoundaryInfo
      = std::conditional_t<
          dynamics
            == BoundaryConditionSet::Type::transport,
          Store_Absorbed,
          Store_Absorbed_Reinjections>;
    
    // Construct with mesh from specified OpenFOAM case time
    template <typename Parameters>
    Geometry
    (DirectoriesOF const& directories,
     Parameters const& parameters = useful::Empty{})
    : mesh{
        Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories.time.timeName(),
          directories.time,
          Foam::IOobject::MUST_READ } }
    {}
    
    // Make a Boundary object for the preset type of dynamics
    // 'transport': No parameters, custom boundary does nothing
    // 'firstpassage': Parameters should hold InitialCondition object,
    // custom boundary reinjects according to initial condition
    template
    <typename TransportParameters,
    typename ReactionParameters,
    typename SolverParameters,
    typename Parameters = useful::Empty>
    auto makeBoundary
    (Directories const& directories,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     Parameters&& parameters = {}) const
    {
      if constexpr (dynamics
                    == BoundaryConditionSet::Type::transport)
        return ptof::Boundary_Cases{
          ptof::get_boundary_conditions(directories.dir_boundaryconditions
                                        + "/boundary_conditions_"
                                        + BoundaryConditionSet::name(dynamics)
                                        + ".dat"),
          mesh_search,
          BoundaryInfo{} };
      if constexpr (dynamics
                    == BoundaryConditionSet::Type::firstpassage)
        return ptof::Boundary_Cases{
          ptof::get_boundary_conditions(directories.dir_boundaryconditions
                                        + "/boundary_conditions_"
                                        + BoundaryConditionSet::name(dynamics)
                                        + ".dat"),
          mesh_search,
          BoundaryInfo{},
          Boundary_Reinject{ parameters }
        };
      throw std::runtime_error{
        std::string{ "Boundary conditions for " }
          + BoundaryConditionSet::name(dynamics)
          + " not supported" };
    }
    
    // Output generic information about object
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
    
    // Output information about current object
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
    
    Mesh mesh;                        // Mesh for fields and boundaries
    MeshSearch mesh_search{ mesh };   // To find positions in mesh, etc
    Locator locator{ mesh_search };   // Wrapper on mesh search
  };
  
  // Geometry class for periodic geometry
  // Holds spatial dimension info,
  // mesh, and mesh search objects
  template
  <std::size_t dim_val,
  BoundaryConditionSet::Type dynamics,
  BoundaryConditionSet::Type periodicity_type>
  struct Geometry_Periodic
  {
    static constexpr std::size_t dim{ dim_val };
    using Mesh = Foam::fvMesh;
    using MeshSearch = Foam::meshSearch;
    using Locator = Locator_Cell;
    using Boundary_Periodic = std::conditional_t<
      periodicity_type
      == BoundaryConditionSet::Type::cartesian,
    boundary::Periodic_WithOutsideInfo,
    boundary::Periodic_SymmetryPlanes_WithOutsideInfo<
      geometry::SymmetryPlanes_Bcc>>;
    
    // Store only absorptions
    using BoundaryInfo = Store_Absorbed;
    
    // Construct with mesh from specified OpenFOAM case time
    template <typename Parameters>
    Geometry_Periodic
    (DirectoriesOF const& directories,
     Parameters const& params)
    : mesh{
        Foam::IOobject{
          Foam::fvMesh::defaultRegion,
          directories.time.timeName(),
          directories.time,
          Foam::IOobject::MUST_READ } }
    , boundary_periodic{
        make_boundary_periodic(useful::Selector<
                                BoundaryConditionSet::Type,
                                periodicity_type>{},
                               params) }
    {}
    
    // Make a Boundary object for the preset type of periodicity
    // Parameters holds periodic boundary
    // 'cartesian': Cartesian periodicity
    // 'symmetryplanes': Periodicity according to symmetry planes
    template
    <typename TransportParameters,
    typename ReactionParameters,
    typename SolverParameters,
    typename Parameters = useful::Empty>
    auto makeBoundary
    (Directories const& directories,
     TransportParameters const& params_transport,
     ReactionParameters const& params_reaction,
     SolverParameters const& params_solvers,
     Parameters&& parameters = {})
    {
      return ptof::Boundary_Cases{
        ptof::get_boundary_conditions(directories.dir_boundaryconditions
                                      + "/boundary_conditions_"
                                      + BoundaryConditionSet::name(dynamics)
                                      + ".dat"),
        mesh_search,
        BoundaryInfo{},
        Boundary_Periodic_OF{
          boundary_periodic,
          mesh_search } };
    }
    
    // Output generic information about object
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
        "\tabsorbing: Absorbing\n"
        "\tinfo: Information upon crossing\n"
        "\tcustom: ";
      output <<
        (BoundaryConditionSet::name(periodicity_type) == "cartesian"
         ? "Periodic in the primitive unit cell\n"
         : "Periodic in the minimal unit cell\n");
      output <<
        "\tempty: No effect (default for unspecified patches in mesh)\n"
        "--------------------------------------------------\n";
    }
    
    // Output information about current object
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
    
    Mesh mesh;                            // Mesh for fields and boundaries
    MeshSearch mesh_search{ mesh };       // To find positions in mesh, etc
    Locator locator{ mesh_search };       // Wrapper on mesh search
    Boundary_Periodic boundary_periodic;  // Boundary object to enforce periodicity
  
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
        midpoint.push_back( (boundary.second - boundary.first)/2. );
      return { params.radius, midpoint };
    }
  };
}

#endif /* Geometry_OF_h */
