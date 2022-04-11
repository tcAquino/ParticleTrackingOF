//
//  Boundary_OF.h
//  ParticleTracking_OpenFOAM
//
//  Created by Tomás Aquino on 17/02/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef Boundary_OF_h
#define Boundary_OF_h

#include <algorithm>
#include <cstddef>
#include <exception>
#include <fstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>
#include <boost/algorithm/string.hpp>
#include "general/useful.h"
#include "ParticleTrackingOF/useful.h"
#include "Stochastic/CTRW/Boundary.h"

namespace ptof
{
  // Keep track of names and
  // types of boundary conditions
  struct ImplementedBoundaryConditions
  {
    // Implemented types
    enum class Type
    {
      reflecting,
      absorbing,
      info,
      custom,
      empty
    };
    
    // Construct with standard 'custom' name
    ImplementedBoundaryConditions()
    {}
    
    // Construct with given name for 'custom' type
    ImplementedBoundaryConditions(std::string const& custom_name)
    { set_custom_name(custom_name); }
    
    // Change the name of 'custom' type
    void set_custom_name(std::string const& custom_name)
    { type_to_name[Type::custom] = custom_name; }
    
    // Type from name
    auto type(std::string const& name) const
    { return name_to_type.at(name); }
    
    // Name from type
    auto name(Type type) const
    { return type_to_name.at(type); }
    
    // Check if name exists
    bool exists(std::string const& name) const
    { return name_to_type.count(name); }
    
  private:
    // Map names to types
    std::unordered_map<std::string, Type> name_to_type
    {
      { "reflecting", Type::reflecting },
      { "absorbing", Type::absorbing },
      { "info", Type::info },
      { "custom", Type::custom },
      { "empty", Type::custom },
    };
    
    // Map types to names
    std::unordered_map<Type, std::string> type_to_name
    {
      { Type::reflecting, "reflecting" },
      { Type::absorbing, "absorbing" },
      { Type::info, "info" },
      { Type::custom, "custom" },
      { Type::empty, "empty" }
    };
  };
  
  // Apply reflecting boundary condition
  // given the contact point at boundary
  // and the unit normal at the contact point
  // pointing inwards
  template <typename State>
  void boundary_reflecting
  (State& state,
   Foam::point const& contact_point,
   Foam::vector const& unit_normal)
  {
    auto leftover_jump = make_point(state.position)
      - contact_point;
    
    // When a particle starts at a reflecting wall,
    // it must not be reflected if it jumps towards the interior
    // and not through the boundary
    // This dot product can only be positive in that case
    auto leftover_dot_normal = leftover_jump&unit_normal;
    state.set_position(contact_point
                       + leftover_jump
                       + (leftover_dot_normal < 0.
                          ? - 2.*leftover_dot_normal*unit_normal
                          : Foam::zero{}));
  }

  // Apply absorbing boundary condition
  template <typename State>
  void boundary_absorbing
  (State& state, Foam::point const& contact_point)
  {
    state.set_position(contact_point);
  }

  // Get the boundary conditions associated
  // with a set of patches from file
  // Each file line should hold a patch name
  // followed by a bc type
  auto get_boundary_conditions
  (std::string const& filename,
   std::size_t header_lines = 0,
   std::string const& delims = "\t,| ")
  {
    std::unordered_map<std::string, std::string>
      boundary_conditions;
    
    auto file = useful::open_read(filename);
    std::string line;
    for (std::size_t ll = 0; ll < header_lines; ++ll)
      getline(file, line);
    
    while (getline(file, line))
    {
      std::vector<std::string> split_line;
      boost::trim_if(line, boost::is_any_of(delims+"\r"));
      boost::algorithm::split(split_line, line,
                              boost::is_any_of(delims),
                              boost::token_compress_on);
      boundary_conditions[split_line[0]] = split_line[1];
    }
    file.close();
    
    return boundary_conditions;
  }
  
  // Verify that patch names exist in mesh
  // and associated boundary condition types are implemented
  // BCs should be a container of pairs of patch names and
  // boundary condition types
  template <typename BCs, typename Mesh, typename Implemented>
  void verify_boundary_conditions
  (BCs const& boundary_conditions, Mesh const& mesh,
   Implemented const& implemented)
  {
    for (auto const& bc : boundary_conditions)
      mesh.boundaryMesh().findPatchID(bc.first, false);
    
    for (auto const& bc : boundary_conditions)
      if (!implemented.exists(bc.second))
        throw std::runtime_error{
          std::string("Boundary condition type ")
          + bc.second
          + " not supported"
        };
  }
  
  // Add default type 'empty'
  // to boundary conditions for unspecified patches in patch_names
  // BCs should be an associative container of pairs of patch names and
  // boundary condition types
  template <typename BCs, typename Container>
  void add_unspecified_patches
  (BCs& boundary_conditions,
   Container const& patch_names)
  {
    for (auto const& name : patch_names)
      if (!boundary_conditions.count(name))
        boundary_conditions[name] = "empty";
  }
  
  // Add default type 'empty'
  // to boundary conditions for unspecified patches in mesh
  // BCs should be a container of pairs of patch names and
  // boundary condition types
  template <typename BCs, typename Mesh>
  void add_unspecified_patches_in_mesh
  (BCs& boundary_conditions, Mesh const& mesh)
  {
    add_unspecified_patches(boundary_conditions,
                            mesh.boundaryMesh().names());
  }
  
  // Boundary object with no effect
  struct Boundary_DoNothing
  {
    // Check if position is out of bounds (it never is)
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return 0; }
    
    // Apply boundary condition (do nothing, return 0 for no effect)
    template <typename State, typename Intersection>
    bool operator() (State& state, Intersection const&) const
    { return 0; }
    
    // Custom name of boundary condition type
    std::string name() const
    { return "empty"; }
  };
  
  // Boundary object to handle implemented boundary types
  template
  <typename MeshSearch, typename Store_Info, typename Boundary_Custom>
  class Boundary_Cases
  {
  public:
    // Container type to hold patch names and associated bc type names
    using BCs = std::unordered_map<std::string, std::string>;
    
    // Construct given
    // boundary condition patches and associated types,
    // Mesh search object,
    // Object to store info when applying different boundary condition types
    // Custom boundary object for boundary condition types
    Boundary_Cases
    (BCs boundary_conditions,
     MeshSearch&& mesh_search,
     Store_Info&& store_info,
     Boundary_Custom&& boundary_custom = Boundary_DoNothing{})
    : boundary_conditions{ boundary_conditions }
    , mesh_search{ std::forward<MeshSearch>(mesh_search) }
    , store_info{ std::forward<Store_Info>(store_info) }
    , boundary_custom{ std::forward<Boundary_Custom>(boundary_custom) }
    {
      add_unspecified_patches(this->boundary_conditions,
                              patch_names());
      verify_boundary_conditions(this->boundary_conditions,
                                 mesh_search.mesh(),
                                 boundary_condition_types);
    }
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      Foam::label cell = mesh_search.findCell(make_point(position));
      return cell < 0;
    }

    // Enforce boundary condition by choosing appropriate type
    // Return 1 if some boundary had an effect, 0 otherwise
    template <typename State>
    bool operator()
    (State& state, State const& state_old = {}) const
    {
      // Find first intersection with boundary patch
      auto intersection =
        mesh_search.intersection(make_point(state_old.position),
                                 make_point(state.position));
      
      // When the start point is on a cell face,
      // sometimes the intersection with it is not found
      // Avoid breaking by checking for intersections with a small
      // backwards offset when the final state is out of bounds
      state.cell = mesh_search.findCell(make_point(state.position),
                                        state.cell);
      if (!intersection.hit() && state.cell == -1)
        intersection =
          mesh_search.intersection(
          offset_backward_cell(make_point(state_old.position),
                               state_old.cell,
                               make_point(state.position)
                                 -make_point(state_old.position)),
                               make_point(state.position));
      
      bool had_effect = 0;
      
      // While some patch is intersected
      while (intersection.hit())
      {
        // Find patch in mesh and associated boundary type name
        auto patch_id = mesh_search.mesh().boundaryMesh().
          whichPatch(intersection.index());
        auto patch = patch_name(patch_id);
        auto type_name = boundary_conditions.at(patch);
        
        // Apply the approriate boundary and store associated info
        switch (boundary_condition_types.type(type_name))
        {
          case ImplementedBoundaryConditions::Type::reflecting:
          {
            store_info(state, state_old, intersection,
                       boundary_condition_types,
                       useful::Selector<
                        ImplementedBoundaryConditions::Type,
                        ImplementedBoundaryConditions::Type::reflecting>{});
            boundary_reflecting(state,
                                intersection.rawPoint(),
                                reflection_normal(intersection.index()));
            had_effect = 1;
            break;
          }
          case ImplementedBoundaryConditions::Type::custom:
          {
            store_info(state, state_old, intersection,
                       boundary_condition_types,
                       useful::Selector<
                        ImplementedBoundaryConditions::Type,
                        ImplementedBoundaryConditions::Type::custom>{});
            had_effect += boundary_custom(state, intersection);
            break;
          }
          case ImplementedBoundaryConditions::Type::info:
          {
            store_info(state, state_old, intersection,
                       boundary_condition_types,
                       useful::Selector<
                        ImplementedBoundaryConditions::Type,
                        ImplementedBoundaryConditions::Type::info>{});
            break;
          }
          case ImplementedBoundaryConditions::Type::absorbing:
          {
            store_info(state, state_old, intersection,
                       boundary_condition_types,
                       useful::Selector<
                        ImplementedBoundaryConditions::Type,
                        ImplementedBoundaryConditions::Type::absorbing>{});
            boundary_absorbing(state, intersection.rawPoint());
            had_effect = 1;
            break;
          }
          case ImplementedBoundaryConditions::Type::empty:
          {}
          default:
            throw std::runtime_error{
              std::string("Boundary condition type ")
              + type_name
              + " not supported"
            };
        }
        
        if (make_point(state.position) == intersection.rawPoint())
          break;
        
        // Find the next intersection with a boundary patch
        intersection = next_intersection(intersection,
                                         make_point(state.position));
      }
      
      return had_effect;
    }
    
    // Output information about current object
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Boundary conditions\n"
        "--------------------------------------------------\n"
        "Patch / Condition\n"
        "\n";
      for (auto const& bc : boundary_conditions)
      {
        output << bc.first << "\t";
        if (bc.second == "custom")
          output << boundary_custom.name();
        else
          output << bc.second;
        output << "\n";
      }
      output <<
        "--------------------------------------------------\n";
    }
    
  private:
    BCs boundary_conditions;          // Patch names and associated bc type names
    MeshSearch mesh_search;           // Mesh search object to find points, intersections, etc
    Store_Info store_info;            // Object to handle boundary info storing in states
    Boundary_Custom boundary_custom;  // Boundary object to handle 'custom' bc type
    
    const ImplementedBoundaryConditions
      boundary_condition_types{
          boundary_custom.name() }; // Implemnted boundary condition names and types
    
    // Names of patches in mesh
    auto patch_names() const
    { return mesh_search.mesh().boundaryMesh().names(); }
    
    // Name of given patch index in mesh
    auto patch_name(std::size_t patch) const
    { return patch_names()[patch]; }
    
    // Wrapper function to get normal to face for reflection
    auto reflection_normal(Foam::label face) const
    {
      return unit_normal_inward(face, mesh_search.mesh());
    }
    
    // Small offset forward given current face and direction
    Foam::vector offset_forward_face
    (Foam::point const& begin,
     Foam::label face,
     Foam::vector const& direction) const
    {
      Foam::label owner_cell = mesh_search.mesh().faceOwner()[face];
      Foam::point const& center = mesh_search.mesh().
        cellCentres()[owner_cell];
      Foam::scalar typDim = Foam::mag(center - begin);
      
      return begin + mesh_search.tol_*typDim*
        direction/Foam::mag(direction);
    }
    
    // Small offset backward
    // given current cell and direction
    Foam::vector offset_backward_cell
    (Foam::point const& begin,
     Foam::label cell,
     Foam::vector const& direction) const
    {
      Foam::point const& center = mesh_search.mesh().cellCentres()[cell];
      Foam::scalar typDim = Foam::mag(center - begin);
      
      return begin - mesh_search.tol_*typDim*
        direction/Foam::mag(direction);
    }
    
    // Find next intersection
    // given current intersection and end of jump
    template <typename Intersection>
    auto next_intersection
    (Intersection const& intersection,
     Foam::point const& end) const
    {
      return mesh_search.
        intersection(offset_forward_face(intersection.rawPoint(),
                                         intersection.index(),
                                         end-intersection.rawPoint()),
                     end);
    }
  };
  template
  <typename MeshSearch,
  typename Store_Info,
  typename Boundary_Custom>
  Boundary_Cases
  (typename Boundary_Cases<
    MeshSearch, Store_Info, Boundary_Custom>::BCs,
   MeshSearch&&,
   Store_Info&&,
   Boundary_Custom&&)->
  Boundary_Cases<MeshSearch, Store_Info, Boundary_Custom>;
  template <typename MeshSearch, typename Store_Info>
  Boundary_Cases
  (typename Boundary_Cases<
    MeshSearch, Store_Info, Boundary_DoNothing>::BCs,
   MeshSearch&&,
   Store_Info&&)->
  Boundary_Cases<MeshSearch, Store_Info, Boundary_DoNothing>;
  
  // Boundary object for
  // Periodic boundaries along each dimension
  // State must define: position (vector type)
  struct Boundary_Periodic
  {
    // Construct given boundary positions along each dimension
    Boundary_Periodic
    (std::vector<std::pair<double, double>> boundaries)
    : boundary{ boundaries }
    {}
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      return boundary.outOfBounds(position);
    }
    
    // Enforce boundary condition
    template <typename State, typename Intersection>
    bool operator() (State& state, Intersection const&) const
    {
      return boundary(state);
    }
    
    // Custom name of boundary condition type
    std::string name() const
    { return "periodic"; }
    
    boundary::Periodic boundary;
  };
  
  // Boundary object for
  // Periodic boundaries along each symmetry plane
  // State must define: position (vector type)
  template <typename SymmetryPlanes>
  struct Boundary_Periodic_SymmetryPlanes
  {
    // Construct given symmetry plane object,
    // overall scale factor, and coordinate origin
    Boundary_Periodic_SymmetryPlanes
    (SymmetryPlanes symmetry_planes, double scale = 1.,
     std::vector<double> origin = {})
    : boundary{ symmetry_planes, scale, origin }
    {}
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      return boundary.outOfBounds(position);
    }
    
    // Enforce boundary condition
    template <typename State, typename Intersection>
    bool operator() (State& state, Intersection const&) const
    {
      return boundary(state);
    }
    
    // Custom name of boundary condition type
    std::string name() const
    { return "periodic"; }
    
    boundary::Periodic_SymmetryPlanes<
      SymmetryPlanes> boundary;
  };
  
  // Boundary object for
  // Periodic boundaries along each dimension
  // with information about where position would be outside domain
  // State must define: position (vector type)
  //                    periodicity (std::vector<int>) counting how many domains
  //                                                   have been traveled along each dimension
  struct Boundary_Periodic_WithOutsideInfo
  {
    // Construct given boundary positions along each dimension
    Boundary_Periodic_WithOutsideInfo
    (std::vector<std::pair<double, double>> boundaries)
    : boundary{ boundaries }
    {}
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      return boundary.outOfBounds(position);
    }
    
    // Enforce boundary condition
    template <typename State, typename Intersection>
    bool operator() (State& state, Intersection const&) const
    {
      return boundary(state);
    }
    
    // Custom name of boundary condition type
    std::string name() const
    { return "periodic"; }
    
    boundary::Periodic_WithOutsideInfo boundary;
  };
  
  template <typename SymmetryPlanes>
  struct Boundary_Periodic_SymmetryPlanes_WithOutsideInfo
  {
    // Construct given symmetry plane object,
    // overall scale factor, and coordinate origin
    Boundary_Periodic_SymmetryPlanes_WithOutsideInfo
    (SymmetryPlanes symmetry_planes, double scale = 1.,
     std::vector<double> origin = {})
    : boundary{ symmetry_planes, scale, origin }
    {}
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position)
    {
      return boundary.outOfBounds(position);
    }
    
    // Enforce boundary condition
    template <typename State, typename Intersection>
    bool operator() (State& state, Intersection const&) const
    {
      return boundary(state);
    }
    
    // Custom name of boundary condition type
    std::string name() const
    { return "periodic"; }
    
    boundary::Periodic_SymmetryPlanes_WithOutsideInfo<
      SymmetryPlanes> boundary;
  };
  
  // Boundary object that reinjects particles according to
  // prescribed initial condition
  template <typename InitialCondition>
  struct Boundary_Reinject
  {
    // Construct given initial condition object
    Boundary_Reinject
    (InitialCondition&& initial_condition)
    : initial_condition{
        std::forward<InitialCondition>(initial_condition) }
    {}
    
    // Enforce boundary condition
    template <typename State, typename Intersection>
    bool operator() (State& state, Intersection const&) const
    {
      state.set_position(initial_condition.make_position());
      return 1;
    }
    
    // Custom name of boundary condition type
    std::string name() const
    { return "reinject"; }
    
    InitialCondition initial_condition;
  };
  template <typename InitialCondition>
  Boundary_Reinject
  (InitialCondition&&)->
  Boundary_Reinject<InitialCondition>;
}


#endif /* Boundary_OF_h */
