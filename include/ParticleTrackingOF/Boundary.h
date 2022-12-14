//
//  Boundary.h
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
#include "general/Operations.h"
#include "ParticleTrackingOF/useful.h"

namespace ptof
{
  // Keep track of names and
  // types of boundary conditions
  struct BoundaryCondition
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
    BoundaryCondition()
    {}
    
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
      { "reflecting", Type::reflecting },
      { "absorbing", Type::absorbing },
      { "info", Type::info },
      { "custom", Type::custom },
      { "empty", Type::custom },
    };
    
    // Map types to names
    inline static const
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
      if (!implemented.contains(bc.second))
        throw std::runtime_error{
          "Boundary condition type "
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
                                 -make_point(state_old.position),
                               mesh_search),
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
          case BoundaryCondition::Type::reflecting:
          {
            store_info(state, state_old, intersection,
                       boundary_condition_types,
                       useful::Selector<
                        BoundaryCondition::Type,
                        BoundaryCondition::Type::reflecting>{});
            boundary_reflecting(state,
                                intersection.point(),
                                reflection_normal(intersection.index()));
            had_effect = 1;
            break;
          }
          case BoundaryCondition::Type::custom:
          {
            store_info(state, state_old, intersection,
                       boundary_condition_types,
                       useful::Selector<
                        BoundaryCondition::Type,
                        BoundaryCondition::Type::custom>{});
            had_effect += boundary_custom(state, intersection);
            break;
          }
          case BoundaryCondition::Type::info:
          {
            store_info(state, state_old, intersection,
                       boundary_condition_types,
                       useful::Selector<
                        BoundaryCondition::Type,
                        BoundaryCondition::Type::info>{});
            break;
          }
          case BoundaryCondition::Type::absorbing:
          {
            store_info(state, state_old, intersection,
                       boundary_condition_types,
                       useful::Selector<
                        BoundaryCondition::Type,
                        BoundaryCondition::Type::absorbing>{});
            boundary_absorbing(state, intersection.point());
            had_effect = 1;
            break;
          }
          case BoundaryCondition::Type::empty:
          {}
          default:
            throw std::runtime_error{
              "Boundary condition type "
              + type_name
              + " not supported"
            };
        }
        
        if (make_point(state.position) == intersection.point())
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
    
    const BoundaryCondition
      boundary_condition_types{}; // Boundary condition names and types
    
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
    
    // Find next intersection
    // given current intersection and end of jump
    template <typename Intersection>
    auto next_intersection
    (Intersection const& intersection,
     Foam::point const& end) const
    {
      return mesh_search.
        intersection(offset_forward_face(intersection.point(),
                                         intersection.index(),
                                         end-intersection.point(),
                                         mesh_search),
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
  
  // Periodic image of intersection point
  template
  <typename State, typename Intersection,
  typename Boundary, typename MeshSearch>
  auto periodic_intersection
  (State state_outside,
   State const& state_image,
   Intersection const& intersection,
   Boundary const& boundary,
   MeshSearch const& mesh_search)
  {
    // Compute displacement from intersection to outside position
    auto displacement = make_point(state_outside.position)-
      intersection.point();
    
    // Slight offset from intersection towards the outside point
    // The goal is to find the periodic image
    // of a point close to the intersection,
    // so we can find its periodic image within the domain
    // and then use it to find the periodic image of the intersection
    // itself. Using a periodic image like this should work also for
    // non-cartesian symmetry planes.
    auto offset_pos
      = offset_forward_face(intersection.point(),
                            intersection.index(),
                            displacement,
                            mesh_search);
    
    // Ensure offset position is closer to intersection than outside position
    // If less than halfway (arbitrary point between the two),
    // set outside position to offset position
    // otherwise, use the halfway point
    if (Foam::magSqr(offset_pos - intersection.point())
        < 0.5*Foam::magSqr(displacement))
      state_outside.set_position(offset_pos);
    else
      state_outside.set_position(intersection.point()+
                                 0.5*displacement);
    
    // Find the periodic image of the offset point
    boundary(state_outside);
    
    // Find the periodic image of the intersection by
    // looking from outside the domain towards the final periodic image
    // Start at the outside point found by taking
    return mesh_search.
      intersection(make_point(state_outside.position)
                   - 2.*(make_point(state_image.position)-
                         make_point(state_outside.position)),
                   make_point(state_outside.position));
  }
  
  // Boundary object for periodic boundaries
  // The templated boundary object should enforce the periodic boundaries
  template <typename Boundary, typename MeshSearch>
  struct Boundary_Periodic_OF
  {
    // Construct given boundary positions along each dimension
    Boundary_Periodic_OF
    (Boundary&& boundary,
     MeshSearch&& mesh_search)
    : boundary_periodic{ std::forward<Boundary>(boundary) }
    , mesh_search{ std::forward<MeshSearch>(mesh_search) }
    {}
    
    // Check if position is out of bounds
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      return boundary_periodic.outOfBounds(position);
    }
    
    // Enforce boundary condition
    // and place intersection at periodic image
    template <typename State, typename Intersection>
    bool operator() (State& state, Intersection& intersection) const
    {
      auto state_outside = state;
      boundary_periodic(state);
      intersection
        = periodic_intersection(state_outside, state,
                                intersection, boundary_periodic,
                                mesh_search);
      if (!intersection.hit())
        throw std::runtime_error{
          "Could not find image of periodic boundary intersection" };
      
      return 1;
    }
    
    // Custom name of boundary condition type
    std::string name() const
    { return "periodic"; }
    
    Boundary boundary_periodic;
    MeshSearch mesh_search;
  };
  template <typename Boundary, typename MeshSearch>
  Boundary_Periodic_OF(Boundary&&, MeshSearch&&)
  ->Boundary_Periodic_OF<Boundary, MeshSearch>;
  
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
