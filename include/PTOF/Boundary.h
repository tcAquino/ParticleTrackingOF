/**
 \file PTOF/Boundary.h
 \author Tomás Aquino
 \date 17/02/2022
*/

#ifndef PTOF_BOUNDARY_H
#define PTOF_BOUNDARY_H

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
#include "General/Useful.h"
#include "General/Operations.h"
#include "PTOF/Reaction.h"
#include "PTOF/Useful.h"

namespace ptof
{
   /** \struct BoundaryCondition PTOF/Boundary.h "PTOF/Boundary.h"
    \brief Names and types of boundary conditions. */
  struct BoundaryCondition
  {
    /** \enum Type
     *  \brief Implemented types. */
    enum class Type
    {
      reflecting,           /**< Reflecting                                          */
      reacting_reflecting,  /**< Reacting and reflecting                     */
      periodic,             /**< Periodic                                            */
      absorbing,            /**< Absorbing                                           */
      info,                 /**< Information upon hitting                      */
      custom,               /**< Enforced by custom Boundary object                  */
      empty                 /**< No effect (default for unspecified patches in mesh) */
    };
    
    /** Constructor.
     \brief Default constructor. */
    BoundaryCondition()
    {}
    
    /** \brief Type from name.
     \param name Boundary condition name.
     \return Boundary condition type.
     */
    static auto type(std::string const& name)
    { return name_to_type.at(name); }
    
    /** \brief Name from type.
     \param type Boundary condition type.
     \return Boundary condition name.
    */
    static auto name(Type type)
    { return type_to_name.at(type); }
    
    /** \brief Check if name exists.
     \param name condition name.
     \return \c true if name exists, \c false otherwise.
     */
    static bool contains(std::string const& name)
    { return name_to_type.count(name); }
    
    /** \brief Map names to types. */
    inline static const
    std::unordered_map<std::string, Type> name_to_type
    {
      { "reflecting", Type::reflecting },
      { "reacting_reflecting", Type::reacting_reflecting },
      { "periodic", Type::periodic },
      { "absorbing", Type::absorbing },
      { "info", Type::info },
      { "custom", Type::custom },
      { "empty", Type::empty },
    };
    
    /** \brief Map types to names. */
    inline static const
    std::unordered_map<Type, std::string> type_to_name
    {
      { Type::reflecting, "reflecting" },
      { Type::reacting_reflecting, "reacting_reflecting" },
      { Type::periodic, "periodic" },
      { Type::absorbing, "absorbing" },
      { Type::info, "info" },
      { Type::custom, "custom" },
      { Type::empty, "empty" }
    };
  };
  
  /**
   \brief Apply reflecting boundary condition.
   \param state Particle state to update.
   \param contact_point Contact point at boundary.
   \param unit_normal Unit normal at the contact point, pointing inwards. */
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

  /**
  \brief Apply absorbing boundary condition.
  \param state Particle state to update.
  \param contact_point Contact point at boundary. */
  template <typename State>
  void boundary_absorbing
  (State& state, Foam::point const& contact_point)
  {
    state.set_position(contact_point);
  }

  /**
   \brief Get the boundary conditions associated with a set of patches from file.
   \param filename Name of file holding a patch name followed by a boundary condition type in each line
   \param header_lines Number of header lines to ignore in file
   \param delims Possible delimiters between patch names and boundary condition types
  */
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
  
  /** \brief Throw if any patch name does not exist in mesh or associated boundary condition type is not implemented.
   \param boundary_conditions Container of pairs of patch names and boundary condition types.
   \param mesh OpenFOAM mesh.
   \param implemented BoundaryCondition object holding information about implemented boundary conditions */
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
  
  /**
   \brief Add default type 'empty' to boundary conditions for unspecified patches
   \param boundary_conditions Associative container of pairs of patch names and boundary condition types, to add to.
   \param patch_names Container with all patch names
  */
  template <typename BCs, typename Container>
  void add_unspecified_patches
  (BCs& boundary_conditions,
   Container const& patch_names)
  {
    for (auto const& name : patch_names)
      if (!boundary_conditions.count(name))
        boundary_conditions[name] = "empty";
  }
  
  /**
   \brief Add default type 'empty' to boundary conditions for unspecified patches in mesh
   \param boundary_conditions Associative container of pairs of patch names and boundary condition types, to add to.
   \param mesh OpenFOAM mesh
  */
  template <typename BCs, typename Mesh>
  void add_unspecified_patches_in_mesh
  (BCs& boundary_conditions, Mesh const& mesh)
  {
    add_unspecified_patches(boundary_conditions,
                            mesh.boundaryMesh().names());
  }
  
  /** \struct Boundary_DoNothing PTOF/Boundary.h "PTOF/Boundary.h"
   * \brief Boundary object with no effect. */
  struct Boundary_DoNothing
  {
    /**
     \brief Check if position is out of bounds.
     \param position Position to check.
     \return \c false (always in bounds).
     */
    template <typename Position>
    bool outOfBounds(Position const& position) const
    { return false; }
    
    /**
     \brief Apply boundary condition (do nothing).
     \param state Current particle state.
     \param state_old Previous particle state.
     \param dummy Dummy argument (unused)
     \return \c false (no effect).
    */
    template <typename State, typename Dummy>
    bool operator()
    (State const& state, State const& state_old,
     Dummy const& dummy) const
    { return false; }
    
    /**
     \brief Apply boundary condition (do nothing).
     \param state Current particle state.
     \param dummy Dummy argument (unused)
     \return \c false (no effect).
    */
    template <typename State, typename Dummy>
    bool operator()
    (State const& state, Dummy const& dummy) const
    { return false; }
    
    /**
     \brief Apply boundary condition (do nothing).
     \param state Current particle state.
     \return \c false (no effect).
    */
    template <typename State>
    bool operator()
    (State const& state) const
    { return false; }
    
    /** \brief Name of boundary condition type. */
    std::string name() const
    { return "empty"; }
  };
  
  /** \class Boundary_Cases PTOF/Boundary.h "PTOF/Boundary.h"
   *  \brief Boundary object to handle implemented boundary types */
  template
  <typename Locator, typename Store_Info,
  typename Boundary_Periodic, typename Boundary_Custom,
  typename SurfaceReaction>
  class Boundary_Cases
  {
  public:
    /** \brief Container type to hold patch names and associated bc type names. */
    using BCs = std::unordered_map<std::string, std::string>;
    
    /** Constructor
     \param boundary_conditions Patch names and associated  boundary condition types.
     \param locator Object to locate positions in mesh.
     \param store_info Object to handle storing of information upon boundary hitting.
     \param boundary_periodic Periodic boundary condition enforcer.
     \param boundary_custom Custom boundary condition enforcer.
     \param surface_reaction Surface Reaction.
    */
    Boundary_Cases
    (BCs boundary_conditions,
     Locator&& locator,
     Store_Info&& store_info,
     Boundary_Periodic&& boundary_periodic = Boundary_DoNothing{},
     Boundary_Custom&& boundary_custom = Boundary_DoNothing{},
     SurfaceReaction&& surface_reaction = SurfaceReaction_DoNothing{})
    : _boundary_conditions{ boundary_conditions }
    , _locator{ std::forward<Locator>(locator) }
    , _store_info{ std::forward<Store_Info>(store_info) }
    , _boundary_periodic{ std::forward<Boundary_Periodic>(boundary_periodic) }
    , _boundary_custom{ std::forward<Boundary_Custom>(boundary_custom) }
    , surface_reaction{ std::forward<SurfaceReaction>(surface_reaction) }
    {
      add_unspecified_patches(_boundary_conditions,
                              patch_names());
      verify_boundary_conditions(_boundary_conditions,
                                 _locator.mesh_search().mesh(),
                                 _boundary_condition_types);
    }
    
    /**
     \brief Check if position is out of bounds.
     \param position Position to check
     \return \c true if out of bounds, \c false otherwise.
    */
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      Foam::label cell = _locator.mesh_search().findCell(make_point(position));
      return outside(cell);
    }

    /**
     \brief Enforce boundary conditions if necessary by choosing appropriate types.
     \param state Current particle state to apply BCs if needed (possibly out of bounds).
     \param state_old Previous particle state (should be in bounds).
     \return \c true if some boundary had an effect, \c false otherwise.
    */
    template <typename State>
    bool operator()
    (State& state, State const& state_old = {})
    {
      /** Find first intersection with boundary patch. */
      auto intersection =
        _locator.mesh_search().intersection(make_point(state_old.position),
                                            make_point(state.position));
      
      // When the start point is on a cell face,
      // sometimes the intersection with it is not found
      // Avoid breaking by checking for intersections with a small
      // backwards offset when the final state is out of bounds
      if (!intersection.hit())
      {
        state.cell = _locator(state);
        if (outside(state.cell))
        {
          intersection =
            _locator.mesh_search().intersection(offset_backward_cell(make_point(state_old.position),
                                                                     state_old.cell,
                                                                     make_point(state.position)
                                                                     - make_point(state_old.position),
                                                                     _locator),
                                                make_point(state.position));
        }
      }
      else
      {
        // Ignore new intersection if offset went beyond final point
        if (((intersection.point() - make_point(state_old.position))
                & (make_point(state.position) - make_point(state_old.position))) < 0.)
          intersection.setMiss();
      }
      
      bool had_effect = 0;
      
      // While some patch is intersected
      while (intersection.hit())
      {
        /** Find patch in mesh and associated boundary type name. */
        auto patch_id = _locator.mesh_search().mesh().boundaryMesh().
          whichPatch(intersection.index());
        auto patch = patch_name(patch_id);
        auto type_name = _boundary_conditions.at(patch);
        
        // Apply the approriate boundary and store associated info
        switch (_boundary_condition_types.type(type_name))
        {
          case BoundaryCondition::Type::reflecting:
          {
           _store_info(state, state_old, intersection,
                       _boundary_condition_types,
                       useful::Selector<
                         BoundaryCondition::Type,
                         BoundaryCondition::Type::reflecting>{});
            boundary_reflecting(state,
                                intersection.point(),
                                reflection_normal(intersection.index()));
            had_effect = 1;
            break;
          }
          case BoundaryCondition::Type::reacting_reflecting:
          {
            _store_info(state, state_old, intersection,
                        _boundary_condition_types,
                        useful::Selector<
                          BoundaryCondition::Type,
                          BoundaryCondition::Type::reflecting>{});
            surface_reaction(state, state_old, intersection.index());
            boundary_reflecting(state,
                                intersection.point(),
                                reflection_normal(intersection.index()));
            had_effect = 1;
            break;
          }
          case BoundaryCondition::Type::periodic:
          {
            _store_info(state, state_old, intersection,
                        _boundary_condition_types,
                        useful::Selector<
                          BoundaryCondition::Type,
                          BoundaryCondition::Type::periodic>{});
            had_effect += _boundary_periodic(state, intersection);
            break;
          }
          case BoundaryCondition::Type::absorbing:
          {
            _store_info(state, state_old, intersection,
                        _boundary_condition_types,
                        useful::Selector<
                          BoundaryCondition::Type,
                          BoundaryCondition::Type::absorbing>{});
            boundary_absorbing(state, intersection.point());
            had_effect = 1;
            break;
          }
          case BoundaryCondition::Type::custom:
          {
            _store_info(state, state_old, intersection,
                        _boundary_condition_types,
                        useful::Selector<
                          BoundaryCondition::Type,
                          BoundaryCondition::Type::custom>{});
            had_effect += _boundary_custom(state, state_old, intersection);
            break;
          }
          case BoundaryCondition::Type::info:
          {
            _store_info(state, state_old, intersection,
                        _boundary_condition_types,
                        useful::Selector<
                          BoundaryCondition::Type,
                          BoundaryCondition::Type::info>{});
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
        
        // Stop if right at the last intersection point
        // This avoids numerical issues for some edge cases
        if (make_point(state.position) == intersection.point())
          break;
        
        // Find the next intersection with a boundary patch
        auto old_intersection_point = intersection.point();
        intersection = next_intersection(intersection,
                                         make_point(state.position),
                                         state.cell);
        // Ignore new intersection if offset went beyond final point
        if (intersection.hit()
            && ((intersection.point() - old_intersection_point)
                & (make_point(state.position) - old_intersection_point)) < 0.)
          intersection.setMiss();
      }
      
      return had_effect;
    }
    
    /** \brief Output information about current object. */
    template <typename OStream>
    void info_runtime(OStream& output) const
    {
      output <<
        "--------------------------------------------------\n"
        "Boundary conditions\n"
        "--------------------------------------------------\n"
        "Patch / Condition\n"
        "\n";
      for (auto const& bc : _boundary_conditions)
      {
        output << bc.first << "\t";
        if (bc.second == "custom")
          output << _boundary_custom.name();
        else
          output << bc.second;
        output << "\n";
      }
      output <<
        "--------------------------------------------------\n";
    }
    
  private:
    BCs _boundary_conditions;              /**< Patch names and associated bc type names.         */
    Locator _locator;                      /**< Object to locate positions in mesh. */
    Store_Info _store_info;                /**< Object to handle boundary info storing in states. */
    Boundary_Periodic _boundary_periodic;  /**< Boundary object to handle 'periodic' bc type.     */
    Boundary_Custom _boundary_custom;      /**< Boundary object to handle 'custom' bc type.       */
    
  public:
    SurfaceReaction surface_reaction;              /**< Surface reaction for reactive boundaries*/
    
  private:
    const BoundaryCondition
      _boundary_condition_types{};         /**< Boundary condition names and types. */
    
    /** \return Names of patches in mesh. */
    auto patch_names() const
    { return _locator.mesh_search().mesh().boundaryMesh().names(); }
    
    /**
     \param patch Index of patch in mesh.
     \return Name of patch. */
    auto patch_name(std::size_t patch) const
    { return patch_names()[patch]; }
    
    /**
     \param face Face index
     \return Inward normal to face (used for reflection). */
    auto reflection_normal(Foam::label face) const
    {
      return unit_normal_inward(face, _locator.mesh_search().mesh());
    }
    
    /** \brief Find next intersection with a boundary.
     \details Find the next intersection, if any, along straight line from previous intersection to endpoint.
     \param intersection previous_intersection.
     \param end Endpoint.
     \param cell Cell index for previous intersection.
     \return Next intersection. */
    template <typename Intersection>
    auto next_intersection
    (Intersection const& intersection,
     Foam::point const& end,
     Foam::label cell) const
    {
      return _locator.mesh_search().
        intersection(offset_forward_face_keep_inside(intersection.point(),
                                                     intersection.index(),
                                                     end-intersection.point(),
                                                     _locator),
                     end);
    }
  };
  template
  <typename Locator,
  typename Store_Info,
  typename Boundary_Periodic,
  typename Boundary_Custom,
  typename Surface_Reaction>
  Boundary_Cases
  (typename Boundary_Cases<Locator, Store_Info,
   Boundary_Periodic, Boundary_Custom, Surface_Reaction>::BCs,
   Locator&&,
   Store_Info&&,
   Boundary_Periodic&&,
   Boundary_Custom&&,
   Surface_Reaction&&)->
  Boundary_Cases<Locator,
    Store_Info, Boundary_Periodic, Boundary_Custom, Surface_Reaction>;
  template
  <typename Locator,
  typename Store_Info,
  typename Boundary_Periodic,
  typename Boundary_Custom>
  Boundary_Cases
  (typename Boundary_Cases<Locator, Store_Info,
   Boundary_Periodic, Boundary_Custom, SurfaceReaction_DoNothing>::BCs,
   Locator&&,
   Store_Info&&,
   Boundary_Periodic&&,
   Boundary_Custom&&)->
  Boundary_Cases<Locator, Store_Info,
    Boundary_Periodic, Boundary_Custom, SurfaceReaction_DoNothing>;
  template
  <typename Locator,
  typename Store_Info,
  typename Boundary_Periodic>
  Boundary_Cases
  (typename Boundary_Cases<Locator, Store_Info,
   Boundary_Periodic, Boundary_DoNothing, SurfaceReaction_DoNothing>::BCs,
   Locator&&,
   Store_Info&&,
   Boundary_Periodic&&)->
  Boundary_Cases<Locator, Store_Info,
    Boundary_Periodic, Boundary_DoNothing, SurfaceReaction_DoNothing>;
  template <typename Locator, typename Store_Info>
  Boundary_Cases
  (typename Boundary_Cases<Locator, Store_Info,
   Boundary_DoNothing, Boundary_DoNothing, SurfaceReaction_DoNothing>::BCs,
   Locator&&,
   Store_Info&&)->
  Boundary_Cases<Locator,
    Store_Info, Boundary_DoNothing, Boundary_DoNothing, SurfaceReaction_DoNothing>;
  
  /** \brief Find periodic image of intersection point.
   \param state_outside Out-of-bounds state.
   \param state_image Periodic image of state, in bounds.
   \param intersection Boundary intersection info.
   \param boundary Periodic boundary enforcer.
   \param locator Object to locate positions in mesh.
   \return Periodic image of intersection. */
  template
  <typename State, typename Intersection,
  typename Boundary, typename Locator>
  auto periodic_intersection
  (State state_outside,
   State const& state_image,
   Intersection const& intersection,
   Boundary const& boundary,
   Locator const& locator)
  {
    // Compute displacement from intersection to outside position.
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
                            locator);
      
    // Ensure offset position is closer to
    // intersection than outside position
    // If less than halfway
    // (arbitrary point between the two; 0.25 of mag squared),
    // set outside position to offset position
    // otherwise, use the halfway point
    if (Foam::magSqr(offset_pos - intersection.point())
        < 0.25*Foam::magSqr(displacement))
      state_outside.set_position(offset_pos);
    else
      state_outside.set_position(intersection.point()+
                                 0.5*displacement);
    
    // Find the periodic image of the offset point
    boundary(state_outside);
    
    // Find the periodic image of the intersection by
    // looking from outside the domain towards the final periodic image
    // Start at the outside point found by taking
    return locator.mesh_search().
      intersection(make_point(state_outside.position)
                   - 2.*(make_point(state_image.position)-
                         make_point(state_outside.position)),
                   make_point(state_outside.position));
  }
  
  /** \class Boundary_Periodic_OF PTOF/Boundary.h "PTOF/Boundary.h"
   *  \brief Interface for periodic boundary enforcer. */
  template <typename Boundary, typename Locator>
  struct Boundary_Periodic_OF
  {
    /** Constructor
     \param boundary Periodic boundary enforcer.
     \param locator Object to locate positions in mesh.
    */
    Boundary_Periodic_OF
    (Boundary&& boundary,
     Locator&& locator)
    : boundary_periodic{ std::forward<Boundary>(boundary) }
    , locator{ std::forward<Locator>(locator) }
    {}
    
    /**
     \brief Check if position is out of bounds.
     \param position Position to check
     \return \c true if out of bounds, \c false otherwise.
    */
    template <typename Position>
    bool outOfBounds(Position const& position) const
    {
      return boundary_periodic.outOfBounds(position);
    }
    
    /** \brief Enforce boundary condition  and place intersection at its periodic image.
     \param state Particle state to apply BC to.
     \param intersection Boundary intersection info, to be replaced by periodic image.
     \return \c true
     */
    template <typename State, typename Intersection>
    bool operator() (State& state, Intersection& intersection) const
    {
      auto state_outside = state;
      boundary_periodic(state);
      intersection
        = periodic_intersection(state_outside, state,
                                intersection, boundary_periodic,
                                locator);
      if (!intersection.hit())
        throw std::runtime_error{
          "Could not find image of periodic boundary intersection" };
      
      return 1;
    }
    
    Boundary boundary_periodic; /**< Periodic boundary type. */
    Locator locator;            /**< Object to locate positions in mesh. */
  };
  template <typename Boundary, typename Locator>
  Boundary_Periodic_OF(Boundary&&, Locator&&)
  ->Boundary_Periodic_OF<Boundary, Locator>;
  
  /** \class Boundary_Reinject PTOF/Boundary.h "PTOF/Boundary.h"
   \brief Boundary object that reinjects particles according to prescribed initial condition. */
  template <typename InitialCondition>
  struct Boundary_Reinject
  {
    /** Constructor.
     \param initial_condition Initial condition object. */
    Boundary_Reinject
    (InitialCondition&& initial_condition)
    : initial_condition{
        std::forward<InitialCondition>(initial_condition) }
    {}
    
    /** \brief Enforce boundary condition
    \param state Particle state to apply BC to.
    \param state_old Old particle state (unused).
    \param intersection Boundary intersection info (unused).
    \return \c true
    */
    template <typename State, typename Intersection>
    bool operator()
    (State& state, State const& state_old, Intersection const& intersection) const
    {
      state.set_position(initial_condition.make_position());
      return 1;
    }
    
    /** \brief Name of boundary condition type. */
    std::string name() const
    { return "reinject"; }
    
    InitialCondition initial_condition; /**< Initial condition according to which to reinject. */
  };
  template <typename InitialCondition>
  Boundary_Reinject
  (InitialCondition&&)->
  Boundary_Reinject<InitialCondition>;
}


#endif /* PTOF_BOUNDARY_H */
