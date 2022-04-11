//
//  State_OF.h
//  ParticleTracking_OpenFOAM
//
//  Created by Tomás Aquino on 21/02/2022.
//  Copyright © 2022 Tomás Aquino. All rights reserved.
//

#ifndef State_OF_h
#define State_OF_h

#include <cstddef>
#include <fieldTypes.H>
#include "ParticleTrackingOF/Locator.h"
#include "general/useful.h"

namespace ptof
{
  // State for 1D positions
  template
  <typename Info_t,
  typename Time_t = useful::Empty,
  typename Mass_t = useful::Empty,
  typename Tag_t = useful::Empty>
  struct State1D
  {
    // Typedefs for state quantities
    using Position = Foam::scalar;
    using Index = Foam::label;
    using Info = Info_t;
    using Time = Time_t;
    using Mass = Mass_t;
    using Tag = Tag_t;
    
    Position position{ 0. };
    Index cell{ -1 };
    Time time{};
    Mass mass{};
    Tag tag{};
    
    // Set position from 3D position
    void set_position(Foam::point const& point)
    { position = point[0]; }
    
    // Make position from 3D position
    static auto make_position(Foam::point const& point)
    { return point[0]; }
    
    Info info{};
  };

  // State for 2D positions
  template
  <typename Info_t,
  typename Time_t = useful::Empty,
  typename Mass_t = useful::Empty,
  typename Tag_t = useful::Empty>
  struct State2D
  {
    // Typedefs for state quantities
    using Position = Foam::Vector2D<Foam::scalar>;
    using Index = Foam::label;
    using Info = Info_t;
    using Time = Time_t;
    using Mass = Mass_t;
    using Tag = Tag_t;
    
    Position position{ 0., 0. };
    Index cell{ -1 };
    Time time{};
    Mass mass{};
    Tag tag{};
    
    // Set position from 3D position
    void set_position(Foam::point const& point)
    {
      position[0] = point[0];
      position[1] = point[1];
    }
    
    // Make position from 3D position
    static auto make_position(Foam::point const& point)
    { return Position{ point[0], point[1] }; }
    
    Info info{};
  };

  // State for 3D positions
  template
  <typename Info_t,
  typename Time_t = useful::Empty,
  typename Mass_t = useful::Empty,
  typename Tag_t = useful::Empty>
  struct State3D
  {
    // Typedefs for state quantities
    using Position = Foam::vector;
    using Index = Foam::label;
    using Info = Info_t;
    using Time = Time_t;
    using Mass = Mass_t;
    using Tag = Tag_t;
    
    Position position{ 0., 0., 0. };
    Index cell{ -1 };
    Time time{};
    Mass mass{};
    Tag tag{};
    
    // Set position from 3D position
    void set_position(Foam::point const& point)
    { position = point; }
    
    // Make position from 3D position
    static auto make_position(Foam::point const& point)
    { return point; }
    
    Info info{};
  };

  // Make state for 1D positions
  // Locate mesh cell using Locator object
  template
  <typename Info,
  typename Locator,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State1D<Info, Time, Mass, Tag> make_state
  (Foam::scalar position,
   Info info,
   Locator const& locator,
   Time time = {},
   Mass mass = {},
   Tag tag = {})
  {
    Foam::label cell_id = locator(position);
    if (cell_id == -1)
      throw std::runtime_error{
        "Particle initialized outside mesh" };
    return { position, cell_id, time, mass, tag };
  }

  // Make state for 2D positions
  // Locate mesh cell using Locator object
  template
  <typename Info,
  typename Locator,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State2D<Info, Time, Mass, Tag> make_state
  (Foam::Vector2D<Foam::scalar> const& position,
   Info info,
   Locator const& locator,
   Time time = {},
   Mass mass = {},
   Tag tag ={})
  {
    Foam::label cell_id = locator(position);
    if (cell_id == -1)
      throw std::runtime_error{
        "Particle initialized outside mesh" };
    return { position, cell_id, time, mass, tag };
  }

  // Make state for 3D positions
  // Locate mesh cell using Locator object
  template
  <typename Info,
  typename Locator,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State3D<Info, Time, Mass, Tag> make_state
  (Foam::vector const& position,
   Info info,
   Locator const& locator,
   Time time = {},
   Mass mass = {},
   Tag tag = {})
  {
    Foam::label cell_id = locator(position);
    if (cell_id == -1)
      throw std::runtime_error{
        "Particle initialized outside mesh" };
    return { position, cell_id, time, mass, tag };
  }
  
  // Make state for 1D positions with given cell index
  template
  <typename Info,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State1D<Info, Time, Mass, Tag> make_state
  (Foam::scalar position,
   Info info,
   Foam::label cell_id,
   Time time = {},
   Mass mass = {},
   Tag tag = {})
  {
    return { position, cell_id, time, mass, tag };
  }

  // Make state for 2D positions with given cell index
  template
  <typename Info,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State2D<Info, Time, Mass, Tag> make_state
  (Foam::Vector2D<Foam::scalar> const& position,
   Foam::label cell_id,
   Time time = {},
   Mass mass = {},
   Tag tag ={})
  {
    return { position, cell_id, time, mass, tag };
  }

  // Make state for 3D positions with given cell index
  template
  <typename Info,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State3D<Info, Time, Mass, Tag> make_state
  (Foam::vector const& position,
   Foam::label cell_id,
   Time time = {},
   Mass mass = {},
   Tag tag = {})
  {
    return { position, cell_id, time, mass, tag };
  }
  
  // Typedef state type templated on spatial dimension
  template
  <std::size_t dim, typename Info,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  using StateDim = typename std::conditional_t<
    dim == 3,
    State3D<Info, Time, Mass, Tag>,
    std::conditional_t<
      dim == 2,
      State2D<Info, Time, Mass, Tag>,
      State1D<Info, Time, Mass, Tag>>>;
  
  // Object to make particles given position
  template <typename Particle, typename Locator>
  struct ParticleMaker
  {
    // Typedefs for particle state quantities
    using State = typename Particle::State;
    using Info = typename State::Info;
    using Time = typename State::Time;
    using Mass = typename State::Mass;
    using Tag = typename State::Tag;
    
    // Construct given Locator to find mesh cell,
    // time, mass, and tag
    ParticleMaker
    (useful::Selector_t<Particle>,
     Locator&& locator,
     Time time = {},
     Mass mass = {},
     Tag tag = {})
    : locator{ std::forward<Locator>(locator) }
    , time{ time }
    , mass{ mass }
    , tag{ tag }
    {}
    
    // Make particle given position
    template <typename Position>
    Particle operator()
    (Position const& position)
    {
      return { make_state(State::make_position(position),
                          Info{}, locator, time, mass, tag) };
    }
    
    Locator const& locator; // Locator to locate particle states in mesh
    Time time;              // Particle state time
    Mass mass;              // Particle state mass
    Tag tag;                // Particle state tag
  };
  template <typename Particle, typename Locator>
  ParticleMaker
  (useful::Selector_t<Particle>,
   Locator&& locator,
   typename Particle::State::Time,
   typename Particle::State::Mass,
   typename Particle::State::Tag) ->
  ParticleMaker<Particle, Locator>;
  template <typename Particle, typename Locator>
  ParticleMaker
  (useful::Selector_t<Particle>,
   Locator&& locator,
   typename Particle::State::Time,
   typename Particle::State::Mass) ->
  ParticleMaker<Particle, Locator>;
  template <typename Particle, typename Locator>
  ParticleMaker
  (useful::Selector_t<Particle>,
   Locator&& locator,
   typename Particle::State::Time) ->
  ParticleMaker<Particle, Locator>;
  template <typename Particle, typename Locator>
  ParticleMaker
  (useful::Selector_t<Particle>,
   Locator&& locator) ->
  ParticleMaker<Particle, Locator>;
}


#endif /* State_OF_h */
