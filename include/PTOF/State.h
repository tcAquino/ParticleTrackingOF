/**
 \file PTOF/State.h
 \author Tomás Aquino
 \date 21/02/2022
*/

#ifndef PTOF_STATE_H
#define PTOF_STATE_H

#include <cstddef>
#include <exception>
#include <fieldTypes.H>
#include <type_traits>
#include <vector>
#include "General/Useful.h"
#include "PTOF/Useful.h"

namespace ptof
{
  /** \class State1D PTOF/State.h "PTOF/State.h"
   \brief State for 1D positions.
   \details Defines:
   - \c position
   - \c cell
   - \c time
   - \c mass
   - \c tag
   - \c info */
  template
  <typename Info_t,
  typename Time_t = useful::Empty,
  typename Mass_t = useful::Empty,
  typename Tag_t = useful::Empty>
  struct State1D
  {
    using Position = Foam::scalar;
    using Index = Foam::label;
    using Time = Time_t;
    using Mass = Mass_t;
    using Tag = Tag_t;
    using Info = Info_t;
    
    Position position{ 0. };
    Index cell{ -1 };
    Time time{};
    Mass mass{};
    Tag tag{};
    Info info{};
    
    /** \brief Set position from 3D position. */
    void set_position(Foam::point const& point)
    { position = point[0]; }
    
    /** \brief Make position from 3D position. */
    static Position make_position(Foam::point const& point)
    { return point[0]; }
  };

  /** \class State2D PTOF/State.h "PTOF/State.h"
   \brief State for 2D positions.
   \details Defines:
   - \c position
   - \c cell
   - \c time
   - \c mass
   - \c tag
   - \c info */
  template
  <typename Info_t,
  typename Time_t = useful::Empty,
  typename Mass_t = useful::Empty,
  typename Tag_t = useful::Empty>
  struct State2D
  {
    using Position = Foam::Vector2D<Foam::scalar>;
    using Index = Foam::label;
    using Time = Time_t;
    using Mass = Mass_t;
    using Tag = Tag_t;
    using Info = Info_t;
    
    Position position{ 0., 0. };
    Index cell{ -1 };
    Time time{};
    Mass mass{};
    Tag tag{};
    Info info{};
    
    /** \brief Set position from 3D position. */
    void set_position(Foam::point const& point)
    {
      position[0] = point[0];
      position[1] = point[1];
    }
    
    /** \return Position from 3D position. */
    static Position make_position(Foam::point const& point)
    { return Position{ point[0], point[1] }; }
  };

  /** \class State3D PTOF/State.h "PTOF/State.h"
   \brief State for 3D positions.
   \details Defines:
   - \c position
   - \c cell
   - \c time
   - \c mass
   - \c tag
   - \c info */
  template
  <typename Info_t,
  typename Time_t = useful::Empty,
  typename Mass_t = useful::Empty,
  typename Tag_t = useful::Empty>
  struct State3D
  {
    // Typedefs for state quantities
    using Position = Foam::point;
    using Index = Foam::label;
    using Time = Time_t;
    using Mass = Mass_t;
    using Tag = Tag_t;
    using Info = Info_t;
    
    Position position{ 0., 0., 0. };
    Index cell{ -1 };
    Time time{};
    Mass mass{};
    Tag tag{};
    Info info{};
    
    /** \brief Set position from 3D position. */
    void set_position(Foam::point const& point)
    { position = point; }
    
    /** \return Position from 3D position. */
    static Position make_position(Foam::point const& point)
    { return point; }
  };

  /** \return State for 1D positions.
   \details Locate mesh cell using \c Locator object. */
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
    if (outside(cell_id))
      throw std::runtime_error{
        "Particle initialized outside mesh" };
    return { position, cell_id, time, mass, tag };
  }

  /** \return State for 2D positions.
    \details Locate initial mesh cell using \c Locator object. */
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
    if (outside(cell_id))
      throw std::runtime_error{
        "Particle initialized outside mesh" };
    return { position, cell_id, time, mass, tag };
  }

  /** \return State for 3D positions.
   \details Locate initial mesh cell using \c Locator object. */
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
    if (outside(cell_id))
      throw std::runtime_error{
        "Particle initialized outside mesh" };
    return { position, cell_id, time, mass, tag };
  }
  
  /** \return State for 1D positions.
  \details Initial cell index to be provided. */
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

  /** \return State for 2D positions.
  \details Initial cell index to be provided. */
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

  /** \return State for 3D positions.
  \details Initial cell index to be provided. */
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
  
  /** \brief State type templated on spatial dimension. */
  template
  <std::size_t dim, typename Info,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  using State = typename std::conditional_t<
    dim == 3,
    State3D<Info, Time, Mass, Tag>,
    std::conditional_t<
      dim == 2,
      State2D<Info, Time, Mass, Tag>,
      State1D<Info, Time, Mass, Tag>>>;
  
  /** \class State1D_Periodic PTOF/State.h "PTOF/State.h"
  \brief State for 1D positions with periodicity info.
  \details State periodicity holds signed number of cells away from unit cell along each basis direction.
  \details Defines:
  - \c position
  - \c periodicity
  - \c cell
  - \c time
  - \c mass
  - \c tag
  - \c info */
  template
  <typename Info_t,
  typename Time_t = useful::Empty,
  typename Mass_t = useful::Empty,
  typename Tag_t = useful::Empty>
  struct State1D_Periodic
  {
    using Position = Foam::scalar;
    using Periodicity = std::vector<int>;
    using Index = Foam::label;
    using Info = Info_t;
    using Time = Time_t;
    using Mass = Mass_t;
    using Tag = Tag_t;
    
    Position position{ 0. };
    Periodicity periodicity{ 0 };
    Index cell{ -1 };
    Time time{};
    Mass mass{};
    Tag tag{};
    
    /** \brief Set position from 3D position. */
    void set_position(Foam::point const& point)
    { position = point[0]; }
    
    /** \return Position from 3D position. */
    static Position make_position(Foam::point const& point)
    { return point[0]; }
    
    Info info{};
  };

  /** \class State2D_Periodic PTOF/State.h "PTOF/State.h"
  \brief State for 2D positions with periodicity info.
  \details State periodicity holds signed number of cells away from unit cell along each basis direction.
  \details Defines:
  - \c position
  - \c periodicity
  - \c cell
  - \c time
  - \c mass
  - \c tag
  - \c info */
  template
  <typename Info_t,
  typename Time_t = useful::Empty,
  typename Mass_t = useful::Empty,
  typename Tag_t = useful::Empty>
  struct State2D_Periodic
  {
    // Typedefs for state quantities
    using Position = Foam::Vector2D<Foam::scalar>;
    using Periodicity = std::vector<int>;
    using Index = Foam::label;
    using Time = Time_t;
    using Mass = Mass_t;
    using Tag = Tag_t;
    using Info = Info_t;
    
    Position position{ 0., 0. };
    Periodicity periodicity{ 0, 0 };
    Index cell{ -1 };
    Time time{};
    Mass mass{};
    Tag tag{};
    Info info{};
    
    /** \brief Set position from 3D position. */
    void set_position(Foam::point const& point)
    {
      position[0] = point[0];
      position[1] = point[1];
    }
    
    /** \return Position from 3D position. */
    static Position make_position(Foam::point const& point)
    { return Position{ point[0], point[1] }; }
  };

  /** \class State3D_Periodic PTOF/State.h "PTOF/State.h"
  \brief State for 3D positions with periodicity info
  \details State periodicity holds signed number of cells away from unit cell along each basis direction.
  \details Defines:
  - \c position
  - \c periodicity
  - \c cell
  - \c time
  - \c mass
  - \c tag
  - \c info */
  template
  <typename Info_t,
  typename Time_t = useful::Empty,
  typename Mass_t = useful::Empty,
  typename Tag_t = useful::Empty>
  struct State3D_Periodic
  {
    // Typedefs for state quantities
    using Position = Foam::point;
    using Periodicity = std::vector<int>;
    using Index = Foam::label;
    using Info = Info_t;
    using Time = Time_t;
    using Mass = Mass_t;
    using Tag = Tag_t;
    
    Position position{ 0., 0., 0. };
    Periodicity periodicity{ 0, 0, 0 };
    Index cell{ -1 };
    Time time{};
    Mass mass{};
    Tag tag{};
    
    /** Set position from 3D position. */
    void set_position(Foam::point const& point)
    { position = point; }
    
    /** Make position from 3D position. */
    static Position make_position(Foam::point const& point)
    { return point; }
    
    Info info{};
  };

  /** \return State for 1D positions.
  \details Locate mesh cell using \c Locator object. */
  template
  <typename Info,
  typename Locator,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State1D_Periodic<Info, Time, Mass, Tag> make_state
  (Foam::scalar position,
   std::vector<int> periodicity,
   Info info,
   Locator const& locator,
   Time time = {},
   Mass mass = {},
   Tag tag = {})
  {
    Foam::label cell_id = locator(position);
    if (outside(cell_id))
      throw std::runtime_error{
        "Particle initialized outside mesh" };
    return { position, periodicity, cell_id, time, mass, tag };
  }

  /** \return State for 2D positions.
  \details Locate mesh cell using \c Locator object. */
  template
  <typename Info,
  typename Locator,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State2D_Periodic<Info, Time, Mass, Tag> make_state
  (Foam::Vector2D<Foam::scalar> const& position,
   std::vector<int> periodicity,
   Info info,
   Locator const& locator,
   Time time = {},
   Mass mass = {},
   Tag tag ={})
  {
    Foam::label cell_id = locator(position);
    if (outside(cell_id))
      throw std::runtime_error{
        "Particle initialized outside mesh" };
    return { position, periodicity, cell_id, time, mass, tag };
  }

  /** \return State for 3D positions.
  \details Locate mesh cell using \c Locator object. */
  template
  <typename Info,
  typename Locator,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State3D_Periodic<Info, Time, Mass, Tag> make_state
  (Foam::vector const& position,
   std::vector<int> periodicity,
   Info info,
   Locator const& locator,
   Time time = {},
   Mass mass = {},
   Tag tag = {})
  {
    Foam::label cell_id = locator(position);
    if (outside(cell_id))
      throw std::runtime_error{
        "Particle initialized outside mesh" };
    return { position, periodicity, cell_id, time, mass, tag };
  }
  
  /** Make state for 1D positions with periodicity info
   * with given cell index. */
  template
  <typename Info,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State1D_Periodic<Info, Time, Mass, Tag> make_state
  (Foam::scalar position,
   std::vector<int> periodicity,
   Info info,
   Foam::label cell_id,
   Time time = {},
   Mass mass = {},
   Tag tag = {})
  {
    return { position, periodicity, cell_id, time, mass, tag };
  }

  /** \return State for 2D positions.
  \details Initial cell index to be provided. */
  template
  <typename Info,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State2D_Periodic<Info, Time, Mass, Tag> make_state
  (Foam::Vector2D<Foam::scalar> const& position,
   std::vector<int> periodicity,
   Foam::label cell_id,
   Time time = {},
   Mass mass = {},
   Tag tag ={})
  {
    return { position, periodicity, cell_id, time, mass, tag };
  }

  /** \return State for 3D positions.
  \details Initial cell index to be provided. */
  template
  <typename Info,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  State3D_Periodic<Info, Time, Mass, Tag> make_state
  (Foam::vector const& position,
   std::vector<int> periodicity,
   Foam::label cell_id,
   Time time = {},
   Mass mass = {},
   Tag tag = {})
  {
    return { position, periodicity, cell_id, time, mass, tag };
  }
  
  /** \brief State type templated on spatial dimension. */
  template
  <std::size_t dim, typename Info,
  typename Time = useful::Empty,
  typename Mass = useful::Empty,
  typename Tag = useful::Empty>
  using State_Periodic = typename std::conditional_t<
    dim == 3,
    State3D_Periodic<Info, Time, Mass, Tag>,
    std::conditional_t<
      dim == 2,
      State2D_Periodic<Info, Time, Mass, Tag>,
      State1D_Periodic<Info, Time, Mass, Tag>>>;
}

#endif /* PTOF_STATE_H */
