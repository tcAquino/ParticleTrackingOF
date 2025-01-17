/**
   \file PTOF/State.h
   \author Tomás Aquino
   \date 21/02/2022
   \brief Particle states.
*/

#ifndef PTOF_STATE_H
#define PTOF_STATE_H

#include "General/Meta.h"
#include "PTOF/Useful.h"
#include <Vector2D.H>
#include <array>
#include <cstddef>
#include <fieldTypes.H>
#include <point.H>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace ptof {
/**
   \class State1D PTOF/State.h "PTOF/State.h"
   \brief State for 1D positions.
   \details Defines:
   - \c position
   - \c cell
   - \c time
   - \c mass
   - \c tag
   - \c info
*/
template <typename Info_t = meta::Empty, typename Time_t = meta::Empty,
          typename Mass_t = meta::Empty, typename Tag_t = meta::Empty>
struct State1D {
  using Position = Foam::scalar;
  using Cell = Foam::label;
  using Info = Info_t;
  using Time = Time_t;
  using Mass = Mass_t;
  using Tag = Tag_t;

  Position position{0.};
  Cell cell{-1};
  Info info{};
  Time time{};
  Mass mass{};
  Tag tag{};

  /** \brief Set position from 3D point. */
  void set_position(Foam::point const &point) { position = point[0]; }

  /** \brief Set position. */
  void set_position(Position const &position) { this->position = position; }

  /** \brief Set position. */
  template <typename Position> void set_position(Position const &position) {
    make_position(position);
  }

  /** \brief Set position and cell. */
  template <typename Position>
  void set_position(Position const &position, Cell cell) {
    set_position(position);
    this->cell = cell;
  }

  /** \brief Make position from 3D point. */
  static Position make_position(Foam::point const &point) { return point[0]; }

  /** \brief Make position. */
  template <typename PositionType>
  static Position make_position(PositionType const &position) {
    if (position.size() > 0)
      return position[0];
    else
      return 0.;
  }
};

/**
   \class State2D PTOF/State.h "PTOF/State.h"
   \brief State for 2D positions.
   \details Defines:
   - \c position
   - \c cell
   - \c time
   - \c mass
   - \c tag
   - \c info
*/
template <typename Info_t = meta::Empty, typename Time_t = meta::Empty,
          typename Mass_t = meta::Empty, typename Tag_t = meta::Empty>
struct State2D {
  using Position = Foam::Vector2D<Foam::scalar>;
  using Cell = Foam::label;
  using Info = Info_t;
  using Time = Time_t;
  using Mass = Mass_t;
  using Tag = Tag_t;

  Position position{0., 0.};
  Cell cell{-1};
  Info info{};
  Time time{};
  Mass mass{};
  Tag tag{};

  /** \brief Set position from 3D point. */
  void set_position(Foam::point const &point) {
    position[0] = point[0];
    position[1] = point[1];
  }

  /** \brief Set position. */
  void set_position(Position const &position) { this->position = position; }

  /** \brief Set position. */
  template <typename Position> void set_position(Position const &position) {
    make_position(position);
  }

  /** \brief Set position and cell. */
  template <typename Position>
  void set_position(Position const &position, Cell cell) {
    set_position(position);
    this->cell = cell;
  }
  /** \return Position from 3D point. */
  static Position make_position(Foam::point const &point) {
    return Position{point[0], point[1]};
  }

  /** \brief Make position. */
  template <typename PositionType>
  static Position make_position(PositionType const &position) {
    if (position.size() > 1)
      return {position[0], position[1]};
    else if (position.size() > 0)
      return {position[0], 0.};
    else
      return {0., 0.};
  }
};

/**
   \class State3D PTOF/State.h "PTOF/State.h"
   \brief State for 3D positions.
   \details Defines:
   - \c position
   - \c cell
   - \c time
   - \c mass
   - \c tag
   - \c info
*/
template <typename Info_t = meta::Empty, typename Time_t = meta::Empty,
          typename Mass_t = meta::Empty, typename Tag_t = meta::Empty>
struct State3D {
  using Position = Foam::point;
  using Cell = Foam::label;
  using Info = Info_t;
  using Time = Time_t;
  using Mass = Mass_t;
  using Tag = Tag_t;

  Position position{0., 0., 0.};
  Cell cell{-1};
  Info info{};
  Time time{};
  Mass mass{};
  Tag tag{};

  /** \brief Set position from 3D point. */
  void set_position(Foam::point const &point) { position = point; }

  /** \brief Set position. */
  template <typename Position> void set_position(Position const &position) {
    make_position(position);
  }

  /** \brief Set position and cell. */
  template <typename Position>
  void set_position(Position const &position, Cell cell) {
    set_position(position);
    this->cell = cell;
  }
  /** \return Position from 3D point. */
  static Position make_position(Foam::point const &point) { return point; }

  /** \brief Make position. */
  template <typename PositionType>
  static Position make_position(PositionType const &position) {
    if (position.size() > 2)
      return {position[0], position[1], position[2]};
    else if (position.size() > 1)
      return {position[0], 0., 0.};
    else
      return {0., 0., 0.};
  }
};

/**
   \return State for 1D positions.
   \details Locate mesh cell using \c Locator object.
*/
template <typename Info, typename Locator, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State1D<Info, Time, Mass, Tag>
make_state(Foam::scalar position, Info const &info, Locator const &locator,
           Foam::label hint = -1, Time time = {}, Mass const &mass = {},
           Tag tag = {}) {
  Foam::label cell_id = locator(position, hint);
  if (outside(cell_id))
    throw std::runtime_error{"Particle initialized outside mesh"};
  return {position, cell_id, info, time, mass, tag};
}

/**
   \return State for 2D positions.
   \details Locate initial mesh cell using \c Locator object.
*/
template <typename Info, typename Locator, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State2D<Info, Time, Mass, Tag>
make_state(Foam::Vector2D<Foam::scalar> const &position, Info const &info,
           Locator const &locator, Foam::label hint = -1, Time time = {},
           Mass mass = {}, Tag tag = {}) {
  Foam::label cell_id = locator(position, hint);
  if (outside(cell_id))
    throw std::runtime_error{"Particle initialized outside mesh"};
  return {position, cell_id, info, time, mass, tag};
}

/**
   \return State for 3D positions.
   \details Locate initial mesh cell using \c Locator object.
*/
template <typename Info, typename Locator, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State3D<Info, Time, Mass, Tag>
make_state(Foam::vector const &position, Info const &info,
           Locator const &locator, Foam::label hint = -1, Time time = {},
           Mass const &mass = {}, Tag tag = {}) {
  Foam::label cell_id = locator(position, hint);
  if (outside(cell_id))
    throw std::runtime_error{"Particle initialized outside mesh"};
  return {position, cell_id, info, time, mass, tag};
}

/**
   \return State for 1D positions.
   \details Initial cell index to be provided.
*/
template <typename Info, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State1D<Info, Time, Mass, Tag> make_state(Foam::scalar position, Info info,
                                          Foam::label cell_id, Time time = {},
                                          Mass mass = {}, Tag tag = {}) {
  return {position, cell_id, info, time, mass, tag};
}

/**
   \return State for 2D positions.
   \details Initial cell index to be provided.
*/
template <typename Info, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State2D<Info, Time, Mass, Tag>
make_state(Foam::Vector2D<Foam::scalar> const &position, Info const &info,
           Foam::label cell_id, Time time = {}, Mass const &mass = {},
           Tag tag = {}) {
  return {position, cell_id, info, time, mass, tag};
}

/**
   \return State for 3D positions.
   \details Initial cell index to be provided.
*/
template <typename Info, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State3D<Info, Time, Mass, Tag>
make_state(Foam::vector const &position, Info const &info, Foam::label cell_id,
           Time time = {}, Mass const &mass = {}, Tag tag = {}) {
  return {position, cell_id, info, time, mass, tag};
}

/** \brief Generic state type templated on spatial dimension. */
template <std::size_t dim, typename Info, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
using State_Generic = typename std::conditional_t<
    dim == 3, State3D<Info, Time, Mass, Tag>,
    std::conditional_t<dim == 2, State2D<Info, Time, Mass, Tag>,
                       State1D<Info, Time, Mass, Tag>>>;

/**
   \class State1D_Periodic PTOF/State.h "PTOF/State.h"
   \brief State for 1D positions with periodicity info.
   \details State periodicity holds signed number of cells away from unit cell
   along each basis direction. Defines:
   - \c position
   - \c periodicity
   - \c cell
   - \c time
   - \c mass
   - \c tag
   - \c info
*/
template <typename Info_t = meta::Empty, typename Time_t = meta::Empty,
          typename Mass_t = meta::Empty, typename Tag_t = meta::Empty>
struct State1D_Periodic {
  using Position = Foam::scalar;
  using Periodicity = int;
  using Cell = Foam::label;
  using Info = Info_t;
  using Time = Time_t;
  using Mass = Mass_t;
  using Tag = Tag_t;

  Position position{0.};
  Periodicity periodicity{};
  Cell cell{-1};
  Info info{};
  Time time{};
  Mass mass{};
  Tag tag{};

  /** \brief Set position from 3D point. */
  void set_position(Foam::point const &point) { position = point[0]; }

  /** \brief Set position. */
  void set_position(Position const &position) { this->position = position; }

  /** \brief Set position. */
  template <typename Position> void set_position(Position const &position) {
    make_position(position);
  }

  /** \brief Set position and cell. */
  template <typename Position>
  void set_position(Position const &position, Cell cell) {
    set_position(position);
    this->cell = cell;
  }

  /** \brief Make position from 3D point. */
  static Position make_position(Foam::point const &point) { return point[0]; }

  /** \brief Make position. */
  template <typename PositionType>
  static Position make_position(PositionType const &position) {
    if (position.size() > 0)
      return position[0];
    else
      return 0.;
  }
};

/**
   \class State2D_Periodic PTOF/State.h "PTOF/State.h"
   \brief State for 2D positions with periodicity info.
   \details State periodicity holds signed number of cells away from unit cell
   along each basis direction.
   \details Defines:
   - \c position
   - \c periodicity
   - \c cell
   - \c time
   - \c mass
   - \c tag
   - \c info
*/
template <typename Info_t = meta::Empty, typename Time_t = meta::Empty,
          typename Mass_t = meta::Empty, typename Tag_t = meta::Empty>
struct State2D_Periodic {
  using Position = Foam::Vector2D<Foam::scalar>;
  using Periodicity = std::array<int, 2>;
  using Cell = Foam::label;
  using Info = Info_t;
  using Time = Time_t;
  using Mass = Mass_t;
  using Tag = Tag_t;

  Position position{0., 0.};
  Periodicity periodicity{};
  Cell cell{-1};
  Info info{};
  Time time{};
  Mass mass{};
  Tag tag{};

  /** \brief Set position from 3D point. */
  void set_position(Foam::point const &point) {
    position[0] = point[0];
    position[1] = point[1];
  }

  /** \brief Set position. */
  void set_position(Position const &position) { this->position = position; }

  /** \brief Set position. */
  template <typename Position> void set_position(Position const &position) {
    make_position(position);
  }

  /** \brief Set position and cell. */
  template <typename Position>
  void set_position(Position const &position, Cell cell) {
    set_position(position);
    this->cell = cell;
  }
  /** \return Position from 3D point. */
  static Position make_position(Foam::point const &point) {
    return Position{point[0], point[1]};
  }

  /** \brief Make position. */
  template <typename PositionType>
  static Position make_position(PositionType const &position) {
    if (position.size() > 1)
      return {position[0], position[1]};
    else if (position.size() > 0)
      return {position[0], 0.};
    else
      return {0., 0.};
  }
};

/**
   \class State3D_Periodic PTOF/State.h "PTOF/State.h"
   \brief State for 3D positions with periodicity info
   \details State periodicity holds signed number of cells away from unit cell
   along each basis direction.
   \details Defines:
   - \c position
   - \c periodicity
   - \c cell
   - \c time
   - \c mass
   - \c tag
   - \c info
*/
template <typename Info_t = meta::Empty, typename Time_t = meta::Empty,
          typename Mass_t = meta::Empty, typename Tag_t = meta::Empty>
struct State3D_Periodic {
  using Position = Foam::point;
  using Periodicity = std::array<int, 3>;
  using Cell = Foam::label;
  using Info = Info_t;
  using Time = Time_t;
  using Mass = Mass_t;
  using Tag = Tag_t;

  Position position{0., 0., 0.};
  Periodicity periodicity{};
  Cell cell{-1};
  Info info{};
  Time time{};
  Mass mass{};
  Tag tag{};

  /** \brief Set position from 3D point. */
  void set_position(Foam::point const &point) { position = point; }

  /** \brief Set position. */
  template <typename Position> void set_position(Position const &position) {
    make_position(position);
  }

  /** \brief Set position and cell. */
  template <typename Position>
  void set_position(Position const &position, Cell cell) {
    set_position(position);
    this->cell = cell;
  }
  /** \return Position from 3D point. */
  static Position make_position(Foam::point const &point) { return point; }

  /** \brief Make position. */
  template <typename PositionType>
  static Position make_position(PositionType const &position) {
    if (position.size() > 2)
      return {position[0], position[1], position[2]};
    else if (position.size() > 1)
      return {position[0], 0., 0.};
    else
      return {0., 0., 0.};
  }
};

/**
   \return State for 1D positions.
   \details Locate mesh cell using \c Locator object.
*/
template <typename Info, typename Locator, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State1D_Periodic<Info, Time, Mass, Tag>
make_state(Foam::scalar position, int periodicity, Info const &info,
           Locator const &locator, Foam::label hint = -1, Time time = {},
           Mass const &mass = {}, Tag tag = {}) {
  Foam::label cell_id = locator(position, hint);
  if (outside(cell_id))
    throw std::runtime_error{"Particle initialized outside mesh"};
  return {position, periodicity, cell_id, info, time, mass, tag};
}

/**
   \return State for 2D positions.
   \details Locate mesh cell using \c Locator object.
*/
template <typename Info, typename Locator, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State2D_Periodic<Info, Time, Mass, Tag>
make_state(Foam::Vector2D<Foam::scalar> const &position,
           std::array<int, 2> const &periodicity, Info const &info,
           Locator const &locator, Foam::label hint = -1, Time time = {},
           Mass const &mass = {}, Tag tag = {}) {
  Foam::label cell_id = locator(position, hint);
  if (outside(cell_id))
    throw std::runtime_error{"Particle initialized outside mesh"};
  return {position, periodicity, cell_id, info, time, mass, tag};
}

/**
   \return State for 3D positions.
   \details Locate mesh cell using \c Locator object.
*/
template <typename Info, typename Locator, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State3D_Periodic<Info, Time, Mass, Tag>
make_state(Foam::vector const &position, std::array<int, 3> const &periodicity,
           Info const &info, Locator const &locator, Foam::label hint = -1,
           Time time = {}, Mass const &mass = {}, Tag tag = {}) {
  Foam::label cell_id = locator(position, hint);
  if (outside(cell_id))
    throw std::runtime_error{"Particle initialized outside mesh"};
  return {position, periodicity, cell_id, info, time, mass, tag};
}

/**
   \return State for 1D positions with periodicity info
   \details Initial cell index to be provided.
*/
template <typename Info, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State1D_Periodic<Info, Time, Mass, Tag>
make_state(Foam::scalar position, int periodicity, Info const &info,
           Foam::label cell_id, Time time = {}, Mass const &mass = {},
           Tag tag = {}) {
  return {position, periodicity, cell_id, info, time, mass, tag};
}

/**
   \return State for 2D positions with periodicity info
   \details Initial cell index to be provided.
*/
template <typename Info, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State2D_Periodic<Info, Time, Mass, Tag>
make_state(Foam::Vector2D<Foam::scalar> const &position,
           std::array<int, 2> const &periodicity, Info const &info,
           Foam::label cell_id, Time time = {}, Mass const &mass = {},
           Tag tag = {}) {
  return {position, periodicity, cell_id, info, time, mass, tag};
}

/**
   \return State for 3D positions with periodicity info
   \details Initial cell index to be provided.
*/
template <typename Info, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
State3D_Periodic<Info, Time, Mass, Tag>
make_state(Foam::vector const &position, std::array<int, 3> const &periodicity,
           Info const &info, Foam::label cell_id, Time time = {},
           Mass const &mass = {}, Tag tag = {}) {
  return {position, periodicity, cell_id, info, time, mass, tag};
}

/** \brief Periodic state type templated on spatial dimension. */
template <std::size_t dim, typename Info, typename Time = meta::Empty,
          typename Mass = meta::Empty, typename Tag = meta::Empty>
using State_Periodic = typename std::conditional_t<
    dim == 3, State3D_Periodic<Info, Time, Mass, Tag>,
    std::conditional_t<dim == 2, State2D_Periodic<Info, Time, Mass, Tag>,
                       State1D_Periodic<Info, Time, Mass, Tag>>>;
} // namespace ptof

#endif /* PTOF_STATE_H */
