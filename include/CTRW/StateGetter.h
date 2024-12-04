/**
   \file CTRW/StateGetter.h
   \author Tomás Aquino
   \date 17/05/2020
   \brief Handlers to extract information from particles and particle states.
*/

#ifndef CTRW_STATEGETTER_H
#define CTRW_STATEGETTER_H

#include "General/Operations.h"
#include "General/Useful.h"
#include <type_traits>
#include <utility>

namespace ctrw {
/**
   \class Get_new_from_particle CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get info given particle.
   \tparam Getter Object to extract quantity from particle's new state.
*/
template <typename Getter> struct Get_new_from_particle {
  Getter get;

  Get_new_from_particle(Getter &&get = {}) : get{std::forward<Getter>(get)} {}

  template <typename Particle> auto operator()(Particle const &particle) const {
    return get(particle.state_new());
  }
};
template <typename Getter>
Get_new_from_particle(Getter &&) -> Get_new_from_particle<Getter>;

/**
   \class Get_old_from_particle CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get info given particle.
   \tparam Getter Object to extract quantity from particle's old state.
*/
template <typename Getter> struct Get_old_from_particle {
  Getter get;

  Get_old_from_particle(Getter &&get = {}) : get{std::forward<Getter>(get)} {}

  template <typename Particle> auto operator()(Particle const &particle) const {
    return get(particle.state_old());
  }
};
template <typename Getter>
Get_old_from_particle(Getter &&) -> Get_old_from_particle<Getter>;

/**
   \class Get_from_particle CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get info given particle.
   \tparam Getter Object to extract quantity from particle's new and old state.
*/
template <typename Getter> struct Get_from_particle {
  Getter get;

  Get_from_particle(Getter &&get = {}) : get{std::forward<Getter>(get)} {}
  template <typename Particle> auto operator()(Particle const &particle) const {
    return get(particle.state_new(), particle.state_old());
  }
};
template <typename Getter>
Get_from_particle(Getter &&) -> Get_from_particle<Getter>;

/**
   \struct Get_position CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get position given state.
*/
struct Get_position {
  template <typename State> auto operator()(State const &state) const {
    return state.position;
  }
};

/**
   \class Get_position_component CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get component of position given state.
   \details Component chosen at compile time.
   \tparam dd Component to extract.
*/
template <std::size_t dd> struct Get_position_component {
  template <typename State> auto operator()(State const &state) const {
    return op::project<dd>(state.position);
  }
};

/**
   \struct Get_position_component_arg CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get component of position given state.
   \details Component chosen at runtime when instantiating object.
*/
struct Get_position_component_arg {
  std::size_t dd;

  Get_position_component_arg(std::size_t dd) : dd{dd} {}

  template <typename State> auto operator()(State const &state) const {
    return op::project(state.position, dd);
  }
};

/**
   \struct Get_position_projection CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get projection of position along direction given state.
   \note Direction vector is normalized upon instatiation.
*/
struct Get_position_projection {
  std::vector<double> direction{};

  Get_position_projection(std::vector<double> const &direction)
      : direction{op::div_scalar(direction, op::abs(direction))} {}

  template <typename State> auto operator()(State const &state) const {
    return op::dot(state.position, direction);
  }
};

/**
   \class Get_position_component_sq CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get component squared of state.position given state.
   \tparam dd Component to extract.
*/
template <std::size_t dd> struct Get_position_component_sq {
  template <typename State> auto operator()(State const &state) const {
    double comp = op::project<dd>(state.position);
    return comp * comp;
  }
};

/**
   \struct Get_dist_sq CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get absolute value squared of position given state.
*/
struct Get_dist_sq {
  template <typename State> auto operator()(State const &state) const {
    return op::abs_sq(state.position);
  }
};

/**
   \struct Get_time CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get time given state.
*/
struct Get_time {
  template <typename State> auto operator()(State const &state) const {
    return state.time;
  }
};

/**
   \struct Get_mass CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get mass given state.
*/
struct Get_mass {
  template <typename State> auto operator()(State const &state) const {
    return state.mass;
  }
};

/**
   \struct Get_velocity CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get velocity given state.
*/
struct Get_velocity {
  template <typename State> auto operator()(State const &state) const {
    return state.velocity;
  }
};

/**
   \class Get_velocity_component CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get component of velocity given state.
   \tparam dd Component to extract.
*/
template <std::size_t dd> struct Get_velocity_component {
  template <typename State> auto operator()(State const &state) const {
    return op::project<dd>(state.velocity);
  }
};

/**
   \class Get_property CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get function of state given state.
   \tparam Property Object to return property given state.
*/
template <typename Property> struct Get_property {
  Property property;

  Get_property(Property &&property = {})
      : property{std::forward<Property>(property)} {}

  template <typename State> auto operator()(State const &state) const {
    return property(state);
  }
};
template <typename Property>
Get_property(Property &&) -> Get_property<Property>;

/**
   \class Get_cell_property CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get function of cell given state.
   \tparam Property Object to return property given cell.
*/
template <typename Property> struct Get_cell_property {
  Property property;

  Get_cell_property(Property &&property = {})
      : property{std::forward<Property>(property)} {}

  template <typename State> auto operator()(State const &state) const {
    return property(state.cell);
  }
};
template <typename Property>
Get_cell_property(Property &&) -> Get_cell_property<Property>;

/**
   \class Get_cell_property_arrayop CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get function of cell given state.
   \tparam Property Object to return property given cell through operator[].
*/
template <typename Property> struct Get_cell_arraystyleproperty {
  Property property;

  Get_cell_arraystyleproperty(Property &&property = {})
      : property{std::forward<Property>(property)} {}

  template <typename State> auto operator()(State const &state) const {
    return property[state.cell];
  }
};
template <typename Property>
Get_cell_arraystyleproperty(Property &&)
    -> Get_cell_arraystyleproperty<Property>;

/**
   \class Get_position_property CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get function of position given state.
   \tparam Property Object to return property given state.
*/
template <typename Property> struct Get_position_property {
  Property property;

  Get_position_property(Property &&property = {})
      : property{std::forward<Property>(property)} {}

  template <typename State> auto operator()(State const &state) const {
    return property(state.position);
  }
};
template <typename Property>
Get_position_property(Property &&) -> Get_position_property<Property>;

/**
   \class Get_velocity_property CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get function of velocity given state.
   \tparam Property Object to return property given velocity.
*/
template <typename Property> struct Get_velocity_property {
  Property property;

  Get_velocity_property(Property &&property = {})
      : property{std::forward<Property>(property)} {}

  template <typename State> auto operator()(State const &state) const {
    return property(state.velocity);
  }
};
template <typename Property>
Get_velocity_property(Property &&) -> Get_velocity_property<Property>;

/**
   \struct Get_position_periodic_cartesian CTRW/StateGetter.h
   "CTRW/StateGetter.h"
   \brief Get absolute position given state from position and periodicity.
   \details <tt>state.periodicity</tt> should hold the (directional) integer
   number of times unit cell has been traversed.
   \note Direction vector is normalized upon construction.
*/
struct Get_position_periodic_cartesian {
  std::vector<double> domain_dimensions;

  Get_position_periodic_cartesian(std::vector<double> domain_dimensions)
      : domain_dimensions{domain_dimensions} {}

  template <typename State> auto operator()(State const &state) const {
    return op::plus(state.position,
                    op::times(domain_dimensions, state.periodicity));
  }
};

/**
   \class Get_position_periodic_cartesian_component CTRW/StateGetter.h
   "CTRW/StateGetter.h"
   \brief Get component of absolute position given state from position and
   periodicity
   \details <tt>state.periodicity</tt> should hold the (directional) integer
   number of times unit cell has been traversed.
   \tparam dd Component to extract.
   \note Unit cell faces must be perpendicular to (cartesian) coordinate axes.
*/
template <std::size_t dd> struct Get_position_periodic_cartesian_component {
  Get_position_periodic_cartesian get_position;

  Get_position_periodic_cartesian_component(
      std::vector<double> domain_dimensions)
      : get_position{domain_dimensions} {}

  template <typename State> auto operator()(State const &state) const {
    return op::project<dd>(get_position(state));
  }
};

/**
   \struct Get_position_periodic_cartesian_projection CTRW/StateGetter.h
   "CTRW/StateGetter.h"
   \brief Get projection of absolute position along given direction, given state
   from position and periodicity.
   \details <tt>state.periodicity</tt> should hold the (directional) integer
   number of times unit cell has been traversed.
   \note
   - Unit cell faces must be perpendicular to (cartesian) coordinate axes.
   - Direction vector is normalized upon construction.
*/
struct Get_position_periodic_cartesian_projection {
  Get_position_periodic_cartesian_projection(
      std::vector<double> domain_dimensions,
      std::vector<double> const &direction)
      : _get_position{domain_dimensions}, direction{op::div_scalar(
                                              direction, op::abs(direction))} {}

  template <typename State> auto operator()(State const &state) const {
    return op::dot(_get_position(state), direction);
  }

private:
  Get_position_periodic_cartesian _get_position;

public:
  std::vector<double> direction;
};

/**
   \class Get_position_periodic CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get absolute position given state from position and periodicity
   \tparam Boundary Object that must define \c translate method given position
   and periodicity, which changes the first according to the latter.
*/
template <typename Boundary> struct Get_position_periodic {
  Get_position_periodic(Boundary &&boundary)
      : _boundary{std::forward<Boundary>(boundary)} {}

  template <typename State> auto operator()(State const &state) const {
    auto position = state.position;
    _boundary.translate(position, state.periodicity);
    return position;
  }

private:
  Boundary _boundary;
};
template <typename Boundary>
Get_position_periodic(Boundary &&) -> Get_position_periodic<Boundary>;

/**
   \class Get_projection CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get projection of property along direction.
   \tparam Getter Object to extract quantity to project from state.
   \note Direction vector is normalized upon construction.
*/
template <typename Getter> struct Get_projection {
  Get_projection(Getter &&get, std::vector<double> const &direction)
      : _get{std::forward<Getter>(get)}, direction{op::div_scalar(
                                             direction, op::abs(direction))} {}

  template <typename State> auto operator()(State const &state) const {
    return op::dot(_get(state), direction);
  }

private:
  Getter _get;

public:
  std::vector<double> direction;
};
template <typename Getter>
Get_projection(Getter &&, std::vector<double> const &)
    -> Get_projection<Getter>;

/**
   \struct Get_tag CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get tag given state.
*/
struct Get_tag {
  template <typename State> auto operator()(State const &state) const {
    return state.tag;
  }
};

/**
   \class Get_new CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get new state given new and old state.
*/
template <typename Getter> struct Get_new {
  Getter get;

  Get_new(Getter &&get = {}) : get{std::forward<Getter>(get)} {}

  template <typename State>
  auto operator()(State const &state_new, State const &) const {
    return get(state_new);
  }
};
template <typename Getter> Get_new(Getter &&) -> Get_new<Getter>;

/**
   \class Get_old CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get old state given new and old state.
*/
template <typename Getter> struct Get_old {
  Getter get;

  Get_old(Getter &&get = {}) : get{std::forward<Getter>(get)} {}

  template <typename State>
  auto operator()(State const &, State const &state_old) const {
    return get(state_old);
  }
};
template <typename Getter> Get_old(Getter &&) -> Get_old<Getter>;

/**
   \class Get_interp CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get linearly interpolated quantity given new state and old state.
   \tparam Getter Object to extract quantity given state.
*/
template <typename Getter> struct Get_interp {
  double time;
  Getter get;

  Get_interp(double time, Getter &&get = {})
      : time{time}, get{std::forward<Getter>(get)} {}

  template <typename State>
  auto operator()(State const &state_new, State const &state_old) const {
    auto dt_state = state_new.time - state_old.time;
    auto value_old = get(state_old);
    if (dt_state == 0.)
      return value_old;

    auto delta_t = time - state_old.time;
    auto value_new = get(state_new);
    return op::plus(
        value_old,
        op::times_scalar(delta_t / dt_state, op::minus(value_new, value_old)));
  }
};
template <typename Getter> Get_interp(double, Getter &&) -> Get_interp<Getter>;

/**
   \class Get_position_interp_velocity CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get linearly interpolated position according to velocity in old state
   given new state and old state.
   \tparam Getter Object to extract position given state.
   \tparam Getter_velocity Object to extract velocity given state.
*/
template <typename Getter = Get_position,
          typename Getter_velocity = Get_velocity>
struct Get_position_interp_velocity {
  double time;
  Getter get;
  Getter_velocity get_velocity;

  Get_position_interp_velocity(double time, Getter &&get = {},
                               Getter_velocity &&get_velocity = {})
      : time{time}, get{std::forward<Getter>(get)}, get_velocity{get_velocity} {
  }

  template <typename State>
  auto operator()(State const &state_new, State const &state_old) const {
    double delta_t = time - state_old.time;
    auto vel = get_velocity(state_old);
    return op::plus(get(state_old), vel * delta_t);
  }
};
template <typename Getter, typename Getter_velocity>
Get_position_interp_velocity(double, Getter &&, Getter_velocity &&)
    -> Get_position_interp_velocity<Getter, Getter_velocity>;

/**
   \class Get_time_interp_velocity CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get linearly interpolated time according to velocity in old state
   given new state and old state.
   \tparam Getter_velocity Object to extract velocity given state.
   \tparam Getter_position Object to extract position given state.
*/
template <typename Getter_velocity = Get_velocity,
          typename Getter_position = Get_position>
struct Get_time_interp_velocity {
  double position;
  Getter_velocity get_velocity;
  Getter_position get_position;

  Get_time_interp_velocity(double position, Getter_velocity &&get_velocity = {},
                           Getter_position &&get_position = {})
      : position{position}, get_velocity{std::forward<Getter_velocity>(
                                get_velocity)},
        get_position{std::forward<Getter_position>(get_position)} {}

  template <typename State>
  auto operator()(State const &state_new, State const &state_old) const {
    double delta_x = std::abs(position - getter_position(state_old));
    double vel = op::abs(get_velocity(state_old.velocity));
    return state_old.time + delta_x / vel;
  }
};
template <typename Getter_velocity, typename Getter_position>
Get_time_interp_velocity(double, Getter_velocity &&, Getter_position &&)
    -> Get_time_interp_velocity<Getter_velocity, Getter_position>;

/**
   \class Get_time_interp CTRW/StateGetter.h "CTRW/StateGetter.h"
   \brief Get linearly interpolated time.
   \tparam Getter_position Object to extract position given state.
*/
template <typename Getter_position = Get_position> struct Get_time_interp {
  double position;
  Getter_position get_position;

  Get_time_interp(double position, Getter_position &&get_position = {})
      : position{position}, get_position{
                                std::forward<Getter_position>(get_position)} {}

  template <typename State>
  auto operator()(State const &state_new, State const &state_old) const {
    auto dt_state = state_new.time - state_old.time;
    if (dt_state == 0.)
      return state_old.time;
    double delta_x = std::abs(position - get_position(state_old));
    double vel =
        std::abs(get_position(state_new) - get_position(state_old)) / dt_state;
    return state_old.time + delta_x / vel;
  }
};
template <typename Getter_position>
Get_time_interp(double, Getter_position &&) -> Get_time_interp<Getter_position>;
} // namespace ctrw

#endif /* CTRW_STATEGETTER_H */
