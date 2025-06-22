/**
   \file CTRW/State.h
   \author Tomas Aquino
   \date 09/26/2017
   \brief Particle states.
*/

#ifndef CTRW_STATE_H
#define CTRW_STATE_H

#include "General/Meta.h"
#include "General/Operation.h"
#include <valarray>
#include <vector>

namespace ctrw {
/**
   \class State_Periodic CTRW/State.h "CTRW/State.h"
   \brief State for positions including information about periodicity.
   \details
   - Defines:
       - \c position
       - \c periodicity
       - \c mass
       - \c time
       - \c tag.
   - Periodicity information keeps track of the periodic cell the position is
   in.
*/
template <typename Position_t, typename Periodicity_t = std::vector<int>,
          typename Mass_t = double, typename Time_t = double,
          typename Tag_t = meta::Empty>
struct State_Periodic {
  using Position = Position_t;
  using Periodicity = Periodicity_t;
  using Mass = Mass_t;
  using Time = Time_t;
  using Tag = Tag_t;

  State_Periodic() {}

  State_Periodic(Position position, Periodicity periodicity, Mass mass = 1,
                 Time_t time = 0, Tag_t tag = {})
      : position{position},
        periodicity{periodicity}, mass{mass}, time{time}, tag{tag} {}

  Position position;
  Periodicity periodicity;
  Mass mass;
  Time time;
  Tag tag;
};

/**
   \class State_Velocity CTRW/State.h "CTRW/State.h"
   \details Defines:
   - \c position
   - \c velocity
   - \c time
   - \c tag
*/
template <typename Position_t, typename Velocity_t, typename Time_t = double,
          typename Tag_t = meta::Empty>
struct State_Velocity {
  using Position = Position_t;
  using Velocity = Velocity_t;
  using Time = Time_t;
  using Tag = Tag_t;

  State_Velocity() {}

  State_Velocity(Position_t position, Velocity_t velocity, Time_t time = {},
                 Tag_t tag = {})
      : position(position), velocity(velocity), time(time), tag(tag) {}

  Position position;
  Velocity velocity;
  Time time;
  Tag tag;
};

/**
   \class State_Velocity_Mass CTRW/State.h "CTRW/State.h"
   \details Defines:
   - \c position
   - \c velocity
   - \c mass
   - \c time
   - \c tag
*/
template <typename Position_t, typename Velocity_t, typename Mass_t = double,
          typename Time_t = double, typename Tag_t = meta::Empty>
struct State_Velocity_Mass {
  using Position = Position_t;
  using Velocity = Velocity_t;
  using Mass = Mass_t;
  using Time = Time_t;
  using Tag = Tag_t;

  State_Velocity_Mass() {}

  State_Velocity_Mass(Position position, Velocity velocity, Mass mass = {},
                      Time time = {}, Tag tag = {})
      : position{position}, velocity{velocity}, mass{mass}, time{time},
        tag{tag} {}

  Position position;
  Velocity velocity;
  Mass mass;
  Time time;
  Tag tag;
};

/**
   \class State_Position CTRW/State.h "CTRW/State.h"
   \details Defines:
   - \c position
   - \c time
   - \c tag
*/
template <typename Position_t, typename Time_t = double,
          typename Tag_t = meta::Empty>
struct State_position {
  using Position = Position_t;
  using Time = Time_t;
  using Tag = Tag_t;

  State_position() {}

  State_position(Position position, Time time = {}, Tag tag = {})
      : position(position), time(time) {}

  Position position{};
  Time time{};
  Tag tag{};
};

/**
   \class State_Position_Mass CTRW/State.h "CTRW/State.h"
   \details Defines:
   - \c position
   - \c mass
   - \c time
   - \c tag
*/
template <typename Position_t, typename Mass_t = double,
          typename Time_t = double, typename Tag_t = meta::Empty>
struct State_Position_Mass {
  using Position = Position_t;
  using Mass = Mass_t;
  using Time = Time_t;
  using Tag = Tag_t;

  State_Position_Mass() {}

  State_Position_Mass(Position position, Mass mass = 1., Time time = {},
                      Tag tag = {})
      : position{position}, mass{mass}, time{time}, tag{tag} {}

  Position position;
  Mass mass;
  Time time;
  Tag tag;
};

/**
   \class State_PTRW_Position CTRW/State.h "CTRW/State.h"
   \details Defines:
   - \c position
*/
template <typename Position_t> struct State_PTRW_Position {
  using Position = Position_t;

  State_PTRW_Position() {}

  State_PTRW_Position(Position position) : position(position) {}

  Position position{};
};

/**
   \class State_RunTumble CTRW/State.h "CTRW/State.h"
   \details Defines:
   - \c position
   - \c orientation
   - \c run
   - \c tag
*/
template <typename Orientation_t, typename Tag_t = meta::Empty>
struct State_RunTumble {
  using Position = std::vector<double>;
  using Orientation = Orientation_t;
  using Time = double;
  using Tag = Tag_t;

  State_RunTumble() {}

  State_RunTumble(Position position, Orientation orientation, bool run = 0,
                  Time time = {}, Tag tag = {})
      : position{position},
        orientation{orientation}, run{run}, time{time}, tag{tag} {}

  Position position{};
  Orientation orientation{};
  bool run;
  Time time{};
  Tag tag;
};

/**
   \class State_RunTumble_PTRW CTRW/State.h "CTRW/State.h"
   \details Defines:
   - \c position
   - \c orientation
   - \c state
   - \c tag
*/
template <typename Orientation_t = double, typename Tag_t = std::size_t,
          typename Position_t = std::vector<double>>
struct State_RunTumble_PTRW {
  using Position = Position_t;
  using Orientation = Orientation_t;
  using Tag = Tag_t;

  State_RunTumble_PTRW() {}

  State_RunTumble_PTRW(Position position, Orientation orientation, int state,
                       Tag tag = {})
      : position{position}, orientation{orientation}, state{state}, tag{tag} {}

  Position position{};
  Orientation orientation{};
  int state;
  Tag tag;
};
} // namespace ctrw

#endif /* CTRW_STATE_H */
