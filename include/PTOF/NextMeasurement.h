/**
   \file PTOF/NextMeasurement.h
   \author Tomás Aquino
   \date 29/09/2024
   \brief Handle setting up next measurement.
*/

#ifndef PTOF_NEXTMEASUREMENT_H
#define PTOF_NEXTMEASUREMENT_H

#include <cmath>
#include <cstddef>
#include <limits>

/**
   \class NextMeasurementTime PTOF/NextMeasurement.h "PTOF/NextMeasurement.h"
   \brief Abstract polymorphic object to handle setting up the next measurement
   time.
*/
template <typename Parameters> struct NextMeasurementTime {
  /** Destructor. */
  virtual ~NextMeasurementTime() = default;

  /** \return Next measurement time. */
  double time() const { return _next_time; }

  /** \brief Set up for next measurement.  */
  void advance() {
    ++_next_measurement;
    set_next_time();
  }

protected:
  /** \brief Constructor. */
  NextMeasurementTime(Parameters const &parameters)
      : _parameters{parameters}, _next_time{parameters.time_min} {}

  /** \brief Deleted move constructor. */
  NextMeasurementTime(Parameters &&) = delete;

  Parameters const &_parameters;    /**< Output parameters.         */
  std::size_t _next_measurement{0}; /**< Index of next measurement. */
  double _next_time;                /**< Time of next measurement.  */

  /** \brief Set the time of the next measurment. */
  virtual void set_next_time() = 0;
};

/**
   \class NextMeasurementTime_linear_step PTOF/Output_Cases.h
   "PTOF/Output_Cases.h" \brief NextMeasurementTime object for constant time
   step.
*/
template <typename Parameters>
struct NextMeasurementTime_linear_step final : NextMeasurementTime<Parameters> {
  NextMeasurementTime_linear_step(Parameters const &parameters)
      : NextMeasurementTime<Parameters>{parameters} {}

  void set_next_time() override {
    NextMeasurementTime<Parameters>::_next_time =
        NextMeasurementTime<Parameters>::_parameters.time_min +
        NextMeasurementTime<Parameters>::_next_measurement *
            NextMeasurementTime<Parameters>::_parameters.time_increment;
  }
};

/**
   \class NextMeasurementTime_log_step PTOF/Output_Cases.h "PTOF/Output_Cases.h"
   \brief NextMeasurementTime object for constant time step.
*/
template <typename Parameters>
struct NextMeasurementTime_log_step final : NextMeasurementTime<Parameters> {
  NextMeasurementTime_log_step(Parameters const &parameters)
      : NextMeasurementTime<Parameters>{parameters} {}

  void set_next_time() override {
    NextMeasurementTime<Parameters>::_next_time *=
        NextMeasurementTime<Parameters>::_parameters.time_increment;
  }
};

/**
   \class NextMeasurementTime_linear PTOF/Output_Cases.h "PTOF/Output_Cases.h"
   \brief NextMeasurementTime object for linearly-spaced times.
*/
template <typename Parameters>
struct NextMeasurementTime_linear final : NextMeasurementTime<Parameters> {
  NextMeasurementTime_linear(Parameters const &parameters)
      : NextMeasurementTime<Parameters>{parameters} {}

  void set_next_time() override {
    NextMeasurementTime<Parameters>::_next_time =
        NextMeasurementTime<Parameters>::_next_time >
                NextMeasurementTime<Parameters>::_parameters.time_max
            ? std::numeric_limits<double>::infinity()
            : NextMeasurementTime<Parameters>::_parameters.time_min +
                  NextMeasurementTime<Parameters>::_next_measurement *
                      NextMeasurementTime<Parameters>::_parameters
                          .time_increment;
    ;
  }
};

/**
   \struct Output_Cases::NextMeasurementTime_log PTOF/Output_Cases.h
   "PTOF/Output_Cases.h"
   \brief NextMeasurementTime object for log-spaced times.
*/
template <typename Parameters>
struct NextMeasurementTime_log final : NextMeasurementTime<Parameters> {
  NextMeasurementTime_log(Parameters const &parameters)
      : NextMeasurementTime<Parameters>{parameters} {}

  void set_next_time() override {
    NextMeasurementTime<Parameters>::_next_time =
        NextMeasurementTime<Parameters>::_next_time >
                NextMeasurementTime<Parameters>::_parameters.time_max
            ? std::numeric_limits<double>::infinity()
            : NextMeasurementTime<Parameters>::_parameters.time_min *
                  std::pow(NextMeasurementTime<Parameters>::_parameters
                               .time_increment,
                           NextMeasurementTime<Parameters>::_next_measurement);
  }
};

#endif /* PTOF_NEXTMEASUREMENT_H */
