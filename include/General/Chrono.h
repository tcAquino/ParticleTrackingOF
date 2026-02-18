/**
   \file General/Chrono.h
   \author Tomas Aquino
   \date 18/02/2026
   \brief Utilities related to clock time.
*/

#ifndef GENERAL_CHRONO_H
#define GENERAL_CHRONO_H

#include "General/IO.h"
#include <chrono>
#include <cstdint>
#include <iomanip>

namespace chrono {
// Adapted from Howard Hinnant's answer here:
// https://stackoverflow.com/questions/22590821/convert-stdduration-to-human-readable-time
/** \brief Display execution time in human-readable format. */
inline std::ostream &display_duration(std::ostream &stream,
                                      std::chrono::nanoseconds ns) {
  io::StreamScopeFormat guard{stream};
  using days = std::chrono::duration<std::int_least64_t, std::ratio<86400>>;
  stream.fill('0');
  stream << std::right;
  auto dd = std::chrono::duration_cast<days>(ns);
  ns -= dd;
  auto hh = std::chrono::duration_cast<std::chrono::hours>(ns);
  ns -= hh;
  auto mm = std::chrono::duration_cast<std::chrono::minutes>(ns);
  ns -= mm;
  auto ss = std::chrono::duration_cast<std::chrono::seconds>(ns);
  ns -= ss;
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(ns);
  stream << std::setw(2) << dd.count() << "d:" << std::setw(2) << hh.count()
         << "h:" << std::setw(2) << mm.count() << "m:" << std::setw(2)
         << ss.count() << "s:" << std::setw(3) << ms.count() << "ms";

  return stream;
};

/** \brief Display execution time in human-readable format. */
template <typename Clock>
std::ostream &display_duration(std::ostream &stream,
                               std::chrono::time_point<Clock> start_time,
                               std::chrono::time_point<Clock> end_time) {
  return display_duration(stream,
                          std::chrono::duration_cast<std::chrono::nanoseconds>(
                              end_time - start_time));
}
} // namespace chrono

#endif /* GENERAL_CHRONO_H */
