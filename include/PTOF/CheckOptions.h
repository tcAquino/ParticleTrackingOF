/**
 * @file   CheckOptions.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Thu Feb 17 00:00:00 2022
 *
 * @brief Helpers for choosing checking and warning level.
 */

#ifndef PTOF_CHECKOPTIONS_H
#define PTOF_CHECKOPTIONS_H

namespace ptof {
/** @brief Options for checking conditions. */
struct CheckOptions {
  /** @brief No bounds checking. */
  struct NoCheck {};

  /** @brief Bounds checking, no warning. */
  struct Check {};

  /** @brief Bounds checking and warning. */
  struct Warn {};
};
} // namespace ptof

#endif /* PTOF_CHECKOPTIONS_H */
