/**
   \file PTOF/CheckOptions.h
   \author Tomas Aquino
   \date 2022/02/17
   \brief Helpers for choosing checking and warning level.
*/

#ifndef PTOF_CHECKOPTIONS_H
#define PTOF_CHECKOPTIONS_H

namespace ptof {
/**
   \struct CheckOptions PTOF/CheckOptions.h "PTOF/CheckOptions.h"
   \brief Options for checking conditions.
*/
struct CheckOptions {
  /**
     \struct CheckOptions::NoCheck PTOF/CheckOptions.h "PTOF/CheckOptions.h"
     \brief No bounds checking.
  */
  struct NoCheck {};

  /**
     \struct CheckOptions::Check PTOF/CheckOptions.h "PTOF/CheckOptions.h"
     \brief Bounds checking, no warning.
  */
  struct Check {};

  /**
     \struct CheckOptions::Warn PTOF/CheckOptions.h "PTOF/CheckOptions.h"
     \brief Bounds checking and warning.
  */
  struct Warn {};
};
} // namespace ptof

#endif /* PTOF_CHECKOPTIONS_H */
