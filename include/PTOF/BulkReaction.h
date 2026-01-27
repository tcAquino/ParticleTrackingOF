/**
   \file PTOF/BulkReaction.h
   \author Tomas Aquino
   \date 10/03/2022
   \brief Objects and utilities for surface and bulk reactions.
*/

#ifndef PTOF_BULKREACTION_H
#define PTOF_BULKREACTION_H

#include "General/IO.h"
#include <cmath>
#include <ostream>

namespace ptof {
/**
   \class BulkReaction_DoNothing PTOF/BulkReaction.h "PTOF/BulkReaction.h"
   \brief Bulk reaction that does nothing.
*/
class BulkReaction_DoNothing {
public:
  /**
     \brief React over exposure time (do nothing).
     \param state State of particle to react (unused).
     \param exposure_time Exposure time over which to react (unused).
  */
  template <typename State>
  void operator()(State &state, double exposure_time) const {}

  /**
     \brief Output generic information about object.
     \param output Output stream.
  */
  inline static std::ostream &info(std::ostream &output) {
    return output << io::line() << "Bulk reaction\n"
                  << io::line() << "None\n"
                  << io::line();
  }
};
} // namespace ptof

#endif /* PTOF_BULKREACTION_H */
