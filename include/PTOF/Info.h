/**
   \file PTOF/Info.h
   \author Tomás Aquino
   \date 03/09/2022
   \brief Objects and utilities to store information in particle states.
*/

#ifndef PTOF_INFO_H
#define PTOF_INFO_H

#include "General/Useful.h"
#include "PTOF/Meta.h"
#include <cstddef>
#include <fieldTypes.H>
#include <point.H>
#include <string>

namespace ptof {
/**
   \struct Info_Type PTOF/Info.h "PTOF/Info.h"
   \brief Information about type.
*/
struct Info_Type {
  std::string type;
};

/**
   \struct Info_Absorbed PTOF/Info.h "PTOF/Info.h"
   \brief Information about absorption.
*/
struct Info_Absorbed {
  bool absorbed;
};

/**
   \struct Info_Absorbed PTOF/Info.h "PTOF/Info.h"
   \brief Information about absorption and boundary patch where absorption
   happened.
*/
struct Info_Absorbed_Patch {
  bool absorbed;
  Foam::label patch_id;
};

/**
   \struct Info_Absorbed_Reinjections PTOF/Info.h "PTOF/Info.h"
   \brief Information about absorption and number of reinjections.
*/
struct Info_Absorbed_Reinjections {
  bool absorbed;
  std::size_t reinjections;
};

/**
   \struct Info_Absorbed_Contact_Point PTOF/Info.h "PTOF/Info.h"
   \brief Information about contact point.
*/
struct Info_Contact_Point {
  Foam::point contact_point;
};
} // namespace ptof

#endif /* PTOF_INFO_H */
