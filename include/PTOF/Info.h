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
   \struct Info_absorbed PTOF/Info.h "PTOF/Info.h"
   \brief Information about absorption.
*/
struct Info_absorbed {
  bool absorbed;
};

/**
   \struct Info_absorbed_adsorbed PTOF/Info.h "PTOF/Info.h"
   \brief Information about absorption and adsorption.
*/
struct Info_absorbed_adsorbed {
  bool absorbed;
  bool adsorbed;
};

/**
   \struct Info_absorbed_face PTOF/Info.h "PTOF/Info.h"
   \brief Information about absorption and boundary face where absorption
   happened.
*/
struct Info_absorbed_face {
  bool absorbed;
  Foam::label face;
};

/**
   \struct Info_absorbed_adsorbed_face PTOF/Info.h "PTOF/Info.h"
   \brief Information about absorption and boundary face where absorption
   happened.
*/
struct Info_absorbed_adsorbed_face {
  bool absorbed;
  bool adsorbed;
  Foam::label face;
};

/**
   \struct Info_absorbed_reinjections PTOF/Info.h "PTOF/Info.h"
   \brief Information about absorption and number of reinjections.
*/
struct Info_absorbed_reinjections {
  bool absorbed;
  std::size_t reinjections;
};

/**
   \struct Info_absorbed_adsorbed_reinjections PTOF/Info.h "PTOF/Info.h"
   \brief Information about absorption, adsorption, and number of reinjections.
*/
struct Info_absorbed_adsorbed_reinjections {
  bool absorbed;
  bool adsorbed;
  std::size_t reinjections;
};

/**
   \struct Info_absorbed_contact_point PTOF/Info.h "PTOF/Info.h"
   \brief Information about absorption and contact point.
*/
struct Info_absorbed_contact_point {
  bool absorbed;
  Foam::point contact_point;
};

/**
   \struct Info_absorbed_adsorbed_contact_point PTOF/Info.h "PTOF/Info.h"
   \brief Information about absorption, adsorption, and contact point.
*/
struct Info_absorbed_adsorbed_contact_point {
  bool absorbed;
  bool adsorbed;
  Foam::point contact_point;
};
} // namespace ptof

#endif /* PTOF_INFO_H */
