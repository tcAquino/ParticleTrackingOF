/**
 \file General/Constants.h
 \author Tomás Aquino
 \date 22/08/2015
*/

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <cmath>

/** \namespace constants Values of physical constants. */
namespace constants {
const double pi = 4. * std::atan(1.); /**> Value of pi */
const double gravity = 9.80665;       /**> Acceleration of gravity in m/s */
const double von_karman = 0.4;        /**> von Karman constant */
} // namespace constants

#endif /* CONSTANTS_H_ */
