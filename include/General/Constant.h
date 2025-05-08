/**
   \file General/Constant.h
   \author Tomás Aquino
   \date 22/08/2015
   \brief Mathematical and physical constants.
*/

#ifndef GENERAL_CONSTANT_H_
#define GENERAL_CONSTANT_H_

#include <cmath>

/** \namespace ct Mathematical and physical constants. */
namespace cnst {
const double pi = 4. * std::atan(1.); /**> Value of pi */
const double gravity = 9.80665;       /**> Acceleration of gravity in m/s */
const double von_karman = 0.4;        /**> von Karman constant */
} // namespace constants

#endif /* GENERAL_CONSTANT_H_ */
