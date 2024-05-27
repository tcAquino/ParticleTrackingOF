/**
 \file General/Serial.h
 \author Tomás Aquino
 \date 14/05/2024
 \brief Algorithms for serial implementations
*/

#ifndef GENERAL_SERIAL_H
#define GENERAL_SERIAL_H

#include "General/Meta.h"

namespace useful {
inline std::size_t get_thread_num(meta::ParallelOptions::Serial) { return 0; }

inline std::size_t get_num_threads(meta::ParallelOptions::Serial) { return 1; }
} // namespace useful

#endif /* GENERAL_SERIAL_H */
