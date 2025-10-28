/**
   \file PTOF/Roots.h
   \author Tomas Aquino
   \date 09/03/2022
   \brief Root-finding algorithms.
*/

#ifndef PTOF_ROOTS_H
#define PTOF_ROOTS_H

#include <cmath>
#include <cstddef>

namespace roots {
class NewtonRaphson {
public:
  template <typename Func, typename Derivative>
  double operator()(Func func, Derivative derivative, double guess,
                    double tolerance, std::size_t max_iter) {
    _func_val = func(guess);
    _converged = false;
    for (std::size_t iter = 0; iter < max_iter && !_converged; ++iter) {
      guess = guess - _func_val / derivative(guess);
      _func_val = func(guess);
      _converged = std::abs(_func_val) < tolerance;
    }
    return guess;
  }

  bool converged() { return _converged; }
  double func_val() { return _func_val; }

private:
  bool _converged;
  double _func_val;
};

class Bisection {
public:
  template <typename Func>
  double operator()(Func func, double lower_bracket, double upper_bracket,
                    double tolerance) {
    int nr_iterations =
        std::log((upper_bracket - lower_bracket) / tolerance) / std::log(2.) -
        1;
    double mid_point = (lower_bracket + upper_bracket) / 2.;
    for (int count = 0; count < nr_iterations; ++count) {
      if (func(lower_bracket) * func(mid_point) < 0.) {
        upper_bracket = mid_point;
      } else {
        lower_bracket = mid_point;
      }
      mid_point = (lower_bracket + upper_bracket) / 2.;
    }
    _func_val = func(mid_point);
    _converged = std::abs(_func_val) < tolerance;
    return mid_point;
  }

  bool converged() { return _converged; }
  double func_val() { return _func_val; }

private:
  bool _converged;
  double _func_val;
};
} // namespace roots

#endif /* PTOF_ROOTS_H */
