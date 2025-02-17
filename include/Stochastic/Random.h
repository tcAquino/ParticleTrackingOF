/**
   \file Stochastic/Random.h
   \author Tomás Aquino
   \date 02/07/2015
   \brief Probability and random number generation.
*/

#ifndef STOCHASTIC_RANDOM_H
#define STOCHASTIC_RANDOM_H

#include "Algebra/Algebra.h"
#include "General/Constants.h"
#include "General/IO.h"
#include "General/Operations.h"
#include "General/Parallel.h"
#include "General/Ranges.h"
#include "General/Useful.h"
#include "Geometry/Coordinates.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <list>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

/**
   \namespace ctrw Objects and methods related to probability and random number
   generation.
*/
namespace stochastic {
/**
   \class RNGDist Stochastic/Random.h "Stochastic/Random.h"
   \brief Wrapper for random number generation with own RNG engine.
*/
template <typename Distribution_t, typename Engine_t = std::mt19937>
struct RNGDist {
  using param_type = typename Distribution_t::param_type;
  using result_type = typename Distribution_t::result_type;
  using distribution_type = Distribution_t;
  using Engine = Engine_t;

  template <typename param_t>
  RNGDist(param_t const &params) : dist{param_type(params)} {}

  RNGDist(Distribution_t dist) : dist{dist} {}

  RNGDist(std::size_t seed) : _rng{seed} {}

  result_type operator()() { return dist(_rng); }

  Distribution_t dist;

private:
  Engine_t _rng{std::random_device{}()};
};

/**
   \class RNGDist_shared_engine Stochastic/Random.h "Stochastic/Random.h"
   \brief Wrapper for random number generation with shared RNG engine.
*/
template <typename Distribution_t, typename Engine_t = std::mt19937>
struct RNGDist_shared_engine {
  using param_type = typename Distribution_t::param_type;
  using result_type = typename Distribution_t::result_type;
  using distribution_type = Distribution_t;
  using Engine = Engine_t;

  template <typename param_t>
  RNGDist_shared_engine(param_t const &params, Engine_t &rng)
      : dist{param_type{params}}, _rng{rng} {}

  RNGDist_shared_engine(Distribution_t dist, Engine_t &rng)
      : dist{dist}, _rng{rng} {}

  result_type operator()() { return dist(_rng); }

  Distribution_t dist;

private:
  Engine_t &_rng;
};

/**
   \class RNGThreaded Stochastic/Random.h "Stochastic/Random.h"
   \brief Set of RNGs for thread-safe random number generation.
*/
template <typename ParallelOption, typename Engine_t = std::mt19937>
struct RNGThreaded {
  using Engine = Engine_t;
  using result_type = typename Engine::result_type;

  RNGThreaded() {
    std::size_t num_threads = get_num_threads(ParallelOption{});
    _rngs.reserve(num_threads);
    for (std::size_t ii = 0; ii < num_threads; ++ii)
      _rngs.emplace_back(std::random_device{}());
  }

  /** \brief Constructor.
      \details Ignored seed constructor to ensure compatibility with standard
      RNGs. */
  RNGThreaded(std::size_t) : RNGThreaded{} {}

  /** \brief Constructor.
      \details Seed each RNG separately. */
  RNGThreaded(std::vector<std::size_t> const &seeds) : RNGThreaded{} {
    std::size_t num_threads = get_num_threads(ParallelOption{});
    _rngs.reserve(num_threads);
    for (std::size_t ii = 0; ii < num_threads; ++ii)
      _rngs.emplace_back(seeds[ii]);
  }

  auto operator()() { return _rngs[get_thread_num(ParallelOption{})](); }

  static constexpr auto min() { return Engine_t::min(); }

  static constexpr auto max() { return Engine_t::max(); }

private:
  std::vector<Engine_t> _rngs;
};

template <typename Engine_t>
struct RNGThreaded<par::ParallelOptions::Serial, Engine_t> {
  using Engine = Engine_t;

  RNGThreaded() : _rng{std::random_device{}()} {}

  RNGThreaded(std::size_t seed) : _rng{seed} {}

  auto operator()() { return _rng(); }

  auto min() { return Engine_t::min(); }

  auto max() { return Engine_t::max(); }

private:
  Engine_t _rng;
};

/**
   \class skewedlevystable_distribution Stochastic/Random.h
   "Stochastic/Random.h"
   \brief Skewed Lévy stable distribution.
*/
template <typename Value_type = double> class skewedlevystable_distribution {
public:
  using param_type = std::array<Value_type, 3>;
  using result_type = Value_type;

  const Value_type alpha;
  const Value_type sigma;
  const Value_type mu;

  const double zeta{std::tan(constants::pi * alpha * 0.5)};
  const double xi{std::atan(zeta) / alpha};
  const double vv{std::pow(1. + zeta * zeta, 0.5 / alpha)};

  /**
     \brief Constructor.
     \param alpha Exponent, pdf <tt>x^{-1-alpha}</tt>.
     \param sigma Scale parameter.
     \param mu location parameter.
     \note Tail <tt>~ -sigma^alpha * tan(alpha*pi/2) *
     x^{-1-alpha}/Gamma(-alpha)</tt>.
  */
  skewedlevystable_distribution(Value_type alpha, Value_type sigma = 1.,
                                Value_type mu = 0.)
      : alpha(alpha), sigma(sigma), mu(mu) {}

  /**
     \brief Constructor.
     \param params <tt>param[0] = alpha</tt>, <tt>param[1] = sigma</tt>,
     <tt>param[2] = mu</tt>.
     \note
      - alpha: Tailing exponent, <tt>PDF ~ x^{-1-alpha}</tt>.
      - sigma: Scale parameter.
      - mu: Location parameter.
      - Tail <tt>~ -sigma^alpha * tan(alpha*pi/2) *
      t^{-1-alpha}/Gamma(-alpha)</tt>.
  */
  skewedlevystable_distribution(param_type const &params)
      : alpha(params[0]), sigma(params[1]), mu(params[2]) {}

  template <typename Generator> Value_type operator()(Generator &rng) {
    double uu = constants::pi * (_uniform_dist(rng) - 0.5);
    double tt =
        std::sin(alpha * (uu + xi)) / std::pow(std::cos(uu), 1. / alpha);
    double ss = std::pow(std::cos((1. - alpha) * uu - alpha * xi) /
                             _exponential_dist(rng),
                         (1. - alpha) / alpha);

    return sigma * vv * tt * ss + mu;
  }

private:
  std::uniform_real_distribution<double> _uniform_dist{0., 1.};
  std::exponential_distribution<double> _exponential_dist{1.};
};

/**
   \class pareto_distribution Stochastic/Random.h "Stochastic/Random.h"
   \brief Pareto distribution.
   \note PDF <tt>(x/min)^{-alpha-1}</tt> for <tt>x >= min</tt>.
*/
template <typename Value_type = double> class pareto_distribution {
public:
  using param_type = std::array<Value_type, 2>;
  using result_type = Value_type;

  const Value_type alpha;
  const Value_type min;

  /**
     \brief Constructor.
     \param alpha Tailing exponent, PDF <tt>x^{-1-alpha}</tt>.
     \param min Minimum value.
  */
  pareto_distribution(Value_type alpha, Value_type min)
      : alpha(alpha), min(min) {}

  /**
     \brief Constructor.
     \param params <tt>param[0] = alpha</tt>, <tt>param[1] = min</tt>.
     \note
      - alpha: Tailing exponent, <tt>PDF ~ x^{-1-alpha}</tt>.
      - min: Minimum value.
  */
  pareto_distribution(param_type const &params)
      : alpha(params[0]), min(params[1]) {}

  template <typename Generator> Value_type operator()(Generator &rng) {
    return min * std::exp(_exponential_dist(rng) / alpha);
  }

private:
  std::exponential_distribution<double> _exponential_dist{1.};
};

/**
   \class lomax_distribution Stochastic/Random.h "Stochastic/Random.h"
   \brief Pareto lomax distribution
   \note PDF <tt>alpha/scale (1+x/scale)^{-alpha-1}</tt>. */
template <typename Value_type = double> class lomax_distribution {
public:
  using param_type = std::array<Value_type, 2>;
  using result_type = Value_type;

  const Value_type alpha;
  const Value_type scale;

  /**
     \brief Constructor.
     \param alpha Tailing exponent, PDF <tt>x^{-1-alpha}</tt>.
     \param scale Scale parameter.
  */
  lomax_distribution(Value_type alpha, Value_type scale)
      : alpha{alpha}, scale{scale} {}

  /**
     \brief Constructor.
     \param params <tt>param[0] = alpha</tt>, <tt>param[1] = scale</tt>.
     \note
      - alpha: Tailing exponent, <tt>PDF ~ x^{-1-alpha}</tt>.
      - scale: Scale parameter.
  */
  lomax_distribution(param_type const &params)
      : alpha(params[0]), scale(params[1]) {}

  template <typename Generator> Value_type operator()(Generator &rng) {
    return scale * (std::pow(1. - _uniform_dist(rng), -1. / alpha) - 1.);
  }

private:
  std::uniform_real_distribution<double> _uniform_dist{};
};

/**
   \class inverse_gaussian_distribution Stochastic/Random.h
   "Stochastic/Random.h"
   \brief Inverse Gaussian distribution.
*/
template <typename Value_type = double> class inverse_gaussian_distribution {
public:
  using param_type = std::array<Value_type, 2>;
  using result_type = Value_type;

  const Value_type mean;
  const Value_type var;

  inverse_gaussian_distribution(Value_type mean, Value_type var)
      : mean(mean), var(var) {}

  inverse_gaussian_distribution(param_type const &params)
      : mean(params[0]), var(params[1]) {}

  template <typename Generator> Value_type operator()(Generator &rng) {
    double rand_normal = _normal_dist(rng);
    double aux = _ratio * rand_normal * rand_normal;

    double xx = 1. + aux - std::sqrt(aux * (aux + 2.));
    double rand_uniform = _uniform_dist(rng);

    return 1. / (1. + xx) < rand_uniform ? mean / xx : mean * xx;
  }

private:
  std::uniform_real_distribution<double> _uniform_dist{0., 1.};
  std::normal_distribution<double> _normal_dist{0., 1.};
  double _ratio{0.5 * var / (mean * mean)};
};

/**
   \class isotropic_unit_vector_distribution Stochastic/Random.h
   "Stochastic/Random.h"
   \brief Distribution for unit vectors in arbitrary dimension with uniformly
   random orientation.
*/
class isotropic_unit_vector_distribution {
public:
  using param_type = std::size_t;
  using result_type = std::vector<double>;

  isotropic_unit_vector_distribution(std::size_t dim) : _dim{dim} {}

  template <typename Generator> result_type operator()(Generator &rng) {
    result_type val;
    for (std::size_t dd = 0; dd < _dim; ++dd)
      val.push_back(_normal_dist(rng));
    op::normalize_inplace(val);
    return val;
  }

private:
  param_type _dim;
  std::normal_distribution<double> _normal_dist{};
};

/**
   \class isotropic_unit_vector_distribution Stochastic/Random.h
   "Stochastic/Random.h"
   \brief Uniform distribution inside hypersphere.
*/
class uniform_sphere_distribution {
public:
  struct param_type {
    std::vector<double> center;
    double radius;
    double inner_radius{0.};
  } const params;
  using result_type = std::vector<double>;

  uniform_sphere_distribution(param_type params) : params{params} {}

  uniform_sphere_distribution(std::vector<double> const &center, double radius,
                              double inner_radius = 0.)
      : uniform_sphere_distribution{param_type{center, radius, inner_radius}} {}

  template <typename Generator> result_type operator()(Generator &rng) {
    return op::plus(params.center,
                    op::times_scalar(_uniform_dist(rng), _isotropic_dist(rng)));
  }

private:
  isotropic_unit_vector_distribution _isotropic_dist{params.center.size()};
  std::uniform_real_distribution<> _uniform_dist{params.inner_radius,
                                                 params.radius};
};

/**
   \class discrete_dist_distribution Stochastic/Random.h "Stochastic/Random.h"
   \brief Discrete distribution over given set of values.
   \tparam Discrete_dist A cdf or pdf object implementing
   <tt>bin_probability(std::size_t)</tt> and <tt>value(Scalar)</tt>.
*/
template <typename Discrete_dist> class discrete_dist_distribution {
public:
  using param_type = std::remove_reference_t<Discrete_dist>;
  using result_type = typename param_type::Value_type;

  discrete_dist_distribution(Discrete_dist &&pdf)
      : _dist_base{std::forward<Discrete_dist>(pdf)},
        _discrete_dist{make_discrete_dist()} {}

  template <typename Generator> result_type operator()(Generator &rng) {
    return _dist_base.value(discrete_dist(rng));
  }

  auto const &dist() const { return _dist_base; }

private:
  Discrete_dist _dist_base;
  std::discrete_distribution<std::size_t> _discrete_dist;

  auto make_discrete_dist() {
    std::vector<double> weights;
    weights.reserve(_dist_base.size());
    for (std::size_t ii = 0; ii < _dist_base.size(); ++ii)
      weights.push_back(_dist_base.bin_probability(ii));

    return std::discrete_distribution<std::size_t>{weights.begin(),
                                                   weights.end()};
  }
};
template <typename Discrete_dist>
discrete_dist_distribution(Discrete_dist &&)
    -> discrete_dist_distribution<Discrete_dist>;

/**
   \class discrete_conditional_dist_distribution Stochastic/Random.h
   "Stochastic/Random.h"
   \brief Discrete conditional distribution over given set of values.
   \tparam Discrete_dist Distribution conditioned on a given value. Can be a
   cdf or pdf implementing <tt>bin_probability(std::size_t, std::size_t)</tt>
   and <tt>value(Scalar, std::size_t)</tt>, where the second argument refers to
   the value being conditioned on.
*/
template <typename Discrete_conditional_dist>
class discrete_conditional_dist_distribution {
public:
  using param_type = std::remove_reference_t<Discrete_conditional_dist>;
  using result_type = typename param_type::Value_type;

  discrete_conditional_dist_distribution(Discrete_conditional_dist &&pdf)
      : _dist_base{std::forward<Discrete_conditional_dist>(pdf)},
        _discrete_dists{make_discrete_dists()} {}

  template <typename Generator>
  result_type operator()(std::size_t jj, Generator &rng) {
    return _dist_base.value(_discrete_dists[jj](rng), jj);
  }

  auto const &dist() const { return _dist_base; }

private:
  Discrete_conditional_dist _dist_base;
  std::vector<std::discrete_distribution<std::size_t>> _discrete_dists;

  auto make_discrete_dists() {
    std::vector<std::discrete_distribution<std::size_t>> discrete_dists;
    discrete_dists.reserve(_dist_base.size());

    for (std::size_t jj = 0; jj < _dist_base.size(); ++jj) {
      std::vector<double> weights;
      weights.reserve(_dist_base.size());
      for (std::size_t ii = 0; ii < _dist_base.size(); ++ii)
        weights.push_back(_dist_base.bin_probability(ii, jj));
      discrete_dists.push_back({weights.begin(), weights.end()});
    }

    return discrete_dists;
  }
};
template <typename Discrete_conditional_dist>
discrete_conditional_dist_distribution(Discrete_conditional_dist &&)
    -> discrete_conditional_dist_distribution<Discrete_conditional_dist>;

/**
   \class uniform_cylinder_distribution Stochastic/Random.h
   "Stochastic/Random.h"
   \brief Uniformly random distribution inside cylinder.
*/
class uniform_cylinder_distribution {
public:
  struct param_type {
    std::vector<double> center_1;
    std::vector<double> center_2;
    double radius;
    double inner_radius{1.};
  } const params;
  using result_type = std::vector<double>;

  uniform_cylinder_distribution(param_type params) : params{params} {}

  uniform_cylinder_distribution(std::vector<double> const &center_1,
                                std::vector<double> const &center_2,
                                double radius, double inner_radius = 1.)
      : uniform_cylinder_distribution{
            param_type{center_1, center_2, radius, inner_radius}} {}

  template <typename Generator> result_type operator()(Generator &rng) {
    std::vector<double> result = geom::polar2cartesian(
        {_uniform_dist_radius(rng), _uniform_dist_angle(rng)});
    result.push_back(_uniform_dist_height(rng));
    result = op::dot(_rotation, result);
    op::plus_inplace(result, params.center_1);
    return result;
  }

private:
  std::vector<double> _height_vector =
      op::minus(params.center_2, params.center_1);
  std::vector<std::vector<double>> _rotation =
      algebra::rotation_matrix_align({0., 0., 1.}, _height_vector);
  std::uniform_real_distribution<> _uniform_dist_height{
      0., op::abs(_height_vector)};
  std::uniform_real_distribution<> _uniform_dist_radius{params.inner_radius,
                                                        params.radius};
  std::uniform_real_distribution<> _uniform_dist_angle{0., 2. * constants::pi};
};

/** \class uniform_rectangle_distribution Stochastic/Random.h
    "Stochastic/Random.h"
    \brief Uniformly random distribution inside parallelipiped in arbitrary
    dimension.
    \note Parallelepiped can be degenerate (sides can be the zero vector).
*/
class uniform_parallelipiped_distribution {
public:
  struct param_type {
    std::vector<double> corner;
    std::vector<std::vector<double>> sides;
  } const params;
  using result_type = std::vector<double>;

  uniform_parallelipiped_distribution(param_type params) : params{params} {}

  uniform_parallelipiped_distribution(
      std::vector<double> const &corner,
      std::vector<std::vector<double>> const &sides)
      : uniform_parallelipiped_distribution{param_type{corner, sides}} {}

  template <typename Generator> result_type operator()(Generator &rng) {
    result_type result;
    result.reserve(params.corner.size());
    for (std::size_t dd = 0; dd < params.corner.size(); ++dd)
      op::plus_inplace(result,
                       op::times_scalar(_uniform_dist(rng), params.sides[dd]));
    op::plus_inplace(result, params.corner);
    return result;
  }

private:
  std::uniform_real_distribution<> _uniform_dist{};
};

/**
   \return PDF of given set of values, normalized to integral equal to number of
   samples.
   \note
   - Output bin values are bin centers.
   - Data smaller than leftmost edge and larger than rightmost edge are included
   in leftmost and rightmost bins, respectively.
*/
template <typename Container>
std::vector<std::pair<double, double>> pdf(std::vector<double> const &bin_edges,
                                           Container const &data) {
  std::vector<std::pair<double, double>> pdf(bin_edges.size() - 1);

  // Populate all bins except last
  auto rt_it = data.cbegin();
  for (std::size_t bin = 0; bin < bin_edges.size() - 2; ++bin) {
    pdf[bin].first = (bin_edges[bin] + bin_edges[bin + 1]) / 2.;
    while (rt_it != data.cend() && *rt_it < bin_edges[bin + 1]) {
      ++pdf[bin].second;
      ++rt_it;
    }
    pdf[bin].second /= bin_edges[bin + 1] - bin_edges[bin];
  }

  // Put everything that is left on last bin
  for (; rt_it != data.cend(); ++rt_it) {
    pdf.back().first =
        (bin_edges[bin_edges.size() - 2] + bin_edges.back()) / 2.;
    ++pdf.back().second;
  }
  pdf.back().second /= bin_edges.back() - bin_edges[bin_edges.size() - 2];

  return pdf;
}

/**
   \brief Pick an event.
   \param probs Cumulative probabilities (not necessarily normalized), <tt>p(0),
   p(0)+p(1), ...</tt>.
   \param rng Random number generator.
   \return Event index \c i with probability <tt>p(i)</tt>.
*/
template <typename Container, typename Engine_t = std::mt19937>
std::size_t pick(Container const &probs, Engine_t &rng) {
  typename Container::value_type rnd =
      probs.back() * std::uniform_real_distribution<double>{}(rng);

  for (std::size_t ev = 0; ev < probs.size(); ++ev)
    if (rnd < probs[ev])
      return ev;

  throw std::runtime_error("Nothing picked!");
}

/**
   \brief Generate a random non-repeating sequence of <tt>subset.size()</tt>
   elements out of <tt>{0, 1, ..., total-1}</tt>.
*/
template <typename Engine_t = std::mt19937>
void randSubset(std::size_t total, std::vector<std::size_t> &subset,
                Engine_t &rng) {
  std::vector<std::size_t> set;
  set.reserve(total);
  for (std::size_t ii = 0; ii < total; ++ii)
    set.push_back(ii);

  std::size_t picked;
  for (auto &ii : subset) {
    picked = std::uniform_int_distribution<std::size_t>{0, set.size() - 1}(rng);
    ii = set[picked];
    useful::swap_erase(set, picked);
  }
}

/**
   \return Random non-repeating sequence of <tt>subset_size</tt> elements out of
   <tt>{0, 1, ..., total-1}</tt>.
*/
template <typename Engine_t = std::mt19937>
std::vector<std::size_t> randSubset(std::size_t total, std::size_t subset_size,
                                    Engine_t &rng) {
  std::vector<std::size_t> subset(subset_size);
  RandSubset(total, subset, rng);

  return subset;
}

/**
   \struct Logspacing Stochastic/Random.h "Stochastic/Random.h"
   \brief Type helper to represent logarithmic spacing.
*/
struct Logspacing {};

/**
   \struct Linspacing Stochastic/Random.h "Stochastic/Random.h"
   \brief Type helper to represent linear spacing.
*/
struct Linspacing {};

/** \return CDF of samples. */
template <typename Spacing, typename Container>
auto cdf(Container &samples, double min_edge, double max_edge,
         std::size_t nr_bins) {
  std::pair<std::vector<double>, std::vector<double>> cdf;
  cdf.first.reserve(nr_bins);
  cdf.second.reserve(nr_bins);
  std::sort(samples);

  std::vector<double> edges;
  if constexpr (std::is_same<Spacing, Logspacing>::value)
    edges = range::logspace(min_edge, max_edge, nr_bins + 1);
  else if constexpr (std::is_same<Spacing, Logspacing>::value)
    edges = range::linspace(min_edge, max_edge, nr_bins + 1);
  else
    throw io::bad_parameters();

  for (std::size_t ii = 1; ii < edges.size(); ++ii)
    cdf.first.push_back(edges[ii] + edges[ii - 1] / 2.);

  auto start =
      std::lower_bound(std::begin(samples), std::end(samples), edges[0]);
  for (std::size_t ii = 1; ii < edges.size(); ++ii) {
    double freq = double(std::lower_bound(std::begin(samples),
                                          std::end(samples), edges[ii]) -
                         start) /
                  samples.size();
    cdf.second.push_back(freq);
  }

  return cdf;
}

/**
   \return CDF of samples.
   \note First and last bin edges determined by a factor of min and max of data.
*/
template <typename Spacing, typename Container>
auto cdf(Container &samples, std::size_t nr_bins, double min_factor = 0.9,
         double max_factor = 1.1) {
  std::pair<std::vector<double>, std::vector<double>> cdf;

  cdf.first.reserve(nr_bins);
  cdf.second.reserve(nr_bins);
  std::sort(std::begin(samples), std::end(samples));

  std::size_t nr_samples = samples.size();
  std::size_t zero_samples = 0;
  if constexpr (std::is_same<Spacing, Logspacing>::value) {
    samples.erase(std::remove(std::begin(samples), std::end(samples), 0.),
                  std::end(samples));
    zero_samples = nr_samples - samples.size();
  }

  double min_edge = min_factor * samples[0];
  double max_edge = max_factor * samples[samples.size() - 1];

  std::vector<double> edges;
  if constexpr (std::is_same<Spacing, Logspacing>::value)
    edges = range::logspace(min_edge, max_edge, nr_bins + 1);
  else if constexpr (std::is_same<Spacing, Linspacing>::value)
    edges = range::linspace(min_edge, max_edge, nr_bins + 1);
  else
    throw io::bad_parameters();

  for (std::size_t ii = 1; ii < edges.size(); ++ii)
    cdf.first.push_back((edges[ii] + edges[ii - 1]) / 2.);

  for (std::size_t ii = 1; ii < edges.size(); ++ii) {
    double freq = double(std::lower_bound(std::begin(samples),
                                          std::end(samples), edges[ii]) -
                         std::begin(samples)) /
                  nr_samples;
    cdf.second.push_back(freq);
  }

  cdf.second[0] += zero_samples / nr_samples;

  return cdf;
}

/** \return Tail CDF of samples, <tt>1 - CDF</tt>. */
template <typename Spacing, typename Container>
auto cdf_tail(Container &samples, double min_edge, double max_edge,
              std::size_t nr_bins) {
  auto cdf_tail = cdf<Spacing>(samples, min_edge, max_edge, nr_bins);
  op::scalar_minus_inplace(1., cdf_tail.second);

  return cdf_tail;
}

/** \return Tail CDF of samples, <tt>1 - CDF</tt>.
    \note First and last bin edges determined by a factor of min and max of
   data.
*/
template <typename Spacing, typename Container>
auto cdf_tail(Container &samples, std::size_t nr_bins, double min_factor = 0.9,
              double max_factor = 1.1) {
  auto cdf_tail = cdf<Spacing>(samples, nr_bins, min_factor, max_factor);
  op::scalar_minus_inplace(1., cdf_tail.second);

  return cdf_tail;
}
} // namespace stochastic

#endif /* STOCHASTIC_RANDOM_H */
