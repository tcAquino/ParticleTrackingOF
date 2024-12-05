/**
   \file Algebra/Algebra.h
   \author Tomás Aquino
   \date 08/06/2019
   \brief Linear algebra algorithms
*/

#ifndef ALGEBRA_ALGEBRA_H
#define ALGEBRA_ALGEBRA_H

#include "General/Operations.h"
#include <cmath>
#include <vector>

/** \namespace useful Linear algebra algorithms. */
namespace algebra {
/**
   \brief Orthogonal basis from Gram-schmidt decomposition, given complete set
   of linearly-independent vectors. First vector is kept after normalization.
*/
inline std::vector<std::vector<double>>
gram_schmidt(std::vector<std::vector<double>> const &input) {
  std::size_t nr_vectors = input.size();

  std::vector<std::vector<double>> output;
  output.reserve(nr_vectors);
  output.push_back(op::div_scalar(input[0], op::abs(input[0])));
  for (std::size_t ii = 1; ii < nr_vectors; ++ii) {
    output.push_back(input[ii]);
    for (std::size_t jj = 0; jj < ii; ++jj)
      op::minus_inplace(output[ii],
                        op::times_scalar(op::dot(output[jj], output[ii]) /
                                         op::abs_sq(output[jj]),
                                         output[jj]));
    op::div_scalar_inplace(output[ii], op::abs(output[ii]));
  }

  return output;
}

/**
   \brief Orthogonal basis from Gram-Schmidt decomposition, preserving given
   direction.
   \note Currently does not work for direction aligned with a Cartesian plane
   but not with a Cartesian direction.
*/
inline std::vector<std::vector<double>>
gram_schmidt(std::vector<double> const &input) {
  std::size_t dim = input.size();

  std::vector<std::vector<double>> input_vecs;
  input_vecs.reserve(dim);
  input_vecs.push_back(op::normalize(input));

  std::vector<double> basis_vector(dim);
  for (std::size_t dd = 0; dd < dim && input_vecs.size() < dim; ++dd) {
    basis_vector[dd] = 1.;
    if (input_vecs.front() != basis_vector)
      input_vecs.push_back(basis_vector);
    basis_vector[dd] = 0.;
  }

  return gram_schmidt(input_vecs);
}

/**
   \return Rotation matrix (non-unique) that aligns \p to_align with \p
   align_with.
*/
template <typename Container1 = std::vector<double>,
          typename Container2 = Container1>
std::vector<std::vector<double>>
rotation_matrix_align(Container1 const &to_align,
                      Container2 const &align_with) {
  auto axis_direction =
      op::plus(op::normalize(to_align), op::normalize(align_with));
  auto axis_norm = op::abs(axis_direction);
  auto rotation = op::outer_product(axis_direction, axis_direction);
  for (auto &row : rotation)
    op::times_scalar_inplace(2. / axis_norm, row);
  for (std::size_t ii = 0; ii < rotation.size(); ++ii)
    rotation[ii][ii] -= 1.;

  return rotation;
}

/** \brief Thomas Algorithm for tridiagonal systems. */
inline void solve_tridiag(std::vector<double> const &l_diag,
                          std::vector<double> const &diag,
                          std::vector<double> const &u_diag,
                          std::vector<double> const &ind,
                          std::vector<double> &sol) {
  std::size_t dim = diag.size();
  sol.resize(dim);

  // Modified coefficients
  std::vector<double> u_diag_new(dim);
  std::vector<double> ind_new(dim);
  u_diag_new[0] = u_diag[0] / diag[0];
  ind_new[0] = ind[0] / diag[0];
  double temp;
  for (size_t ii = 1; ii < dim; ++ii) {
    temp = diag[ii] - l_diag[ii] * u_diag_new[ii - 1];
    u_diag_new[ii] = u_diag[ii] / temp;
    ind_new[ii] = (ind[ii] - l_diag[ii] * ind_new[ii - 1]) / temp;
  }

  // Back-substitution
  sol.back() = ind_new.back();

  for (size_t ii = dim - 2; ii != 0; ii--)
    sol[ii] = ind_new[ii] - u_diag_new[ii] * sol[ii + 1];
}

/**
   \brief Uses Sherman-Morrison and the Thomas Algorithm (see solve_tridiag())
   to solve a cyclic tridiagonal system <tt>Ax = ind</tt>, with
   <tt>alpha = A_{ n-1, 0 }, beta = A_{ 0, n-1 }</tt>.
*/
inline void solve_cyclic_tridiag(double alpha, double beta,
                                 std::vector<double> const &l_diag,
                                 std::vector<double> const &diag,
                                 std::vector<double> const &u_diag,
                                 std::vector<double> const &ind,
                                 std::vector<double> &sol) {
  std::size_t nr_points = diag.size();
  sol.resize(nr_points);
  double gamma = -diag[0];

  // Sherman-Morrison v = (1 0 0 ... 0 beta / gamma)
  std::vector<double> zz, uu(nr_points, 0.), diag_mod;
  diag_mod = diag;
  diag_mod[0] -= gamma;
  diag_mod.back() -= alpha * beta / gamma;
  uu[0] = gamma;
  uu.back() = alpha;

  solve_tridiag(l_diag, diag_mod, u_diag, ind, sol);
  solve_tridiag(l_diag, diag_mod, u_diag, uu, zz);
  double factor = (sol[0] + beta * sol.back() / gamma) /
                  (1 + zz[0] + beta * zz.back() / gamma);

  for (size_t ii = 0; ii < nr_points; ++ii)
    sol[ii] -= factor * zz[ii];
}

/**
   \brief Find the determinant of a tridiagonal matrix using the continuant
   recurrence relation.
*/
inline double det_tridiag(std::vector<double> const &l_diag,
                          std::vector<double> const &diag,
                          std::vector<double> const &u_diag) {
  double cont_0 = 1.;
  double cont_1 = diag[0];
  double cont_2 = 0.;

  for (size_t ii = 2; ii <= diag.size(); ++ii) {
    cont_2 = diag[ii - 1] * cont_1 - u_diag[ii - 2] * l_diag[ii - 1] * cont_0;
    cont_0 = cont_1;
    cont_1 = cont_2;
  }

  return cont_2;
}
} // namespace algebra

#endif /* ALGEBRA_ALGEBRA_H */
