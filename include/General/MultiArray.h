/**
 \file General/MultiArray.h
 \author Tomás Aquino
 \date 28/09/2017
*/

#ifndef MultiArray_h
#define MultiArray_h

#include "General/Modular.h"
#include <cassert>
#include <cstddef>
#include <functional>
#include <numeric>
#include <vector>

namespace useful {
/** \class MultiArray General/MultiArray.h "General/MultiArray.h"
 \brief Array container with arbitrary rank.
 \details Different dimensions may have different sizes. The elements are
 indexed, traversing the container in the order (0,...,0), (0,...,1), ...,
 (0,..,1,0) in either 1D or using a multi-dimensional index. The elements are
 continguous in memory in the order of the 1D index.
*/
template <typename value_t> class MultiArray {
public:
  using value_type = value_t; /**< Type of held values. */

  /** Constructor.
   \param dimensions Number of elements along each dimension.
  */
  explicit MultiArray(std::vector<std::size_t> dimensions)
      : _dimensions(dimensions),
        _values(std::accumulate(dimensions.begin(), dimensions.end(), 1,
                                std::multiplies<std::size_t>())) {}

  /** Constructor.
   \param dimensions Number of elements along each dimension.
   \param value Value to assign to all entries.
  */
  explicit MultiArray(std::vector<std::size_t> dimensions, value_type value)
      : _dimensions(dimensions),
        _values(std::accumulate(dimensions.begin(), dimensions.end(), 1,
                                std::multiplies<std::size_t>()),
                value) {}

  /** \return Value at multi-dimensional index \p indexes. */
  value_type const &operator[](std::vector<std::size_t> const &indexes) const {
    return _values[computeIndex(indexes)];
  }

  /** \return Value at multi-dimensional index \p indexes. */
  value_type &operator[](std::vector<std::size_t> const &indexes) {
    return _values[computeIndex(indexes)];
  }

  /** \return Value at 1D index \p index, traversing the container in the order
   * (0,...,0), (0,...,1), ..., (0,..,1,0), etc.  */
  value_type const &operator[](std::size_t index) const {
    return _values[index];
  }

  /** \return Value at 1D index \p index, which traverses the container in the
   * order (0,...,0), (0,...,1), ..., (0,..,1,0), ....  */
  value_type &operator[](std::size_t index) { return _values[index]; }

  /** \return Convert multi-dimensional index \p indexes to 1D index, which
   * traverses the container in the order (0,...,0), (0,...,1), ..., (0,..,1,0),
   * ....  */
  std::size_t computeIndex(std::vector<std::size_t> const &indexes) const {
    return modular::positionToIndex(indexes, _dimensions);
  }

  /** \return Convert 1D index \p index, which traverses the container in the
   * order (0,...,0), (0,...,1), ..., (0,..,1,0), ..., to multi-dimensional
   * index.  */
  std::vector<std::size_t> computeIndexes(std::size_t index) const {
    std::vector<std::size_t> indexes;
    computeIndexes(index, indexes);

    return indexes;
  }

  /** \brief Convert 1D index \p index to multi-dimensional index.
   \param index Index to convert, which traverses the container in the order
   (0,...,0), (0,...,1), ..., (0,..,1,0). \param indexes Multidimensional index
   is returned here. \note: \p indexes must be passed with the right size. */
  void computeIndexes(std::size_t index,
                      std::vector<std::size_t> &indexes) const {
    modular::indexToPosition(index, _dimensions, indexes);
  }

  /** \return Total number of elements. */
  std::size_t size() const { return _values.size(); }

  /** \return Number of elements along dimension \p dd. */
  std::size_t size(std::size_t dd) const { return _dimensions[dd]; }

  /** \return Rank (number of dimensions). */
  std::size_t rank() const { return _dimensions.size(); }

  /** \return Number of elements along each dimension. */
  std::vector<std::size_t> const &sizes() const { return _dimensions; }

  /** \return Iterator to beginning. Container is traverse in the order
   * (0,...,0), (0,...,1), ..., (0,..,1,0). */
  auto begin() { return _values.begin(); }

  /** \return Iterator to end. Container is traverse in the order (0,...,0),
   * (0,...,1), ..., (0,..,1,0). */
  auto end() { return _values.end(); }

  /** \return Iterator to beginning. Container is traverse in the order
   * (0,...,0), (0,...,1), ..., (0,..,1,0). */
  auto begin() const { return _values.begin(); }

  /** \return Iterator to end. Container is traverse in the order (0,...,0),
   * (0,...,1), ..., (0,..,1,0). */
  auto end() const { return _values.end(); }

  /** \return Constant iterator to beginning. Container is traverse in the order
   * (0,...,0), (0,...,1), ..., (0,..,1,0). */
  auto cbegin() const { return _values.cbegin(); }

  /** \return Constant iterator to end. Container is traverse in the order
   * (0,...,0), (0,...,1), ..., (0,..,1,0). */
  auto cend() const { return _values.cend(); }

private:
  std::vector<std::size_t> _dimensions; /**< Size along each dimension. */
  std::vector<value_type> _values;      /**< Values in container. */
};
} // namespace useful

#endif /* MultiArray_h */
