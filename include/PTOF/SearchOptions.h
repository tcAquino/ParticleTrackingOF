/**
 * @file   SearchOptions.h
 * @author Tomás Aquino <tomas.aquino@csic.es>
 * @date   Thu Feb 17 00:00:00 2022
 *
 * @brief Helpers for choosing neighbor search options.
 */

#ifndef PTOF_SEARCHOPTIONS_H
#define PTOF_SEARCHOPTIONS_H

namespace ptof {
/** @brief Options for mesh searching. */
struct SearchOptions {
  /** @brief Level of neighbors to check before searching mesh for position. */
  template <int level> struct NeighborPrecheck {
    static constexpr int neighbor_check_level{level};
  };

  /**
   * @brief Dot not check neighbors, even self, before searching mesh for
   *        position.
   */
  using NoNeighborPrecheck = NeighborPrecheck<-1>;

  /**
   * @brief Check if position is in given cell before searching mesh for
   *        position.
   */
  using ZeroNeighborPrecheck = NeighborPrecheck<0>;

  /**
   * @brief Check if position is in given cell and in first orthogonal neighbors
   *        before searching mesh for position.
   */
  using FirstNeighborPrecheck = NeighborPrecheck<1>;

  /**
   * @brief Check if position is in given cell and in first and second
   *        orthogonal neighbors before searching mesh for position.
   */
  using SecondNeighborPrecheck = NeighborPrecheck<2>;
};
} // namespace ptof

#endif /* PTOF_SEARCHOPTIONS_H */
