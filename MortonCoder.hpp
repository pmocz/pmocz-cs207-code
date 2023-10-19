#ifndef MORTON_CODER_HPP
#define MORTON_CODER_HPP
/** @file MortonCoder.hpp
 * @brief Define the MortonCoder class for Z-order-curve values, aka Morton
 *   codes.
 */

#include <cassert>

#include "Point.hpp"
#include "BoundingBox.hpp"

/** @class MortonCoder
 * @brief Class representing Z-order-curve values, aka Morton codes.
 *
 * The Z-order curve is a space-filling curve: a one-dimensional curve that
 * fills a multi-dimensional space. Space-filling curves offer advantages for
 * representing points in 3D space. Points near each other in 3D space are
 * likely to have close Morton codes! So we can store points in a map with
 * Z-order value as key, then iterate over nearby Z-order values to fetch
 * points nearby in space.
 *
 * Unfortunately, it's impossible to reduce 3D space to 1D space perfectly:
 * there are some discontinuities in any space-filling curve mapping. But the
 * mapping is still an awesome tool, and some optimizations (see BIGMIN on the
 * Wikipedia page linked below) make them very effective.
 *
 * The MortonCoder class encapsulates a BoundingBox and can be used to translate
 * between spatial Points and Morton codes relative to that BoundingBox.
 *
 * A single Morton code corresponds to a rectangular volume within that
 * BoundingBox, called its <em>cell</em>. Each side of the BoundingBox is
 * divided into 2^L equal-sized cells, for a total of 8^L cells.
 *
 * Read more about the Z-order curve here:
 * http://en.wikipedia.org/wiki/Z-order_curve
 *
 * This class computes maps box numbers to point and visa-versa
 * with respect to a bounding box and the number of equal-volume boxes (8^L).
 * These mappings are performed in O(1) time.
 */
template <int L = 5>
class MortonCoder {
  // Using a 32-bit unsigned int for the code_type
  // means we can only resolve 10 3D levels
  static_assert(L >= 1 && L <= 10, "L (LEVELS) must be between 1 and 10");

 public:
  /** The type to use for the Morton codes -- allows 30-bit codes */
  typedef unsigned code_type;

  /** The number of bits per dimension [octree subdivisions]. #cells = 8^L. */
  static constexpr int levels = L;
  /** The number of cells per side of the bounding box (2^L). */
  static constexpr code_type cells_per_side = code_type(1) << L;
  /** One more than the largest code (8^L). */
  static constexpr code_type end_code = code_type(1) << (3*L);

  /** Construct a MortonCoder with a bounding box. */
  MortonCoder(const BoundingBox& bb)
    : pmin_(bb.min()),
      cell_size_((bb.max() - bb.min()) / cells_per_side) {
    assert(!bb.empty());
  }

  /** Return the MortonCoder's bounding box. */
  BoundingBox bounding_box() const {
    return BoundingBox(pmin_, pmin_ + (cell_size_ * cells_per_side));
  }

  /** Return the bounding box of the cell with Morton code @a c.
   * @pre c < end_code */
  BoundingBox cell(code_type c) const {
    assert(c < end_code);
    Point p = deinterleave(c);
    p *= cell_size_;
    p += pmin_;
    return BoundingBox(p, p + cell_size_);
  }

  /** Return the Morton code of Point @a p.
   * @pre bounding_box().contains(@a p)
   * @post cell(result).contains(@a p) */
  code_type code(const Point& p) const {
    Point s = (p - pmin_) / cell_size_;
    s.x = std::min(std::max(0.0, s.x), double(cells_per_side-1));
    s.y = std::min(std::max(0.0, s.y), double(cells_per_side-1));
    s.z = std::min(std::max(0.0, s.z), double(cells_per_side-1));
    return interleave((unsigned) s.x, (unsigned) s.y, (unsigned) s.z);
  }

  // 0x09249249 = 0b001001001001001001001001001001
  static constexpr code_type coordinate_mask = 0x09249249;
  static constexpr code_type x_mask = coordinate_mask << 0;
  static constexpr code_type y_mask = coordinate_mask << 1;
  static constexpr code_type z_mask = coordinate_mask << 2;

  /** True if min <= idx <= max and idx is inside the box defined
   * by the Morton codes @a min and @a max
   */
  bool is_in_box(code_type idx, code_type min, code_type max) const {
    return (min & x_mask) <= (idx & x_mask) && (idx & x_mask) <= (max & x_mask)
        && (min & y_mask) <= (idx & y_mask) && (idx & y_mask) <= (max & y_mask)
        && (min & z_mask) <= (idx & z_mask) && (idx & z_mask) <= (max & z_mask);
  }

  /** Advance idx to the next box contained in the bounding box defined
   * by the Morton codes @a min and @a max
   *
   * @return result = @a min if @a idx <= @a min
   *                  @a idx if @a idx >= @a max
   *                  smallest i >= @a idx such that is_in_box(i,@a min,@a max)
   *                     and for all j in [@a idx,i), !is_in_box(j,@a min,@a max)
   * @post result > @a max || is_in_box(result,@a min,@a max)
   * @post result >= @a idx
   */
  code_type advance_to_box(code_type idx, code_type min, code_type max) const {
    if (idx >= max) return idx;

    // If outside the box in some coord, record the difference
    code_type delta = 0;
    if      ((idx & x_mask) > (max & x_mask))  delta |= (idx ^ max) & x_mask;
    else if ((idx & x_mask) < (min & x_mask))  delta |= (idx ^ min) & x_mask;
    if      ((idx & y_mask) > (max & y_mask))  delta |= (idx ^ max) & y_mask;
    else if ((idx & y_mask) < (min & y_mask))  delta |= (idx ^ min) & y_mask;
    if      ((idx & z_mask) > (max & z_mask))  delta |= (idx ^ max) & z_mask;
    else if ((idx & z_mask) < (min & z_mask))  delta |= (idx ^ min) & z_mask;

    // Delta is only zero if idx is in the box
    if (delta == 0) return idx;

    // Smear into a low bit mask, i.e. 0000111111111111
    delta = smear_low_1(delta >> 1);
    //if ((delta+1) & idx) {   // The idx bit above delta mask is one, need zero
      // Chi masks high bits we cannot carry into
      code_type chi = ~smear_low_3(idx ^ max);
      // The first 0 in idx and chi that is higher than delta
      delta = ~(idx | chi | delta);
      delta = (delta & -delta) - 1;
      //}

    // Flip zero bit of idx and zero all lower bits
    idx = (idx | delta) + 1;

    // For each coordinate, if idx is low set to min
    if ((idx & x_mask) < (min & x_mask))  idx |= min & x_mask;
    if ((idx & y_mask) < (min & y_mask))  idx |= min & y_mask;
    if ((idx & z_mask) < (min & z_mask))  idx |= min & z_mask;

    return idx;
  }

 private:

  /** The minimum of the MortonCoder bounding box. */
  Point pmin_;
  /** The extent of a single cell. */
  Point cell_size_;

  /** Spreads the bits of a 10-bit number so that there are two 0s
   *  in between each bit.
   * @param x 10-bit integer
   * @return 28-bit integer of form 0b0000X00X00X00X00X00X00X00X00X00X,
   * where the X's are the original bits of @a x
   */
  inline code_type spread_bits(code_type x) const {
    //......................9876543210
    x = (x | (x << 10)) & 0x000f801f; //............98765..........43210
    x = (x | (x <<  4)) & 0x00e181c3; //........987....56......432....10
    x = (x | (x <<  2)) & 0x03248649; //......98..7..5..6....43..2..1..0
    x = (x | (x <<  2)) & 0x09249249; //....9..8..7..5..6..4..3..2..1..0
    return x;
  }

  /** Interleave the bits of n into x, y, and z.
   * @pre x = [... x_2 x_1 x_0]
   * @pre y = [... y_2 y_1 y_0]
   * @pre z = [... z_2 z_1 z_0]
   * @post n = [... z_1 y_1 x_1 z_0 y_0 x_0]
   */
  inline code_type interleave(unsigned x, unsigned y, unsigned z) const {
    return spread_bits(x) | (spread_bits(y) << 1) | (spread_bits(z) << 2);
  }

  /** Does the inverse of spread_bits, extracting a 10-bit number from
   * a 28-bit number.
   * @param x 28-bit integer of form 0bYYYYXYYXYYXYYXYYXYYXYYXYYXYYXYYX
   * @return 10-bit integer of form 0b00...000XXXXXXXXXX,
   * where the X's are every third bit of @a x
   */
  inline code_type compact_bits(code_type x) const {
    x = ( x ) & 0x09249249;  //....9..8..7..5..6..4..3..2..1..0
    x = (x | (x >>  2)) & 0x03248649;  //......98..7..5..6....43..2..1..0                                          
    x = (x | (x >>  2)) & 0x00e181c3;  //........987....56......432....10                                       
    x = (x | (x >>  4)) & 0x000f801f;  //............98765..........43210                                          
    x = (x | (x >> 10)) & 0x000003FF;  //......................9876543210   
    return x;
  }

  /** Deinterleave the bits from n into a Point.
   * @pre n = [... n_2 n_1 n_0]
   * @post result.x = [... n_6 n_3 n_0]
   * @post result.y = [... n_7 n_4 n_1]
   * @post result.z = [... n_8 n_5 n_2]
   */
  inline Point deinterleave(code_type c) const {
    return Point(compact_bits(c), compact_bits(c >> 1), compact_bits(c >> 2));
  }


  /** Smears the bits in c into the low bits by steps of one
   *
   * Example: 00011100100 -> 000111111111
   */
  inline code_type smear_low_1(code_type c) const {
    c |= c >>  1;
    c |= c >>  2;
    c |= c >>  4;
    c |= c >>  8;
    c |= c >> 16;
    return c;
  }

  /** Smears the bits in c into the low bits by steps of three
   *
   * Example: 0000010000000000 -> 0000010010010010
   */
  inline code_type smear_low_3(code_type c) const {
    c |= c >>  3;
    c |= c >>  6;
    c |= c >> 12;
    c |= c >> 24;
    return c;
  }
};

#endif
