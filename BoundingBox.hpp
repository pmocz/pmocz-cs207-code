#ifndef CS207_BOUNDINGBOX_HPP
#define CS207_BOUNDINGBOX_HPP 1
/** @file BoundingBox.hpp
 * @brief Define the BoundingBox class for 3D bounding boxes. */

#include "Point.hpp"

/** @class BoundingBox
 * @brief Class representing 3D bounding boxes.
 *
 * A BoundingBox is a 3D volume. Its fundamental operations are contains(),
 * which tests whether a point is in the volume, and operator+=(), which
 * extends the volume as necessary to ensure that the volume contains a point.
 *
 * BoundingBoxes are implemented as boxes -- 3D rectangular cuboids -- whose
 * sides are aligned with the X, Y, and Z axes.
 */
class BoundingBox {
 public:

  static constexpr double min_dimension = 0.000001;

  /** Construct an empty bounding box. */
  BoundingBox()
    : empty_(true) {
  }
  /** Construct the minimal bounding box containing @a p.
   * @post contains(@a p) && min() == @a p && max() == @a p */
  explicit BoundingBox(const Point& p)
    : empty_(false), min_(p), max_(p) {
  }
  /** Construct the minimal bounding box containing a given sphere.
   * @param[in] center center of the sphere
   * @param[in] radius radius of the sphere */
  BoundingBox(const Point& center, double radius)
    : empty_(false) {
    Point vec = Point(1, 1, 1) * fabs(radius);
    min_ = center - vec;
    max_ = center + vec;
  }
  /** Construct the minimal bounding box containing @a p1 and @a p2.
   * @post contains(@a p1) && contains(@a p2) */
  BoundingBox(const Point& p1, const Point& p2)
    : empty_(false), min_(p1), max_(p1) {
    *this |= p2;
  }
  /** Construct a bounding box containing the points in [first, last). */
  template <typename IT>
  BoundingBox(IT first, IT last)
    : empty_(true) {
    insert(first, last);
  }

  /** Test if the bounding box is empty (contains no points). */
  bool empty() const {
    return empty_;
  }
  /** Test if the bounding box is nonempty.
   *
   * This function lets you write code such as "if (b) { ... }" or
   * "if (box1 & box2) std::cout << "box1 and box2 intersect\n". */
  operator const void*() const {
    return empty_ ? 0 : this;
  }

  /** Return the minimum corner of the bounding box.
   * @post empty() || contains(min())
   *
   * The minimum corner has minimum x, y, and z coordinates of any corner.
   * An empty box has min() == Point(). */
  const Point& min() const {
    return min_;
  }
  /** Return the maximum corner of the bounding box.
   * @post empty() || contains(max()) */
  const Point& max() const {
    return max_;
  }
  /** Return the dimensions of the bounding box.
   * @return max() - min(), adjusted as follows: If !empty(), then all
   *   coordinates of dimensions() are at least min_dimension. */
  Point dimensions() const {
    Point d = max_ - min_;
    return Point(std::max(d.x, min_dimension),
		 std::max(d.y, min_dimension),
		 std::max(d.z, min_dimension));
  }
  /** Return the center of the bounding box. */
  Point center() const {
    return (min_ + max_) / 2;
  }

  /** Test if point @a p is in the bounding box. */
  bool contains(const Point& p) const {
    return !empty()
      && p.x >= min_.x && p.y >= min_.y && p.z >= min_.z
      && p.x <= max_.x && p.y <= max_.y && p.z <= max_.z;
  }
  /** Test if @a box is entirely within this bounding box.
   *
   * Returns false if @a box.empty(). */
  bool contains(const BoundingBox& box) const {
    return !box.empty() && contains(box.min()) && contains(box.max());
  }
  /** Test if @a box intersects this bounding box. */
  bool intersects(const BoundingBox& box) const {
    return !empty() && !box.empty()
      && box.min_.x <= max_.x && box.max_.x >= min_.x
      && box.min_.y <= max_.y && box.max_.y >= min_.y
      && box.min_.z <= max_.z && box.max_.z >= min_.z;
  }

  /** Extend the bounding box to contain @a p.
   * @post contains(@a p) is true
   * @post if old contains(@a x) was true, then new contains(@a x) is true */
  BoundingBox& operator|=(const Point& p) {
    if (empty_) {
      empty_ = false;
      min_ = max_ = p;
    } else {
      min_.x = std::min(min_.x, p.x);
      min_.y = std::min(min_.y, p.y);
      min_.z = std::min(min_.z, p.z);
      max_.x = std::max(max_.x, p.x);
      max_.y = std::max(max_.y, p.y);
      max_.z = std::max(max_.z, p.z);
    }
    return *this;
  }
  /** Extend the bounding box to contain @a box.
   * @post contains(@a box) is true
   * @post if old contains(@a x) was true, then new contains(@a x) is true */
  BoundingBox& operator|=(const BoundingBox& box) {
    if (!box.empty())
      (*this |= box.min()) |= box.max();
    return *this;
  }
  /** Extend the bounding box to contain the points in [first, last). */
  template <typename IT>
  BoundingBox& insert(IT first, IT last) {
    while (first != last) {
      *this |= *first;
      ++first;
    }
    return *this;
  }

  /** Intersect this bounding box with @a box. */
  BoundingBox& operator&=(const BoundingBox& box) {
    if (empty() || box.empty()) {
      clear();
      return *this;
    }
    for (int i = 0; i < 3; ++i) {
      if (min_[i] > box.max_[i] || max_[i] < box.min_[i]) {
        clear();
        return *this;
      }
      if (min_[i] < box.min_[i])
        min_[i] = box.min_[i];
      if (max_[i] > box.max_[i])
        max_[i] = box.max_[i];
    }
    return *this;
  }

  /** Clear the bounding box.
   * @post empty() */
  void clear() {
    empty_ = true;
    min_ = max_ = Point();
  }

 private:
  bool empty_;
  Point min_;
  Point max_;
};

/** Write a BoundingBox to an output stream.
 *
 * An empty BoundingBox is written as "[]". A nonempty BoundingBox is
 * written as "[minx miny minz:maxx maxy maxz]", where all coordinates
 * are double-precision numbers.
 */
inline std::ostream& operator<<(std::ostream& s, const BoundingBox& box) {
  if (box.empty())
    return (s << '[' << ']');
  else
    return (s << '[' << box.min() << ':' << box.max() << ']');
}

/** Return a bounding box that contains @a box and @a p. */
BoundingBox operator|(BoundingBox box, const Point& p) {
  return box |= p;
}
/** Return the union of @a box1 and @a box2. */
BoundingBox operator|(BoundingBox box1, const BoundingBox& box2) {
  return box1 |= box2;
}
/** Return a bounding box that contains @a p1 and @a p2. */
BoundingBox operator|(const Point& p1, const Point& p2) {
  return BoundingBox(p1, p2);
}

/** Return the intersection of @a box1 and @a box2. */
BoundingBox operator&(BoundingBox box1, const BoundingBox& box2) {
  return box1 &= box2;
}

#endif
