/**
 * @file hw1q5_cool_image.cpp
 * Script for generating a cool image
 * please run on the large data set:
 * please run: ./hw1q5_cool_image data/large.nodes  data/large.tets 
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"

#include <vector>
#include <deque>
#include <fstream>

#include <math.h>       /* fabs */

/** An iterator that skips over elements of another iterator based on whether
 * those elements satisfy a predicate.
 *
 * Given an iterator range [@a first, @a last) and a predicate @a pred,
 * this iterator models a filtered range such that all i with
 * @a first <= i < @a last and @a pred(*i) appear in order of the original range.
 */
template <typename Pred, typename It>
class filter_iterator
    : private equality_comparable<filter_iterator<Pred,It>> {
 public:
  // Get all of the iterator traits and make them our own
  typedef typename std::iterator_traits<It>::value_type        value_type;
  typedef typename std::iterator_traits<It>::pointer           pointer;
  typedef typename std::iterator_traits<It>::reference         reference;
  typedef typename std::iterator_traits<It>::difference_type   difference_type;
  typedef typename std::input_iterator_tag                     iterator_category;

  typedef filter_iterator<Pred,It> self_type;

  // Constructor
  filter_iterator(const Pred& p, const It& first, const It& last)
      : p_(p), it_(first), end_(last) {
        while (it_!=end_ && !p_(*it_))
          ++it_;
  }

  /** dereference a filter iterator to return dereferenced (original) iterator 
   * @pre original iterator iterator must refer to a valid object
   * @post the returned object satisfies the predicate
   */
  value_type operator*() const {return *it_;};
  /** increment the filter iterator to refer to the next iterator that satisfies the predicate*/
  self_type& operator++() {
    ++it_;
    // 'fix' step. increment further if predicate is not met
    while (it_!=end_ && !p_(*it_))
      ++it_;
    return *this;
  }
  /** check whether two filer iterators refer to the same iterator*/
  bool operator==(const self_type& fi) const {return it_ == fi.it_;};
  
 private:
  Pred p_;
  It it_;
  It end_;
};

/** Helper function for constructing filter_iterators.
 *
 * Usage:
 * // Construct an iterator that filters odd values out and keeps even values.
 * std::vector<int> a = ...;
 * auto it = make_filtered(a.begin(), a.end(), [](int k) {return k % 2 == 0;});
 */
template <typename Pred, typename Iter>
filter_iterator<Pred,Iter> make_filtered(const Iter& it, const Iter& end,
                                         const Pred& p) {
  return filter_iterator<Pred,Iter>(p, it, end);
}


/** My own predicate for HW1 #5
  * takes a slice of the data
  */
struct MySlice1Predicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return n.position().x > 0. && n.position().z < 0.1;
  }
};

struct MySlice2Predicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return n.position().x > 0. && n.position().z > 0.2;
  }
};

struct MySlice3Predicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return n.position().x < -0.1;
  }
};



/** Comparator that compares the distance from a given point p.
 */
struct MyComparator {
   Point p_;
   MyComparator(const Point& p) : p_(p) {
   };

   template <typename NODE>
   bool operator()(const NODE& node1, const NODE& node2) const {
    // Return true if node1 is closer to p_ than node2
    if ( norm(node1.position()-p_) < norm(node2.position()-p_) )
      return true;
    return false;
  }
};


/** Calculate shortest path lengths in @a g from the nearest node to @a point.
 * @param[in,out] g Input graph
 * @param[in] point Point to find the nearest node to.
 * @post Graph has modified node values indicating the minimum path length
 *           to the nearest node to @a point
 * @post Graph nodes that are unreachable to the nearest node to @a point have
 *           the value() -1.
 * @return The maximum path length found.
 *
 * Finds the nearest node to @a point and treats that as the root node for a
 * breadth first search.
 * This sets node's value() to the length of the shortest path to
 * the root node. The root's value() is 0. Nodes unreachable from
 * the root have value() -1.
 */
int shortest_path_lengths(Graph<int>& g, const Point& point) {
  // initialize all node values to -1
  for(Graph<int>::node_iterator ni=g.node_begin(); ni!=g.node_end(); ++ni)
    (*ni).value() = -1;
  // find the root node, set its value() to 0
  Graph<int>::node_type root = *std::min_element(g.node_begin(),g.node_end(),MyComparator(point));
  root.value() = 0;
  // breadth-first search
  int max_path_length = 0;
  std::deque<Graph<int>::node_type> Q;
  Graph<int>::node_type trial_node;
  Graph<int>::incident_iterator ii;
  Q.push_back(root);
  while (Q.size()>0) {
    trial_node = Q.front();  
    if (trial_node.value() < 1) {
      for (ii = trial_node.edge_begin(); ii!= trial_node.edge_end(); ++ii ) {     
        if ((*ii).node2().value() == -1) {
          Q.push_back((*ii).node2());
        } 
        else {
          if (trial_node.value() == -1 || (*ii).node2().value() + 1 < trial_node.value())
            trial_node.value() = (*ii).node2().value() + 1;
        }
      }
      if (trial_node.value() > max_path_length)
        max_path_length = trial_node.value();
    }
    Q.pop_front();    
  }

  return max_path_length;
}




/** A color functor that returns heatmap colors for values it receives
 * long paths lengths in appear in blue-purple and short path lengths in red
 */
struct HeatMap {
  template <typename NODE>
  CS207::Color operator()(const NODE& n) {
    if (n.value() == -1)
      return 0.;
    return color_.make_heat(1.-((float) n.value()/norm_));
  }
    HeatMap(const float norm)  : norm_(norm) {};
  private:
    float norm_;
    CS207::Color color_;
};

struct Purple {
  template <typename NODE>
  CS207::Color operator()(const NODE& n) {
    return color_.make_rgb(0.9,0.,1.0);
  }
    CS207::Color color_;
};

struct RadialColor {
  template <typename NODE>
  CS207::Color operator()(const NODE& n) {
    return color_.make_rgb(0.7+0.3*fabs(n.position().x/norm(n.position())),0.7+0.3*fabs(n.position().y/norm(n.position())),0.7+0.3*fabs(n.position().z/norm(n.position())));
  }
    CS207::Color color_;
};


int main(int argc, char* argv[])
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  // Construct a Graph
  typedef Graph<int> GraphType;
  GraphType graph;
  std::vector<GraphType::node_type> nodes;

  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,4> t;
  while (CS207::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);

  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // select subgraph using predicate
  // Use shortest_path_lengths to set the node values to the path lengths
  // Construct a Color functor and view with the SDLViewer    
  auto fi1_begin = make_filtered(graph.node_begin(), graph.node_end(), MySlice1Predicate());
  auto fi1_end = make_filtered(graph.node_end(), graph.node_end(), MySlice1Predicate());
  auto node_map = viewer.empty_node_map(graph);
  Point root = Point();
  root.x = 1.0;
  root.y = -0.1;
  root.z = 1.0;
  int max_path_length = shortest_path_lengths(graph, root);
  auto fi2_begin = make_filtered(graph.node_begin(), graph.node_end(), MySlice2Predicate());
  auto fi2_end = make_filtered(graph.node_end(), graph.node_end(), MySlice2Predicate());
  auto fi3_begin = make_filtered(graph.node_begin(), graph.node_end(), MySlice3Predicate());
  auto fi3_end = make_filtered(graph.node_end(), graph.node_end(), MySlice3Predicate()); 
  
  viewer.add_nodes(fi1_begin,fi1_end, RadialColor(), node_map);
  viewer.add_nodes(fi2_begin,fi2_end, HeatMap((float) max_path_length), node_map);
  viewer.add_nodes(fi3_begin,fi3_end, Purple(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map );
  
  viewer.center_view();  
  
  return 0;
}
