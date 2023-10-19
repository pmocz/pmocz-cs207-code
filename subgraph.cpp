/**
 * @file subgraph.cpp
 * Test script for making a subgraph from our Graph
 * my predicate is designed for the large data set
 * please run: ./subgraph data/large.nodes  data/large.tets 
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>
#include <iterator>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"

#include "Graph.hpp"

#include <math.h>       /* fmod */

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

/** Test predicate for HW1 #4 */
struct SlicePredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return n.position().x < 0;
  }
};

/** My own predicate for HW1 #4 
  * uses fmod to take cool slices out of the data
  * designed for the 'large' data set
  * please run: ./subgraph data/large.nodes  data/large.tets 
  */
struct ModPredicate {
  template <typename NODE>
  bool operator()(const NODE& n) {
    return fmod(n.position().x,0.1) < 0.05 && 
           fmod(n.position().y,0.1) < 0.05 && 
           fmod(n.position().z,0.1) < 0.05;
  }
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

  // Use the filter_iterator to plot an induced subgraph.
  // SlicePredicate pred; // slice predicate
  ModPredicate pred; // my own predicate
  auto fi_begin = make_filtered(graph.node_begin(), graph.node_end(), pred);
  auto fi_end = make_filtered(graph.node_end(), graph.node_end(), pred);
  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(fi_begin,fi_end,node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map );
  viewer.center_view();

  return 0;
}
