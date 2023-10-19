#ifndef CS207_GRAPH_HPP
#define CS207_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include "CS207/Util.hpp"
#include "Point.hpp"
#include "MortonCoder.hpp"

#include <algorithm>
#include <vector>
#include <cassert>

using namespace std;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {

 private:

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  
  struct internal_node;
  struct internal_edge;
  struct code_uid_pair;
  struct sorted_morton_cont_comp;
  
 public:

  // PUBLIC TYPE DEFINITIONS

  /** Type of the node value. */
  typedef V node_value_type;
 
  /** Type of the edge value. */
  typedef E edge_value_type;
 
  /** Type of this graph. */
  typedef Graph graph_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Type of indexes and sizes. Return type of Node::index() and
      Graph::num_nodes(), argument type of Graph::node. */
  typedef unsigned size_type;

  /**Morton code type **/
  typedef MortonCoder<5>::code_type code_type;  
  
  /** Type of node iterators, which iterate over all graph nodes. */
  class node_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class edge_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class incident_iterator;

  // CONSTRUCTOR AND DESTRUCTOR

  /** Construct an empty graph. 
   * nodes_ is a vector of nodes of size num_nodes()
   * uid2idx_ is a vector to convert from node unique id's (which are always fixed) to indices (which go from 1,...,num_nodes())
   * idx2uid_ is a vector to convert from node indices to unique id's
   * next_node_uid_ is the next unique id to assign to a new node added to the graph. Initialized to 0
   * edges_ is an adjacency list container. It has size num_nodes(), and each entry of edges_ is a vector containing the unique id's of the node's neighbors
   * num_edges_ keeps track of the number of edges in the graph. Initialized to 0
   * domain_bbox_ bounding box of domain
   * morton_coder_ is Morton coder object that generates morton codes
   * morton_is_valid_ is a boolean for whether Morton codes are sorted
   * sorted_morton_cont_ is sorted container for Morton codes of Nodes
   */
  Graph() 
    : nodes_(), uid2idx_(), idx2uid_(), next_node_uid_(0), edges_(),  num_edges_(0), domain_bbox_(BoundingBox(Point(-5,-5,-5),Point(5,5,5))), morton_coder_(MortonCoder<5>(domain_bbox_)), morton_is_valid_(false), sorted_morton_cont_() {
  }
  /** Default destructor */
  ~Graph() = default;

  // NODES

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node : private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Node x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */
    Node() {
    }

    /** Return this node's position. */
    Point position() const {
      return graph_->nodes_[index()].position;
    }
    
    /** Set this node's position. 
     * @pre node is valid
     */
    void set_position(const Point& p) {
      code_type prev_code = graph_->morton_coder_.code(position());
      // set new node position
      graph_->nodes_[index()].position = p;
      // check whether Morton code container is invalidated
      code_type new_code = graph_->morton_coder_.code(p);
      if (new_code != prev_code)
        graph_->morton_is_valid_ = false;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_->uid2idx_[uid_];
    }

    /** return a reference to the node value of a valid node*/
    node_value_type& value() {
      return graph_->nodes_[index()].node_value;
    }
    /** return a reference to the node value of a valid node*/
    const node_value_type& value() const {
      return graph_->nodes_[index()].node_value;
    }
    /** comparator (==) for the Node
     * @pre Node @a n is a valid node
     * @post returns true iff this node is part of same Graph as n and has same index()
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (index() == n.index());
    }
    /** comparator (<) for the Node
     * @pre Node n is a valid node
     * @post returns true iff this node is part of same Graph as n and has smaller index() than Node n
     * or if address of this node's graph is smaller than that of Node n
     */
    bool operator< (const Node& n) const {
      return ((graph_ == n.graph_) && (index() < n.index())) || (graph_ < n.graph_);
    }

    
    /** return the degree of this node (the number of nodes this node is connected to by an edge)
     * Complexity: O(1).
     */
    size_type degree() const {
      return graph_->edges_[index()].size();
    }
    /** return the beginning position of an iterator to traverse the edges that are connected to this node
     * @pre this node must be a valid node
     * @post (*edge_begin()).node1() == *this (this node)
     */
    incident_iterator edge_begin() const {
      return incident_iterator(graph_,uid_,graph_->edges_[index()].begin());
    }
    /** return iterator referring to the past-the-end element of edges that are connected to this node
     * if incident_iterator ii refers to last edge, then ++i == edge_end()
     * *(edge_end()) is invalid
     */
    incident_iterator edge_end() const {
      return incident_iterator(graph_,uid_,graph_->edges_[index()].end());
    }


   private:
    // Only Graph can access our private members
    friend class Graph;
    // Use this space declare private data members for Node objects
    // Graph needs a way to construct valid Node objects
    graph_type* graph_;   // pointer back to Graph container 
    size_type uid_;  // node's index to identify itself
    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node_value A value associated with the node (defaults to node_value_type() if nothing specified)
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    internal_node new_node = internal_node();
    new_node.position = position;
    new_node.node_value = node_value;
    nodes_.push_back(new_node);   
    uid2idx_.push_back(size()-1);         
    idx2uid_.push_back(next_node_uid_);
    ++next_node_uid_;
    edges_.push_back(vector<internal_edge>());  
    morton_is_valid_ = false;    
    sorted_morton_cont_.push_back(code_uid_pair());
    return Node(this, next_node_uid_-1);
  }

  /** Determine if this Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return this == n.graph_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < size()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes());
    return Node(this,idx2uid_[i]);
  }

  /** Remove a node from the graph.
   * @param[in] n Node to be removed
   * @pre @a n is a valid node of this graph.
   * @post new size() == old size() - 1
   *
   * Can invalidate outstanding iterators. @a n becomes invalid, as do any
   * other Node objects equal to @a n. All other Node objects remain valid.
   *
   * Complexity: Polynomial in size().
   */
  void remove_node(const Node& n) {
    assert(has_node(n));  
    
    // remove edges because of node
    size_type n_idx = n.index();
    size_type n2_idx;
    typename vector<internal_edge>::iterator it, it2;
    it = edges_[n_idx].begin();
    while(it!=edges_[n_idx].end()) {
      n2_idx = uid2idx_[(*it).node2_uid];
      it2 = edges_[n2_idx].begin();
      while(it2!=edges_[n2_idx].end()) { 
        if (uid2idx_[(*it2).node2_uid] == n_idx) {
          it2 = edges_[n2_idx].erase(it2);
          --num_edges_;
        }
        else {
         ++it2;
        } 
      }
      ++it;
    }
    edges_.erase(edges_.begin() + n_idx);
    
    // invalidate Morton code container
    morton_is_valid_ = false;
    sorted_morton_cont_.pop_back();
    
    // remove node
    nodes_.erase(nodes_.begin() + n_idx);
    idx2uid_.erase(idx2uid_.begin() + n_idx);
    for (size_type i = n_idx; i < size(); ++i)
      --uid2idx_[idx2uid_[i]];
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    uid2idx_.clear();
    idx2uid_.clear();
    next_node_uid_ = 0;
    edges_.clear();
    num_edges_ = 0;
    morton_is_valid_ = false;
    sorted_morton_cont_.clear();
  }


  // EDGES

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_,node1_uid_);//graph_->node(graph_->uid2idx_[node1_uid_]);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_,node2_uid_);//graph_->node(graph_->uid2idx_[node2_uid_]);
    }
    
    /** Return the length of this Edge */
    double length() const {
      return norm(node1().position() - node2().position());
    }    

    /** return a reference to the edge value of a valid edge */
    edge_value_type& value() {
      // store value at (smaller_node_uid_, larger_node_uid_)
      size_type smaller_node_uid_ = node1_uid_;
      size_type larger_node_uid_ = node2_uid_;
      if (larger_node_uid_<smaller_node_uid_) {
        smaller_node_uid_ = node2_uid_;
        larger_node_uid_ = node1_uid_;
      }
      typename vector<internal_edge>::iterator it = graph_->edges_[graph_->uid2idx_[smaller_node_uid_]].begin();
      while(it != graph_->edges_[graph_->uid2idx_[smaller_node_uid_]].end() && (*it).node2_uid != larger_node_uid_ )
        ++it;
      return (*it).edge_value;
    }
    /** return a reference to the edge value of a valid edge */
    const edge_value_type& value() const {
      // store value at (smaller_node_uid_, larger_node_uid_)
      size_type smaller_node_uid_ = node1_uid_;
      size_type larger_node_uid_ = node2_uid_;
      if (larger_node_uid_<smaller_node_uid_) {
        smaller_node_uid_ = node2_uid_;
        larger_node_uid_ = node1_uid_;
      }
      typename vector<internal_edge>::iterator it = graph_->edges_[graph_->uid2idx_[smaller_node_uid_]].begin();
      while(it != graph_->edges_[graph_->uid2idx_[smaller_node_uid_]].end() && (*it).node2_uid != larger_node_uid_ )
        ++it;
      return (*it).edge_value;
    }

    /** comparator (==) for the Edge
     * @pre Edge @a e is a valid edge
     * @post returns true iff this edge is part of same Graph as @a e and connects the same node1() and node2() (order is irrelevant)
     */
    bool operator==(const Edge& e) const {
      return (graph_ == e.graph_) && (node1()==e.node1() && node2()==e.node2());
    }
    /** comparator (<) for the Edge
     * @pre Edge @a e is a valid edge
     * @post returns true iff exactly one of the following hold. 
     * (1) this edge is part of same Graph as Edge @a e and has (node1()<e.node1())
     * or ( ( node1()==e.node1() ) and (node2()<e.node2())  ) 
     * (2) the address of this edge's Graph is smaller than that of Edge @a e
     */
    bool operator<(const Edge& e) const {
      return ( (graph_ == e.graph_) && ( (node1()<e.node1()) || (node1()==e.node1() && node2()<e.node2()) ) ) || (graph_ < e.graph_);
    }

   private:
    // Only Graph can access our private members
    friend class Graph;
    // Use this space declare private data members for Edge objects
    graph_type* graph_;   // pointer back to Graph container 
    size_type node1_uid_;  // uid of node1() of edge
    size_type node2_uid_;  // uid of node2() of edge
    /** Private Constructor */
    Edge(const Graph* graph, size_type node1_uid, size_type node2_uid)
        : graph_(const_cast<Graph*>(graph)), node1_uid_(node1_uid), node2_uid_(node2_uid) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {                // deprecated -- use iterators to access edges
    assert(i < num_edges());
    size_type node1_idx = 0;
    size_type edge_count = 0;
    typename vector<internal_edge>::const_iterator it;
    while (node1_idx < size()) {
      size_type node1_uid = idx2uid_[node1_idx];
      for (it = edges_[node1_idx].begin(); it != edges_[node1_idx].end(); ++it) {
        if (node1_uid < (*it).node2_uid) {
          if (edge_count == i)
            return Edge(this,node1_uid,(*it).node2_uid);
          ++edge_count;
        }
      }
      ++node1_idx;
    }
    return Edge(); // never reaches this
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(this == a.graph_ && this == b.graph_);
    
    typename vector<internal_edge>::const_iterator it;
    for (it = edges_[uid2idx_[a.uid_]].begin(); it != edges_[uid2idx_[a.uid_]].end(); ++it)
      if ((*it).node2_uid == b.uid_)
        return true;
    
    return false;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {

    assert( this == a.graph_ && this == b.graph_ && a.index() != b.index() );
  
    size_type node1_uid = min(a.uid_,b.uid_);
    size_type node2_uid = max(a.uid_,b.uid_);
    size_type node1_idx = uid2idx_[node1_uid];
    size_type node2_idx = uid2idx_[node2_uid];
    
    bool already_exists = false;

    typename vector<internal_edge>::const_iterator it;
    for(it = edges_[node1_idx].begin(); it!=edges_[node1_idx].end(); ++it)
      if ((*it).node2_uid == node2_uid)
        already_exists = true;

    if (!already_exists) {
      internal_edge new_edge1 = internal_edge();
      internal_edge new_edge2 = internal_edge();
      new_edge1.node2_uid = node1_uid;
      new_edge2.node2_uid = node2_uid;
      edges_[node1_idx].push_back(new_edge2);
      edges_[node2_idx].push_back(new_edge1); 
      ++num_edges_;
    }

    return Edge(this, a.uid_, b.uid_);

  }

  /** Remove an edge, if any, returning the number of edges removed.
   * @param[in] a,b The nodes potentially defining an edge to be removed.
   * @return 1 if old has_edge(@a a, @a b), 0 otherwise
   * @pre @a a and @a b are valid nodes of this graph
   * @post !has_edge(@a a, @a b)
   * @post new num_edges() == old num_edges() - result
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Can invalidate all edge and incident iterators.
   * Invalidates any edges equal to Edge(@a a, @a b). Must not invalidate
   * other outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Node& a, const Node& b) {
    assert(has_node(a) && has_node(b));
    
    size_type n_removed = 0;
    size_type node1_uid = min(a.uid_,b.uid_);
    size_type node2_uid = max(a.uid_,b.uid_);
    size_type node1_idx = uid2idx_[node1_uid];
    size_type node2_idx = uid2idx_[node2_uid];
    
    typename vector<internal_edge>::iterator it;
    it = edges_[node1_idx].begin();
    while(it!=edges_[node1_idx].end()) {
      if ((*it).node2_uid == node2_uid) {
        it = edges_[node1_idx].erase(it);
        --num_edges_;
        n_removed = 1;
      }
      else {
        ++it;
      }
    }
 
    it = edges_[node2_idx].begin();
    while(it!=edges_[node2_idx].end()) {
      if ((*it).node2_uid == node1_uid) {
        it = edges_[node2_idx].erase(it);
      }
      else {
        ++it;
      }
    } 
    
    return n_removed;
  }

  /** Remove an edge, if any, returning the number of edges removed.
   * @param[in] e The edge to remove
   * @pre @a e is a valid edge of this graph
   * @pre has_edge(@a e.node1(), @a e.node2())
   * @post !has_edge(@a e.node1(), @a e.node2())
   * @post new num_edges() == old num_edges() - 1
   *
   * This is a synonym for remove_edge(@a e.node1(), @a e.node2()), but its
   * implementation can assume that @a e is definitely an edge of the graph.
   * This might allow a faster implementation.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Can invalidate all edge and incident iterators.
   * Invalidates any edges equal to Edge(@a a, @a b). Must not invalidate
   * other outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  // ITERATORS

  /** @class Graph::node_iterator
   * @brief Iterator class for nodes. A forward iterator. */
  class node_iterator : private totally_ordered<node_iterator> {
  
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef ptrdiff_t difference_type;

    /** Construct an invalid node_iterator. */
    node_iterator() {
    }

    /** dereference node iterator, returns the the node object 
     * @pre the node iterator must refer to a valid node
     */
    Node operator*() const {return graph_->node(node_idx_);};
    /** increment the node iterator to refer to the next node in the graph */
    node_iterator& operator++() {++node_idx_; return *this;};
    /** compare whether two node iterators are the same (i.e., refer to the same node) */
    bool operator==(const node_iterator& ni) const {return graph_==ni.graph_ && node_idx_==ni.node_idx_;};

   private:
    friend class Graph;
    friend class edge_iterator;  
    graph_type* graph_;
    size_type node_idx_;
    /** private constructor */
    node_iterator(graph_type* g, size_type node_idx) : graph_(g), node_idx_(node_idx) {
    }
  };

  /** return node iterator that refers to first node in graph
   * *node_end() = node(0)
   */
  node_iterator node_begin() const {
    return node_iterator(const_cast<Graph*>(this),0);
  }
  /** return node iterator that refers to the past-the-end element of the nodes in the graph
  *  incrementing a node_iterator (ni) pointing to the last node goes to node_end()
  *  that is, ++ni==node_end()
  */
  node_iterator node_end() const {
    return node_iterator(const_cast<Graph*>(this),size());
  }


  /** @class Graph::edge_iterator
   * @brief Iterator class for edges. A forward iterator. */
  class edge_iterator : private totally_ordered<edge_iterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef ptrdiff_t difference_type;

    /** Construct an invalid edge_iterator. */
    edge_iterator() {
    }

    /** dereference edge iterator to return edge object */
    Edge operator*() const { return Edge(graph_, graph_->idx2uid_[node1_idx_], (*node2_it_).node2_uid);}
    /** increment edge iterator to refer to the next edge of the graph*/
    edge_iterator& operator++() {
      ++edge_idx_;
      if (edge_idx_ < graph_->num_edges()) {
          ++node2_it_;
          size_type node1_uid = graph_->idx2uid_[node1_idx_];
          while(node2_it_ == graph_->edges_[node1_idx_].end()  || (*node2_it_).node2_uid < node1_uid ) {
            if (node2_it_ == graph_->edges_[node1_idx_].end()) {
              ++node1_idx_;
              node1_uid = graph_->idx2uid_[node1_idx_];
              node2_it_ = graph_->edges_[node1_idx_].begin();
            }
            else {
              ++node2_it_;
            }
          }
      }
      return *this;
    }
    /** test whether two edge iterators are the same (i.e., refer to the same edge)*/
    bool operator==(const edge_iterator& ei) const {return graph_==ei.graph_ && edge_idx_==ei.edge_idx_;}

   private:
    friend class Graph;
    graph_type* graph_;
    size_type edge_idx_; // between 0 and num_edges()
    size_type node1_idx_;
    typename vector<internal_edge>::const_iterator node2_it_;
    /** private constructor */
    edge_iterator(graph_type* graph, size_type edge_idx, size_type node1_idx, typename vector<internal_edge>::const_iterator node2_it)
      : graph_(graph), edge_idx_(edge_idx), node1_idx_(node1_idx), node2_it_(node2_it) {
    }
  };

  /** return an edge iterator that refers to the first edge in the graph
  */
  edge_iterator edge_begin() const {   
    size_type node1_idx = 0;
    // a 'fix' step to find an actual edge in the adjacency list
    typename vector<internal_edge>::const_iterator it = edges_[node1_idx].begin(); 
    if(num_edges()>0) {
      size_type node1_uid = idx2uid_[node1_idx];
      while( it == edges_[node1_idx].end()  || (*it).node2_uid < node1_uid ) {
        if (it == edges_[node1_idx].end()) {
          ++node1_idx;
          node1_uid = idx2uid_[node1_idx];
          it = edges_[node1_idx].begin();
        }
        else {
          ++it;
        }
      }
    }
    return edge_iterator(const_cast<Graph*>(this), 0, node1_idx, it);
  }
  /* returns an iterator referring the past the end of the edges in the graph
  */
  edge_iterator edge_end() const {
    size_type node1_idx = size();
    typename vector<internal_edge>::const_iterator it = edges_[node1_idx].end(); 
    return edge_iterator(const_cast<Graph*>(this), num_edges(), node1_idx, it);
  }


  /** @class Graph::incident_iterator
   * @brief Iterator class for edges incident to a given node. A forward
   * iterator. */
  class incident_iterator : private totally_ordered<incident_iterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef ptrdiff_t difference_type;

    /** Construct an invalid incident_iterator. */
    incident_iterator() {
    }

    /** dereference an incident iterator to return the edge it refers to 
     * @pre incident iterator must point to a valid edge
     * @post for ii an incident iterator, (*ii).node1() == *this (this node)
     */
    Edge operator*() const { return Edge(graph_, node1_uid_, (*node2_it_).node2_uid); };
    incident_iterator& operator++() {++node2_it_; return *this;};
    bool operator==(const incident_iterator& ii ) const {return node2_it_==ii.node2_it_;};

   private:
    friend class Graph;
    graph_type* graph_;
    size_type node1_uid_;
    typename vector<internal_edge>::const_iterator node2_it_;
    /** Private constructor */
    incident_iterator(graph_type* graph,size_type node1_uid,typename vector<internal_edge>::const_iterator node2_it) 
      : graph_(graph), node1_uid_(node1_uid), node2_it_(node2_it) {
    }
  };

 
 
 /** @class Graph::neighborhood_iterator
  * @brief Iterator class for nodes within a bounding box. A forward
  * iterator. */
  class neighborhood_iterator : private totally_ordered<neighborhood_iterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef ptrdiff_t difference_type;
 
    /** Construct an invalid neighborhood_iterator. */
    neighborhood_iterator() {
    }
 
    /** Dereference a neighborhood_iterator to return the node it is pointing to
     * @pre neighborhood_iterator must point to a valid node
     * iterator remains valid as long as the node it refers to is not erased
     */
    Node operator*() const { return Node(graph_,(*mi_).uid); };
    /** increment the iterator to next node in bounding box */
    neighborhood_iterator& operator++() {
      ++mi_;
      // increment further if possible
      if (mi_ != mi_end_) {
        code_type next_code = graph_->morton_coder_.advance_to_box((*mi_).code, min_code_, max_code_);
        while(mi_ != mi_end_ && (*mi_).code < next_code)
          ++mi_;
      }
      return *this;
    }
    /** compare whether two iterators point to the same node */
    bool operator==(const neighborhood_iterator& nbi) const {
      return (mi_ == nbi.mi_);//(graph_ == nbi.graph_) && (mi_ == nbi.mi_) && (bbox_.contains(nbi.bbox_) && nbi.bbox_.contains(bbox_));
    }
 
   private:
    friend class Graph;
    graph_type* graph_;
    typename vector<code_uid_pair>::iterator mi_; // Morton code iterator
    BoundingBox bbox_;
    code_type min_code_;
    code_type max_code_;
    typename vector<code_uid_pair>::iterator mi_end_;
    // private constructor
    neighborhood_iterator(graph_type* g, typename vector<code_uid_pair>::iterator mi, BoundingBox b) 
      : graph_(g), mi_(mi), bbox_(b), min_code_(graph_->morton_coder_.code(bbox_.min())), max_code_(graph_->morton_coder_.code(bbox_.max())) {
    code_uid_pair comp_val = code_uid_pair();
    comp_val.code = max_code_;
    mi_end_ = upper_bound(graph_->sorted_morton_cont_.begin(), graph_->sorted_morton_cont_.end(), comp_val, sorted_morton_cont_comp());
    }
};
 
  /* a helper function to make and sort Morton codes of nodes*/
  void make_morton_codes() {
    for (size_type n = 0; n != num_nodes(); ++n ) {
      sorted_morton_cont_[n].code = morton_coder_.code(node(n).position());
      sorted_morton_cont_[n].uid = idx2uid_[n];
    }
    sort(sorted_morton_cont_.begin(),sorted_morton_cont_.end(), sorted_morton_cont_comp());
    morton_is_valid_ = true;
  }
 
  /* Returns a neighborhood_iterator pointing to the first valid node in the neighborhood
  *@param[in] b is the BoundingBox that defines the neighborhood
  */
  neighborhood_iterator node_begin(const BoundingBox& b) {
    code_type min_code = morton_coder_.code(b.min());
    //code_type max_code = morton_coder_.code(b.max());
    if (!morton_is_valid_)
      make_morton_codes();
    code_uid_pair comp_val = code_uid_pair();
    comp_val.code = min_code;
    typename vector<code_uid_pair>::iterator mi = lower_bound(sorted_morton_cont_.begin(), sorted_morton_cont_.end(), comp_val, sorted_morton_cont_comp());
    return neighborhood_iterator(const_cast<graph_type*>(this), mi, b);
  }
   
  
  /* Returns a neighborhood_iterator pointing to past the last valid node in the neighborhood
  *@param[in] b is the BoundingBox that defines the neighborhood
  */
  neighborhood_iterator node_end(const BoundingBox& b) {
    //code_type min_code = morton_coder_.code(b.min());
    code_type max_code = morton_coder_.code(b.max());
    if (!morton_is_valid_)
      make_morton_codes();
    code_uid_pair comp_val = code_uid_pair();
    comp_val.code = max_code;
    typename vector<code_uid_pair>::iterator mi = upper_bound(sorted_morton_cont_.begin(), sorted_morton_cont_.end(), comp_val, sorted_morton_cont_comp());
    return neighborhood_iterator(const_cast<graph_type*>(this), mi, b);
  }

  
 private:

  // Use this space for your Graph class's internals:
  //   helper functions you might need, data members, and so forth.
  
  // internal type for nodes in graph
  struct internal_node {
    Point position;
    node_value_type node_value;
  };
  
  // internal type for edges in graph
  struct internal_edge {
    size_type node2_uid;
    edge_value_type edge_value;
  };
  
  vector<internal_node> nodes_; // vector of nodes in graph. Indexed by node index(). Size of vector is num_edges()
  vector<size_type> uid2idx_;   // map of node uids to indices
  vector<size_type> idx2uid_;   // map of node indices to uids
  size_type next_node_uid_;
  
  // edges in graph stored as adjacency list of node uids
  vector<vector<internal_edge>> edges_; // vector of vector of edge data (node2 uids and edge values) indexed by of node1 index()
  size_type num_edges_;
  
  // Morton Coder 
  struct code_uid_pair {
    code_type code;
    size_type uid;
  };
  BoundingBox domain_bbox_;
  MortonCoder<5> morton_coder_;
  bool morton_is_valid_; // bool for sorted_morton_cont is valid
  vector<code_uid_pair>  sorted_morton_cont_; // construct container on demand if invalid
  struct sorted_morton_cont_comp {
    bool operator()(code_uid_pair p1, code_uid_pair p2) const {
      return (p1.code < p2.code);
    }
  };
  
};

#endif
