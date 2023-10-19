#ifndef CS207_MESH_HPP
#define CS207_MESH_HPP

/** @file Mesh.hpp
 * @brief A triangle mesh type
 */

#include "CS207/Util.hpp"
#include "Point.hpp"
#include "Graph.hpp"
#include <math.h>       // fabs,sqrt

#include <algorithm>
#include <vector>
#include <cassert>

using namespace std;

/** @class Mesh
 * @brief A template for 2D triangle meshes.
 *
 * Users can add and retrieve nodes and triangles.
 * V is user-defined node value type
 * E is user-defined edge value type
 * T is user-defined triangle value type
 */
template<typename V, typename E, typename T>
class Mesh{

 private:

  // Use this space for declarations of important internal types you need
  // later in the Mesh's definition.
  
  struct internal_node;
  struct internal_edge;
  struct internal_triangle;
  struct internal_triangle_pair;  
  
  typedef Graph<internal_node,internal_edge> primal_graph_type;
  typedef Graph<internal_triangle,internal_triangle_pair> dual_graph_type;
  typedef typename primal_graph_type::node_type pg_node_type;
  typedef typename primal_graph_type::edge_type pg_edge_type;
  typedef typename dual_graph_type::node_type dg_node_type;
  typedef typename dual_graph_type::edge_type dg_edge_type;
  typedef typename primal_graph_type::node_iterator pg_node_iterator_type;
  typedef typename primal_graph_type::edge_iterator pg_edge_iterator_type;
  typedef typename primal_graph_type::incident_iterator pg_incident_iterator_type;
  typedef typename dual_graph_type::node_iterator dg_node_iterator_type;
  typedef typename dual_graph_type::edge_iterator dg_edge_iterator_type;  
  
 public:

  // PUBLIC TYPE DEFINITIONS 
  
  /** Type of the node value. */
  typedef V node_value_type;
 
  /** Type of the edge value. */
  typedef E edge_value_type; 
 
  /** Type of the triangle value. */
  typedef T triangle_value_type;
  
  /** Type of this mesh. */
  typedef Mesh mesh_type;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;

  /** Predeclaration of Triangle type. */
  class Triangle;
  /** Synonym for Triangle (following STL conventions). */
  typedef Triangle triangle_type;  
  
  /** Type of indexes and sizes. */
  typedef unsigned size_type;
  
  /** Type of node iterators, which iterate over all mesh nodes. */
  class node_iterator;

  /** Type of edge iterators, which iterate over all mesh edges. */
  class edge_iterator;

  /** Type of triangle iterators, which iterate over all mesh triangles. */
  class triangle_iterator;
  
  /** Type of incident iterators, which iterate incident triangles to a node. */
  class incident_iterator;

  /** Type of edge incident iterators, which iterate incident edges to a node. */
  class edge_incident_iterator;
  
  // CONSTRUCTOR AND DESTRUCTOR

  /** Construct an empty Mesh. 
   * primal_graph_ is a Graph of triangle vertices and edges
   * dual_graph_ is a Graph of triangles and triangle connectivity
   */
  Mesh() 
    : primal_graph_(), dual_graph_() {
  }
  /** Default destructor */
  ~Mesh() = default;
  

  // NODES

  /** @class Mesh::Node
   * @brief Class representing the mesh's nodes.
   *
   * Node objects are used to access information about the mesh's nodes.
   */
  class Node : private totally_ordered<Node> {
   public:
    /** Construct an invalid node.*/
    Node() {
    }

    /** Return this node's position. */
    Point position() const {
      return pgn_.position();
    }

    /** Return this node's index, a number in the range [0, num_nodes()). */
    size_type index() const {
      return pgn_.index();
    }
    
    /** Return this node's degree. */
    size_type degree() const {
      return pgn_.degree();
    }

    /** return a reference to the node value of a valid node*/
    node_value_type& value() {
      return pgn_.value().node_value;
    }
    /** return a reference to the node value of a valid node*/
    const node_value_type& value() const {
      return pgn_.value().node_value;
    }
    /** comparator (==) for the Node
     * @pre Node @a n is a valid node
     * @post returns true iff this node is part of same Mesh as n and has same index()
     */
    bool operator==(const Node& n) const {
      return pgn_ == n.pgn_;
    }
    /** comparator (<) for the Node
     * @pre Node n is a valid node
     * @post returns true iff this node is part of same Mesh as n and has smaller index() than Node n
     * or if address of this node's mesh is smaller than that of Node n
     */
    bool operator< (const Node& n) const {
      return pgn_ < n.pgn_;
    }

    /** return the beginning position of an iterator to traverse the triangles that are adjacent to this node
     * @pre this node must be a valid node
     */
    incident_iterator triangle_begin() const {
      pg_incident_iterator_type pgii_ = pgn_.edge_begin();
      while( pgii_ != (*pgii_).node1().edge_end() && 
            ((*pgii_).value().num_triangles == 1) && 
            (cross((*pgii_).node2().position()-(*pgii_).node1().position(),
                    (*pgii_).value().dgn1.position()-(*pgii_).node1().position()).z <= 0)      ) {
        ++pgii_; 
      }
      return incident_iterator(mesh_,pgii_);
    }
    /** return iterator referring to the past-the-end element of triangles that are adjacent to this node
     */
    incident_iterator triangle_end() const {
      return incident_iterator(mesh_,pgn_.edge_end());
    }
    
    /** return the beginning position of an iterator to traverse the edges that are adjacent to this node
     * @pre this node must be a valid node
     */
    edge_incident_iterator edge_begin() const {
      return edge_incident_iterator(mesh_,pgn_.edge_begin());
    }
    /** return iterator referring to the past-the-end element of edges that are adjacent to this node
     */
    edge_incident_iterator edge_end() const {
      return edge_incident_iterator(mesh_,pgn_.edge_end());
    }    


   private:
    friend class Mesh;
    mesh_type* mesh_;   // pointer back to Mesh container 
    pg_node_type pgn_;
    /** Private Constructor */
    Node(const Mesh* mesh, pg_node_type pgn)
        : mesh_(const_cast<Mesh*>(mesh)), pgn_(pgn) {
    }
  };
  
  /** Add a node to the mesh, returning the added node.
   * @param[in] position The new node's position
   * @param[in] node_value A value associated with the node (defaults to node_value_type() if nothing specified)
   * @post new size() == old size() + 1
   * @post result_node.index() == old size()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
    internal_node new_node = internal_node();
    new_node.node_value = node_value;
    pg_node_type pgn = primal_graph_.add_node(position,new_node);
    return Node(this, pgn);
  }
  
  /** Return the number of nodes in the mesh.
   *
   * Complexity: O(1).
   */
  size_type num_nodes() const {
    return primal_graph_.num_nodes();
  }
  
  // EDGES

  /** @class Mesh::Edge
   * @brief Class representing the mesh's edges.
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
      return Node(mesh_,pge_.node1());
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(mesh_,pge_.node2());
    }
    
    /** Return the length of this Edge */
    double length() const {
      return norm(node1().position() - node2().position());
    } 
    
    /** Return the center of mass of this Edge */
    Point center_of_mass() const {
      return 0.5*(node1().position() + node2().position());
    } 
    
    /** return a reference to the node value of a valid node*/
    edge_value_type& value() {
      return pge_.value().edge_value;
    }
    /** return a reference to the node value of a valid node*/
    const edge_value_type& value() const {
      return pge_.value().edge_value;
    }
    
    /** Return a triangle of this Edge */
    Triangle triangle1() const {
      return Triangle(mesh_,pge_.value().dgn1);
    }

    /** Return the other triangle of this Edge */
    Triangle triangle2() const {
      return Triangle(mesh_,pge_.value().dgn2);
    }
    
    /** Return the number of triangles */
    size_type num_triangles() const {
      return pge_.value().num_triangles;
    }
    
    /** Return the normal vector from triangle1() to triangle2() */
    Point norm_vector_12() const {
      Point nv;
      nv = cross(node2().position()-node1().position(),triangle1().center_of_mass()-node1().position());
      nv = cross(node2().position()-node1().position(),nv);
      nv /= norm(nv);  // normalize
      nv *= pge_.length();
      return nv;
    } 
    
    /** comparator (==) for the Edge
     * @pre Edge @a e is a valid edge
     * @post returns true iff this edge is part of same Mesh as @a e and connects the same triangle1() and triangle2() (order is irrelevant)
     */
    bool operator==(const Edge& e) const {
      return pge_==e.pge_;
    }
    /** comparator (<) for the Edge
     * @pre Edge @a e is a valid edge
     * @post returns true iff exactly one of the following hold. 
     * (1) this edge is part of same Mesh as Edge @a e and has (node1()<e.node1())
     * or ( ( node1()==e.node1() ) and (node2()<e.node2())  ) 
     * (2) the address of this edge's Mesh is smaller than that of Edge @a e
     */
    bool operator<(const Edge& e) const {
      return pge_<e.pge_;
    }

   private:
    friend class Mesh;
    mesh_type* mesh_;   // pointer back to Mesh container 
    pg_edge_type pge_;
    /** Private Constructor */
    Edge(const Mesh* mesh, pg_edge_type pge)
        : mesh_(const_cast<Mesh*>(mesh)), pge_(pge) {
    }
  };
  
  /** Return the number of edges in the mesh.
   *
   * Complexity: O(1).
   */
  size_type num_edges() const {
    return primal_graph_.num_edges();
  }
  

  // TRIANGLES

  /** @class Mesh::Triangle
   * @brief Class representing the mesh's nodes.
   *
   * Triangle objects are used to access information about the mesh's triangles.
   */
  class Triangle : private totally_ordered<Triangle> {
   public:
    /** Construct an invalid node.*/
    Triangle() {
    }
    
    /** Return this triangle's index, a number in the range [0, num_triangles()). */
    size_type index() const {
      return dgn_.index();
    }
    
    /** return a reference to the triangle value of a valid triangle*/
    triangle_value_type& value() {
      return dgn_.value().triangle_value;
    }
    /** return a reference to the triangle value of a valid triangle*/
    const triangle_value_type& value() const {
      return dgn_.value().triangle_value;
    }

    /** Return node i of this Triangle 
     * @pre i \in {1,2,3}
     */
    Node node(size_type i) const {
      assert(1 <= i && i <= 3);
      if(i == 2) {
        return Node(mesh_,dgn_.value().pgn2);
      }
      else if(i == 3) {
        return Node(mesh_,dgn_.value().pgn3);
      }
      else {
        return Node(mesh_,dgn_.value().pgn1);
      }
    } 
 
    /** Return edge i of this Triangle 
     * @pre i \in {1,2,3}
     */
    Edge edge(size_type i) const {
      assert(1 <= i && i <= 3);
      if(i == 2) {
        return Edge(mesh_,mesh_->primal_graph_.add_edge(dgn_.value().pgn3,dgn_.value().pgn1));
      }
      else if(i == 3) {
        return Edge(mesh_,mesh_->primal_graph_.add_edge(dgn_.value().pgn1,dgn_.value().pgn2));
      }
      else {      
        return Edge(mesh_,mesh_->primal_graph_.add_edge(dgn_.value().pgn2,dgn_.value().pgn3));
      }
    } 
    
    /** Return the center of mass of this Triangle */
    Point center_of_mass() const {
      return dgn_.position();
    } 
    
    /** Return the outward normal vector of edge i of this Triangle */
    Point norm_vector(size_type i) const {
      Point nv;
      assert(1 <= i && i <= 3);
      if(i == 1) {
        nv = node(3).position()-node(2).position();
      }
      else if(i == 2) {
        nv = node(1).position()-node(3).position();
      }
      else {
        nv = node(2).position()-node(1).position();
      }
      nv = Point(nv.y,-nv.x,0.0); // rotate 90 deg clockwise (y,-x)
      //nv /= norm(nv);  // normalize
      return nv;
    } 
    
    /** Return the area of this Triangle */
    double area() const {
      return 0.5*fabs(cross(node(3).position()-node(1).position(),node(2).position()-node(1).position()).z);
    }  
    
    /** comparator (==) for the Triangle
     * @pre Triangle @a t is a valid triangle
     * @post returns true iff this triangle is part of same Mesh as t and has same index()
     */
    bool operator==(const Triangle& t) const {
      return dgn_==t.dgn_;
    }
    /** comparator (<) for the Triangle
     * @pre Triangle t is a valid triangle
     * @post returns true iff this triangle is part of same Mesh as t and has smaller index() than Triangle t
     * or if address of this triangle's mesh is smaller than that of Triangle t
     */
    bool operator< (const Triangle& t) const {
      return dgn_<t.dgn_;
    }

   private:
    friend class Mesh;
    mesh_type* mesh_;   // pointer back to Mesh container 
    dg_node_type dgn_;  // dual_graph_ node
    /** Private Constructor */
    Triangle(const Mesh* mesh, dg_node_type dgn)
        : mesh_(const_cast<Mesh*>(mesh)), dgn_(dgn) {
    }
  };
  
  /** Check whether Mesh has triangle defined by 3 nodes
  */
  bool has_triangle(const Node& n1, const Node& n2, const Node& n3) {
    bool already_exists = false;
    // XXX implement if necessary, but for this HW we can give the precondition that add_triangle does not add an existing triangle
    return already_exists;
  }  
  
  /** Add a triangle to the mesh, returning the added triangle.
   * @param[in] n1,n2,n3 The nodes of the triangle
   * @param[in] triangle_value A value associated with the triangle (defaults to triangle_value_type() if nothing specified)
   * @pre n1,n2,n3 are distinct
   * @pre triangle has not been added before (XXX may wish to change this by implementing has_triangle)
   * @post 
   *
   * Complexity: O(1) amortized operations.
   */
  Triangle add_triangle(Node& n1, Node& n2, Node& n3, const triangle_value_type& triangle_value = triangle_value_type()) {
    assert(n1 != n2 && n1 != n3 && n2 != n3);
    //if has_triangle(const Node& n1, const Node& n2, const Node& n3) //  XXX no need to do has_triangle. use precondition
    //  return Triangle(this, /*XXX*/);
    internal_triangle new_triangle = internal_triangle();
    new_triangle.triangle_value = triangle_value;
    // find lowest node index to be node 1
    if(n1.index() < n2.index() && n1.index() < n3.index()) {
      new_triangle.pgn1 = n1.pgn_;
      new_triangle.pgn2 = n2.pgn_;
      new_triangle.pgn3 = n3.pgn_;
    }
    else if(n2.index() < n3.index()) {
      new_triangle.pgn1 = n2.pgn_;
      new_triangle.pgn2 = n3.pgn_;
      new_triangle.pgn3 = n1.pgn_;
    }
    else {
      new_triangle.pgn1 = n3.pgn_;
      new_triangle.pgn2 = n1.pgn_;
      new_triangle.pgn3 = n2.pgn_;
    }
    // orient nodes counter-clockwise (2-1 cross 3-2 should be positive, else switch nodes 2 and 3)
    if (cross(new_triangle.pgn2.position()-new_triangle.pgn1.position(),
              new_triangle.pgn3.position()-new_triangle.pgn1.position()).z < 0) {
      pg_node_type temp = new_triangle.pgn2;
      new_triangle.pgn2 = new_triangle.pgn3;
      new_triangle.pgn3 = temp;
    }
    // add node to dual graph
    auto dgn = dual_graph_.add_node((n1.position()+n2.position()+n3.position())/3.0,new_triangle); // triangle node position is CM
    // add edges to primal graph and add edges to dual graph
    pg_edge_type pge;   
    pge = primal_graph_.add_edge(n1.pgn_,n2.pgn_);
    if (pge.value().num_triangles == 0) {
      pge.value().dgn1 = dgn;
   }
    else {
      pge.value().dgn2 = dgn;
      dual_graph_.add_edge(pge.value().dgn1,pge.value().dgn2);
    }
    ++pge.value().num_triangles;
    
    pge = primal_graph_.add_edge(n2.pgn_,n3.pgn_);
    if (pge.value().num_triangles == 0) {
      pge.value().dgn1 = dgn;
   }
    else {
      pge.value().dgn2 = dgn;
      dual_graph_.add_edge(pge.value().dgn1,pge.value().dgn2);
    }
    ++pge.value().num_triangles;
    
    pge = primal_graph_.add_edge(n3.pgn_,n1.pgn_);
    if (pge.value().num_triangles == 0) {
      pge.value().dgn1 = dgn;
   }
    else {
      pge.value().dgn2 = dgn;
      dual_graph_.add_edge(pge.value().dgn1,pge.value().dgn2);
    }
    ++pge.value().num_triangles;
    
    
    return Triangle(this, dgn);
  }
  
  /** Return the number of triangles in the mesh.
   *
   * Complexity: O(1).
   */
  size_type num_triangles() const {
    return dual_graph_.num_nodes();
  }  
  
  
 
  // ITERATORS

  /** @class Mesh::node_iterator
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
    Node operator*() const { return Node(mesh_,*pgni_); };
    /** increment the node iterator to refer to the next node in the mesh */
    node_iterator& operator++() {++pgni_; return *this;};
    /** compare whether two node iterators are the same (i.e., refer to the same node) */
    bool operator==(const node_iterator& ni) const {return mesh_==ni.mesh_ && pgni_ == ni.pgni_;};

   private:
    friend class Mesh;
    mesh_type* mesh_;
    pg_node_iterator_type pgni_;
    /** private constructor */
    node_iterator(mesh_type* mesh, pg_node_iterator_type pgni) : mesh_(const_cast<Mesh*>(mesh)), pgni_(pgni) {
    }
  };

  /** return node iterator that refers to first node in mesh */
  node_iterator node_begin() const {
    return node_iterator(const_cast<Mesh*>(this),primal_graph_.node_begin());
  }
  /** return node iterator that refers to the past-the-end element of the nodes in the mesh */
  node_iterator node_end() const {
    return node_iterator(const_cast<Mesh*>(this),primal_graph_.node_end());
  }

 
  /** @class Mesh::edge_iterator
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

    /** dereference edge iterator, returns the the edge object 
     * @pre the edge iterator must refer to a valid edge
     */
    Edge operator*() const { return Edge(mesh_,*pgei_); };
    /** increment the edge iterator to refer to the next edge in the mesh */
    edge_iterator& operator++() {++pgei_; return *this;};
    /** compare whether two edge iterators are the same (i.e., refer to the same edge) */
    bool operator==(const edge_iterator& ei) const {return mesh_==ei.mesh_ && pgei_ == ei.pgei_;};

   private:
    friend class Mesh;
    mesh_type* mesh_;
    pg_edge_iterator_type pgei_;
    /** private constructor */
    edge_iterator(mesh_type* mesh, pg_edge_iterator_type pgei) : mesh_(const_cast<Mesh*>(mesh)), pgei_(pgei) {
    }
  };

  /** return edge iterator that refers to first edge in mesh */
  edge_iterator edge_begin() const {
    return edge_iterator(const_cast<Mesh*>(this),primal_graph_.edge_begin());
  }
  /** return edge iterator that refers to the past-the-end element of the edges in the mesh */
  edge_iterator edge_end() const {
    return edge_iterator(const_cast<Mesh*>(this),primal_graph_.edge_end());
  }
  

  /** @class Mesh::triangle_iterator
   * @brief Iterator class for triangles. A forward iterator. */
  class triangle_iterator : private totally_ordered<triangle_iterator> {
  
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef ptrdiff_t difference_type;

    /** Construct an invalid triangle_iterator. */
    triangle_iterator() {
    }

    /** dereference triangle iterator, returns the the triangle object 
     * @pre the triangle iterator must refer to a valid triangle
     */
    Triangle operator*() const { return Triangle(mesh_,*dgni_); };
    /** increment the triangle iterator to refer to the next triangle in the mesh */
    triangle_iterator& operator++() {++dgni_; return *this;};
    /** compare whether two triangle iterators are the same (i.e., refer to the same triangle) */
    bool operator==(const triangle_iterator& ti) const {return mesh_==ti.mesh_ && dgni_ == ti.dgni_;}; 

   private:
    friend class Mesh;
    mesh_type* mesh_;
    dg_node_iterator_type dgni_;
    /** private constructor */
    triangle_iterator(mesh_type* mesh, dg_node_iterator_type dgni) : mesh_(const_cast<Mesh*>(mesh)), dgni_(dgni) {
    }
  };

  /** return triangle iterator that refers to first triangle in mesh */
  triangle_iterator triangle_begin() const {
    return triangle_iterator(const_cast<Mesh*>(this),dual_graph_.node_begin());
  }
  /** return triangle iterator that refers to the past-the-end element of the triangles in the mesh */
  triangle_iterator triangle_end() const {
    return triangle_iterator(const_cast<Mesh*>(this),dual_graph_.node_end());
  }
 
 
  /** @class Mesh::incident_iterator
   * @brief Iterator class for triangles incident to a given node. A forward
   * iterator. */
  class incident_iterator : private totally_ordered<incident_iterator> {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Triangle value_type;
    /** Type of pointers to elements. */
    typedef Triangle* pointer;
    /** Type of references to elements. */
    typedef Triangle& reference;
    /** Iterator category. */
    typedef input_iterator_tag iterator_category;
    /** Difference between iterators */
    typedef ptrdiff_t difference_type;

    /** Construct an invalid incident_iterator. */
    incident_iterator() {
    }

    /** dereference an incident iterator to return the triangle it refers to 
     * @pre incident iterator must point to a valid triangle
     */
    Triangle operator*() const { 
      if (cross((*pgii_).node2().position()-(*pgii_).node1().position(),
            (*pgii_).value().dgn1.position()-(*pgii_).node1().position()).z > 0) {
        return Triangle(mesh_,(*pgii_).value().dgn1);  
      } 
      else {
        return Triangle(mesh_,(*pgii_).value().dgn2);        
      }
    };
    incident_iterator& operator++() {
      ++pgii_; 
      while( pgii_ != (*pgii_).node1().edge_end() && 
            ((*pgii_).value().num_triangles == 1) && 
            (cross((*pgii_).node2().position()-(*pgii_).node1().position(),
                    (*pgii_).value().dgn1.position()-(*pgii_).node1().position()).z <= 0)      ) {
        ++pgii_; 
      }
      return *this;
    };
    bool operator==(const incident_iterator& ii ) const {return mesh_==ii.mesh_ && pgii_ == ii.pgii_;};

   private:
    friend class Mesh;
    mesh_type* mesh_;
    pg_incident_iterator_type pgii_;
    /** private constructor */
    incident_iterator(mesh_type* mesh, pg_incident_iterator_type pgii) : mesh_(const_cast<Mesh*>(mesh)), pgii_(pgii) {
    }
  };
  
 
 
  /** @class Mesh::edge_incident_iterator
   * @brief Iterator class for edges incident to a given node. A forward
   * iterator. */
  class edge_incident_iterator : private totally_ordered<edge_incident_iterator> {
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

    /** Construct an invalid edge_incident_iterator. */
    edge_incident_iterator() {
    }

    /** dereference an incident iterator to return the edge it refers to 
     * @pre incident iterator must point to a valid triangle
     */
    Edge operator*() const { 
      return Edge(mesh_,*pgii_);
    };
    edge_incident_iterator& operator++() {
      ++pgii_; 
      return *this;
    };
    bool operator==(const edge_incident_iterator& ii ) const {return mesh_==ii.mesh_ && pgii_ == ii.pgii_;};

   private:
    friend class Mesh;
    mesh_type* mesh_;
    pg_incident_iterator_type pgii_;
    /** private constructor */
    edge_incident_iterator(mesh_type* mesh, pg_incident_iterator_type pgii) : mesh_(const_cast<Mesh*>(mesh)), pgii_(pgii) {
    }
  };
   
 
 
 
 private:
  // Use this space for your Mesh class's internals:
  //   helper functions you might need, data members, and so forth.
  friend class Node;
  // internal type for nodes in mesh
  struct internal_node {
    node_value_type node_value;
  };
 
  // internal type for edge in mesh
  struct internal_edge {
    edge_value_type edge_value;
    dg_node_type dgn1;
    dg_node_type dgn2;  // represents 2 triangles adjacent to edge
    size_type num_triangles;
  };  
  
  // internal type for triangles in mesh
  struct internal_triangle {
    triangle_value_type triangle_value;
    pg_node_type pgn1;
    pg_node_type pgn2;
    pg_node_type pgn3; // 3 nodes of triangle, oriented counter-clockwise, pgn1 has smallest index
    };
  
  // internal type for triangle pairs in mesh
  struct internal_triangle_pair {
  }; 
  
  primal_graph_type primal_graph_; // graph for nodes and edges and data associated with them
  dual_graph_type dual_graph_;     // graph for triangle (each node is a triangle, edges are triangle connectivity)
 
};

#endif