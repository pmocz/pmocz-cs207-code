/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include <fstream>

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"
#include "MortonCoder.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

// data held by a node
struct node_data {
  double m;  // mass
  Point v;   // velocity
  double c;  // damping coefficient
};

// data held by an edge
struct edge_data {
  double k;  // spring constant
  double l;  // spring rest length
};

typedef Graph<node_data,edge_data> GraphType;
typedef GraphType::Node Node;
typedef GraphType::Edge Edge;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g Graph
 * @param[in] t The current time (useful for time-dependent forces)
 * @param[in] dt The time step
 * @param[in] force Function object defining the force per node
 * @param[in] constrain Functor object enforcing constraints
 * @pre G::node_value_type supports node_data structure (has members m and v)
 * @return the next time step (usually @a t + @a dt)
 *
 * @a force is called as @a force(n, @a t), where n is a node of the graph
 * and @a t is the current time parameter. @a force must return a Point
 * representing the node's force at time @a t.
 * @a constraint takes the graph @a g and the time @a t, searches through 
 * nodes in the graph that violate some constraints, and resets their positions
 * and velocities before the forces are calculated
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {

  for(typename G::node_iterator ni=g.node_begin(); ni!=g.node_end(); ++ni)
    (*ni).set_position( (*ni).position() + dt * (*ni).value().v);
  
  constraint(g,t);
  
  for(typename G::node_iterator ni=g.node_begin(); ni!=g.node_end(); ++ni)
    (*ni).value().v += force(*ni,t) * dt / (*ni).value().m;
  
  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;     // silence compiler warnings
    
    Point f = Point();
    f.z = -grav*n.value().m;
    
    Point xi = n.position();
    for (GraphType::incident_iterator ii = n.edge_begin(); ii!= n.edge_end(); ++ii ) {
      Point xj = (*ii).node2().position();
      f += -(*ii).value().k * ((xi-xj)/norm(xi-xj))*(norm(xi-xj) - (*ii).value().l);
    }
    
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0) )
      f = Point(0,0,0);
    
    return f;
  } 
};


/** Gravity Force */
struct GravityForce {
  /** Return the gravity force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;     // silence compiler warnings  
    Point f = Point();
    f.z = -grav*n.value().m;
    return f;
  } 
};


/** Mass Spring Force */
struct MassSpringForce {
  /** Return the gravity force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;     // silence compiler warnings  
    Point f = Point();
    Point xi = n.position();
    for (GraphType::incident_iterator ii = n.edge_begin(); ii!= n.edge_end(); ++ii ) {
      Point xj = (*ii).node2().position();
      f += -(*ii).value().k * ((xi-xj)/norm(xi-xj))*(norm(xi-xj) - (*ii).value().l);
    }
    return f;
  } 
};


/** Damping Force */
struct DampingForce {
  /** Return the damping force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;     // silence compiler warnings  
    return -n.value().c*n.value().v;
  } 
};


/** Null Force */
struct NullForce {
  /** Return 0 force . */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n, (void) t;     // silence compiler warnings  
    return Point();
  } 
};


/** Combined Force */
template <typename F1, typename F2>
struct CombinedForce {
  /** Return the combined force applying to @a n at time @a t. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n,t) + f2_(n,t);
  } 
  // constructor 2 forces
  CombinedForce(const F1 f1, const F2 f2) : f1_(f1), f2_(f2) {
  };
  private:
    F1 f1_;
    F2 f2_;
};

/** function that combines 2 forces*/
template <typename F1, typename F2>
CombinedForce<F1,F2> make_combined_force(F1 force1, F2 force2) {
  return CombinedForce<F1,F2>(force1,force2);
}

/** function that combines 3 forces*/
template <typename F1, typename F2, typename F3>
CombinedForce<CombinedForce<F1,F2>,F3> make_combined_force(F1 force1, F2 force2, F3 force3) {
  return CombinedForce<CombinedForce<F1,F2>,F3>(make_combined_force(force1,force2),force3);
}

/** fixed nodes Constraint */
struct ConstantNodeConstraint {
  /**  Fix the points at (0, 0, 0) and (1, 0, 0) in the graph. */

  void operator()(GraphType& g, double t) {
    (void)g, (void) t;     // silence compiler warnings   
    for(unsigned i = 0; i< fixed_nodes000_.size(); ++i) {
      (*fixed_nodes000_[i]).set_position(Point(0,0,0));
      (*fixed_nodes000_[i]).value().v == Point(0,0,0); 
    }
    for(unsigned i = 0; i< fixed_nodes100_.size(); ++i) {
      (*fixed_nodes100_[i]).set_position(Point(1,0,0));
      (*fixed_nodes100_[i]).value().v == Point(0,0,0); 
    }
  }
    
  // constructor
  ConstantNodeConstraint(GraphType& g) {
    fixed_nodes000_ = std::vector<GraphType::node_iterator>();
    fixed_nodes100_ = std::vector<GraphType::node_iterator>();
    for(GraphType::node_iterator ni=g.node_begin(); ni!=g.node_end(); ++ni) {
      if ((*ni).position() == Point(0,0,0))
        fixed_nodes000_.push_back(ni);
      if ((*ni).position() == Point(1,0,0))
        fixed_nodes100_.push_back(ni);
    }
  }
  
  private:
    std::vector<GraphType::node_iterator> fixed_nodes000_; 
    std::vector<GraphType::node_iterator> fixed_nodes100_;
};

/** Null Constraint */
struct NullConstraint {
  /**  do nothing */
  void operator()(GraphType& g, double t) {
    (void) g, (void) t;     // silence compiler warnings   
  }
};

/** Plane Constraint */
struct PlaneConstraint {
  /**  plane at z = ..., initialized by constructor 
   * reset position to nearest point on surface of plane
   * reset normal component of velocity to 0
  */
  void operator()(GraphType& g, double t) {
    (void) t;     // silence compiler warnings   
    BoundingBox b = BoundingBox(Point(-5,-5,plane_z_-0.1),Point(5,5,plane_z_));
    GraphType::neighborhood_iterator end = g.node_end(b);
    //for(GraphType::node_iterator ni=g.node_begin(); ni!=g.node_end(); ++ni)
    for(GraphType::neighborhood_iterator ni=g.node_begin(b); ni!=end; ++ni)
      if ((*ni).position().z < plane_z_) {
        Point new_pos = (*ni).position();
        new_pos.z = plane_z_;
        (*ni).value().v.z = 0;
        (*ni).set_position(new_pos);
      }
  }
 
  //constructor
  PlaneConstraint(double plane_z) : plane_z_(plane_z) {
  }
  
  private:
    double plane_z_;
};


/** Sphere Constraint */
struct SphereConstraint {
  /**  sphere with center and radius, initialized by constructor 
   * reset position to nearest point on surface of sphere
   * reset normal component of velocity to 0
   */
  void operator()(GraphType& g, double t) {
    (void) t;     // silence compiler warnings   
    BoundingBox b = BoundingBox(c_-Point(r_,r_,r_),c_+Point(r_,r_,r_));
    GraphType::neighborhood_iterator end = g.node_end(b);
    //for(GraphType::node_iterator ni=g.node_begin(); ni!=g.node_end(); ++ni)
    for(GraphType::neighborhood_iterator ni=g.node_begin(b); ni!=end; ++ni) {
      Point r_from_ctr = (*ni).position()-c_;
      double d_from_ctr = norm(r_from_ctr);
      if ( d_from_ctr < r_ && (*ni).index() != ctr_node_idx_) {
        Point new_pos = r_ * r_from_ctr/d_from_ctr + c_;
        Point Ri = r_from_ctr/d_from_ctr;
        (*ni).value().v -= dot((*ni).value().v,Ri) * Ri;
        (*ni).set_position(new_pos);
      }
    }
  }
  
  //constructor
  SphereConstraint(Point center, double radius) : c_(center), r_(radius) {
    ctr_node_idx_ = -1;
  }
  SphereConstraint(Point center, double radius, unsigned int ctr_node_idx) : c_(center), r_(radius), ctr_node_idx_(ctr_node_idx) {
  }
  
  private:
    Point c_;
    double r_;
    unsigned int ctr_node_idx_;
};

/** Avoid Self-Collision Constraint*/
struct AvoidSelfCollisionConstraint {
  /**  avoid self-collisions by treating each node as a sphere or radius = minimum of connected edge lengths */
  void operator()(GraphType& g, double t) {
    for(GraphType::node_iterator ni=g.node_begin(); ni!=g.node_end(); ++ni) {
      //find radius of node
      double radius = 0;
      GraphType::incident_iterator ii = (*ni).edge_begin();
      if (ii!= (*ni).edge_end()) {
        radius = (*ii).length();
        ++ii;
      }
      while( ii != (*ni).edge_end() ) {
        if ((*ii).length() < radius)
          radius = (*ii).length();
        ++ii;
      }
      radius *= 0.99; // make threshold radius slightly smaller than minimum edge length
      //apply spherical constraint
      SphereConstraint sc = SphereConstraint((*ni).position(), radius, (*ni).index());
      sc(g,t);
    }
  }
};


/** Combined Constraint */
template <typename C1, typename C2, typename C3=NullConstraint, typename C4=NullConstraint>
struct CombinedConstraint {
  /** Return the combined constraint applying to @a g at time @a t. */
  void operator()(GraphType& g, double t) {
    c1_(g,t);
    c2_(g,t);
  } 
  // constructor 2 constraints
  CombinedConstraint(const C1 c1, const C2 c2) : c1_(c1), c2_(c2) {
  };
  private:
    C1 c1_;
    C2 c2_;
};


/** function that combines 2 constraints*/
template <typename C1, typename C2>
CombinedConstraint<C1,C2> make_combined_constraint(C1 constraint1, C2 constraint2) {
  return CombinedConstraint<C1,C2>(constraint1,constraint2);
}

/** function that combines 3 constraints*/
template <typename C1, typename C2, typename C3>
CombinedConstraint<CombinedConstraint<C1,C2>,C3> make_combined_constraint(C1 c1, C2 c2, C3 c3) {
  return CombinedConstraint<CombinedConstraint<C1,C2>,C3>(make_combined_constraint(c1,c2),c3);
}

/** function that combines 4 constraints*/
template <typename C1, typename C2, typename C3, typename C4>
CombinedConstraint<CombinedConstraint<C1,C2>,CombinedConstraint<C3,C4>> make_combined_constraint(C1 c1, C2 c2, C3 c3, C4 c4) {
  return CombinedConstraint<CombinedConstraint<C1,C2>,CombinedConstraint<C3,C4>>(make_combined_constraint(c1,c2),make_combined_constraint(c3,c4));
}


int main(int argc, char* argv[]) {
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TETS_FILE\n";
    exit(1);
  }

  GraphType graph;
  std::vector<typename GraphType::node_type> nodes;

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
  while (CS207::getline_parsed(tets_file, t)) {
    for (unsigned i = 1; i < t.size(); ++i) {
      graph.add_edge(nodes[t[0]], nodes[t[1]]);
      graph.add_edge(nodes[t[0]], nodes[t[2]]);
//#if 0
      // Diagonal edges: include as of HW2 #2
      graph.add_edge(nodes[t[0]], nodes[t[3]]);
      graph.add_edge(nodes[t[1]], nodes[t[2]]);
//#endif
      graph.add_edge(nodes[t[1]], nodes[t[3]]);
      graph.add_edge(nodes[t[2]], nodes[t[3]]);
    }
  }

  // Set initial conditions for your nodes, if necessary.
  for(GraphType::node_iterator ni=graph.node_begin(); ni!=graph.node_end(); ++ni) {
    (*ni).value().m = 1.0 / graph.num_nodes();
    (*ni).value().c = 1.0 / graph.num_nodes();
  }
  for(GraphType::edge_iterator ei=graph.edge_begin(); ei!=graph.edge_end(); ++ei) {
    (*ei).value().k = 100.;
    (*ei).value().l = (*ei).length();
  }
  
  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  auto node_map = viewer.empty_node_map(graph);
  viewer.launch();

  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.0001;
  double t_start = 0;
  double t_end = 5.0;

  for (double t = t_start; t < t_end; t += dt) {
    //std::cout << "t = " << t << std::endl;
    //symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(),DampingForce()), make_combined_constraint(ConstantNodeConstraint(graph),PlaneConstraint(-0.75),SphereConstraint(Point(0.5,0.5,-0.5),0.15)));
    symp_euler_step(graph, t, dt, make_combined_force(GravityForce(), MassSpringForce(),DampingForce()), make_combined_constraint(ConstantNodeConstraint(graph),PlaneConstraint(-0.75),SphereConstraint(Point(0.5,0.5,-0.5),0.15),AvoidSelfCollisionConstraint()) );
    // Update viewer with nodes' new positions
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    //if (graph.size() < 100)
    //  CS207::sleep(0.001);
  }

  return 0;
}
