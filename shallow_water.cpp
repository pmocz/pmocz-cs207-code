/**
 * @file shallow_water.cpp
 * Implementation of a shallow water system using Mesh
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "Point.hpp"
#include <fstream>
#include <math.h>       /* exp */

#include "Mesh.hpp"

// Standard gravity (average gravity at Earth's surface) in meters/sec^2
static constexpr double grav = 9.80665;

/** Water column characteristics */
struct QVar {
  double h;	  // Height of fluid
  double hu;	// Height times average x velocity of column
  double hv;	// Height times average y velocity of column

  /** Default constructor.
   *
   * A default of (0,0,0) *///A default water column is 1 unit high with no velocity. */
  QVar()
    : h(0), hu(0), hv(0) {
  }
  /** Construct the given water column. */
  QVar(double h_, double hu_, double hv_)
    : h(h_), hu(hu_), hv(hv_) {
  }
  // More operators
  QVar& operator+=(QVar q2) {
    h += q2.h;
    hu += q2.hu;
    hv += q2.hv;
    return *this;
  }
  QVar& operator-=(QVar q2) {
    h -= q2.h;
    hu -= q2.hu;
    hv -= q2.hv;
    return *this;
  }
  QVar& operator*=(double a) {
    h *= a;
    hu *= a;
    hv *= a;
    return *this;
  }
};

// operators
QVar operator*(QVar q, double a) {
  return q *= a;
}
QVar operator*(double a, QVar q) {
  return q *= a;
}
QVar operator+(QVar q1, QVar q2) {
  return q1 += q2;
}
QVar operator/(QVar q, double a) {
  return q *= (1.0/a);
}

/** Initial Condition Structures */
std::string ic_options_string = "{DamBreak,Pebble,FastWave}";
double H(const double& x) {
  if (x < 0)
    return 1;
  else 
    return 0;
}
struct DamBreak {
  QVar operator()(const Point& p) {
     return QVar(1.0+0.75*H(p.x),0,0);
  }
};
struct Pebble {
  QVar operator()(const Point& p) {
    return QVar(1.0-0.75*exp(-80*((p.x-0.75)*(p.x-0.75)+p.y*p.y)),0,0);
  }
};
struct FastWave {
  QVar operator()(const Point& p) {
     return QVar(1.0+0.75*H((p.x-0.75)*(p.x-0.75)+p.y*p.y-0.15*0.15),0,0);
  }
};

struct NodeData {
  QVar q;
};

struct EdgeData {
};

struct TriData {
  QVar q;
  QVar flux_sum;
};

/** Function object for calculating shallow-water flux.
 *          |e
 *   T_k    |---> n = (nx,ny)   T_m
 *   QBar_k |                   QBar_m
 *          |
 * @param[in] nx, ny Defines the 2D outward normal vector n = (@a nx, @a ny)
 *            from triangle T_k to triangle T_m. The length of n is equal to the
 *            the length of the edge, |n| = |e|.
 * @param[in] dt The time step taken by the simulation. Used to compute the
 *               Lax-Wendroff dissipation term.
 * @param[in] qk The values of the conserved variables on the left of the edge.
 * @param[in] qm The values of the conserved variables on the right of the edge.
 * @return The flux of the conserved values across the edge e
 */
struct EdgeFluxCalculator {
  QVar operator()(double nx, double ny, double dt,
                  const QVar& qk, const QVar& qm) {
    double e_length = sqrt(nx*nx + ny*ny);
    nx /= e_length;
    ny /= e_length;

    // The velocities normal to the edge
    double wm = (qm.hu*nx + qm.hv*ny) / qm.h;
    double wk = (qk.hu*nx + qk.hv*ny) / qk.h;

    // Lax-Wendroff local dissipation coefficient
    double vm = sqrt(grav*qm.h) + sqrt(qm.hu*qm.hu + qm.hv*qm.hv) / qm.h;
    double vk = sqrt(grav*qk.h) + sqrt(qk.hu*qk.hu + qk.hv*qk.hv) / qk.h;
    double a  = dt * std::max(vm*vm, vk*vk);

    // Helper values
    double scale = 0.5 * e_length;
    double gh2   = 0.5 * grav * (qm.h*qm.h + qk.h*qk.h);

    // Simple flux with dissipation for stability
    return QVar(scale * (wm*qm.h  + wk*qk.h)           - a * (qm.h  - qk.h),
                scale * (wm*qm.hu + wk*qk.hu + gh2*nx) - a * (qm.hu - qk.hu),
                scale * (wm*qm.hv + wk*qk.hv + gh2*ny) - a * (qm.hv - qk.hv));
  }
};

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) {
    return Point(n.position().x,n.position().y,n.value().q.h);
  }
};

// Mesh Type!
typedef Mesh<NodeData,EdgeData,TriData> MeshType;


/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time.
 */
template <typename MESH, typename FLUX>
double hyperbolic_step(MESH& m, FLUX& flux_function, double t, double dt) {
  // Step the finite volume model in time by dt.
  // Implements Equation 7 from HW4B handout
  // loop over all faces to compute and store the fluxes
  QVar f;
  for (auto eit = m.edge_begin(); eit != m.edge_end(); ++eit) { 
    // calculate flux for boundary triangles
    if ((*eit).num_triangles() == 1) {
      f = flux_function( (*eit).norm_vector_12().x, (*eit).norm_vector_12().y, dt, (*eit).triangle1().value().q, QVar( (*eit).triangle1().value().q.h,0.0,0.0) );
    }
    else {  // calculate flux for triangles in the inner region
      f = flux_function( (*eit).norm_vector_12().x, (*eit).norm_vector_12().y, dt, (*eit).triangle1().value().q, (*eit).triangle2().value().q );
      (*eit).triangle2().value().flux_sum -= f;
   }
    (*eit).triangle1().value().flux_sum += f;
  }

  // loop over all the triangles to add fluxes to q
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit) {
    (*tit).value().q -= dt/(*tit).area()*(*tit).value().flux_sum;
    (*tit).value().flux_sum *= 0;  // clear flux_sum for next iteration
  }
  return t + dt;
}

/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
void post_process(MESH& m) {
  // HW4B: Post-processing step
  // Translate the triangle-averaged values to node-averaged values
  // Implements Equation 8 in HW4B handout
  // update the node values by averaging over the adjacent triangles
  for (auto nit = m.node_begin(); nit != m.node_end(); ++nit) {
    double totalarea = 0;
    QVar acc = QVar(0,0,0);
    for (auto tit = (*nit).triangle_begin(); tit != (*nit).triangle_end(); ++tit) {
      totalarea += (*tit).area();
      acc += (*tit).area()*(*tit).value().q;
    }
    (*nit).value().q = acc*(1.0/totalarea);
  }
}

/** print out lots of info about the mesh for debugging purposes*/
template <typename MESH>
void debug_print(MESH& m, double t) {
  for(auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit ) {
    std::cout << "Triangle " << (*tit).index() << " @" << t << std::endl;
    std::cout << "  Area " << (*tit).area() << std::endl;
    std::cout << "  Node positions (" << (*tit).node(1).position() << ") (" << (*tit).node(2).position() << ") (" << (*tit).node(3).position() << ") ("  << std::endl;
    std::cout << "  Triangle QVar [Q_bar, water column characteristics] h=" << (*tit).value().q.h << " hu=" << (*tit).value().q.hu << " hv=" << (*tit).value().q.hv  << std::endl;
    std::cout << "  Edge 1 (" << (*tit).node(2).position() << ") (" << (*tit).node(3).position() << ")" << std::endl;
    std::cout << "    Normal (" << (*tit).norm_vector(1) << ")" << std::endl;
  } 
}


int main(int argc, char* argv[])
{
  std::cout << "\n== Shallow Water Solver by Xinyi Guo and Philip Mocz ==\n\n\n";

  // Check arguments
  if (argc < 4) {
    std::cerr << "Usage: shallow_water NODES_FILE TRIS_FILE IC_STRING\n\n         IC_STRING options " << ic_options_string <<"\n\n";
    exit(1);
  }

  MeshType mesh;
  std::vector<typename MeshType::node_type> mesh_node;

  // Read all Points and add them to the Mesh
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (CS207::getline_parsed(nodes_file, p)) {
    mesh_node.push_back(mesh.add_node(p));
  }

  // Read all mesh triangles and add them to the Mesh
  std::ifstream tris_file(argv[2]);
  std::array<int,3> t;
  while (CS207::getline_parsed(tris_file, t)) {
    mesh.add_triangle(mesh_node[t[0]], mesh_node[t[1]], mesh_node[t[2]]);
  }

  // Print out the stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;

  // HW4B Initialization
  // Set the initial conditions
  // Perform any needed precomputation
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                   CS207::DefaultColor(), NodePosition(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);


  /* Set the initial conditions */ 
  // Set the initial values of the nodes and get the maximum height double
  double max_h = 0; 
  std::string ic_string = argv[3];
  if (ic_string.compare("DamBreak")==0) {
    std::cout << "IC = DamBreak\n";
  }
  else if (ic_string.compare("Pebble")==0) {
    std::cout << "IC = Pebble\n";
  }
  else if (ic_string.compare("FastWave")==0) {
    std::cout << "IC = FastWave\n";
  }
  else {
    std::cerr << "IC_STRING error. \n";
  }
  auto init_cond_DamBreak = DamBreak();
  auto init_cond_Pebble = Pebble();
  auto init_cond_FastWave = FastWave();
  for (auto nit = mesh.node_begin(); nit != mesh.node_end(); ++nit) { 
    auto n = *nit; 
    if (ic_string.compare("DamBreak")==0) {
      n.value().q = init_cond_DamBreak(n.position()); 
    }
    else if (ic_string.compare("Pebble")==0) {
      n.value().q = init_cond_Pebble(n.position()); 
    }
    else if (ic_string.compare("FastWave")==0) {
      n.value().q = init_cond_FastWave(n.position()); 
    }
    max_h = std::max(max_h, n.value().q.h); 
  } 

  // Set the initial values of the triangles to the average of their nodes
  for (auto tit = mesh.triangle_begin(); tit != mesh.triangle_end(); ++tit) {
    auto tri = *tit; 
    tri.value().q = (tri.node(1).value().q + 
                     tri.node(2).value().q + 
                     tri.node(3).value().q) / 3.0; 
  }
  
  // HW4B: Timestep
  // CFL stability condition requires dt <= dx / max|velocity|
  // For the shallow water equations with u = v = 0 initial conditions
  //   we can compute the minimum edge length and maximum original water height
  //   to set the time-step
  // Compute the minimum edge length and maximum water height for computing dt
  double min_edge_length =  (*mesh.edge_begin()).length();
  for (auto eit = mesh.edge_begin(); eit != mesh.edge_end(); ++eit)
    if( (*eit).length() < min_edge_length ) 
      min_edge_length = (*eit).length();
  double max_height = max_h;
  double dt = 0.25 * min_edge_length / (sqrt(grav * max_height));
  std::cout << "dt = " << dt << std::endl;
  double t_start = 0;
  double t_end = 10;

  // Preconstruct a Flux functor
  EdgeFluxCalculator f;
  
  
  // Begin the time stepping
  for (double t = t_start; t < t_end; t += dt) {
    // debugging
    //debug_print(mesh,t);
    
    // Step forward in time with forward Euler
    hyperbolic_step(mesh, f, t, dt);

    // Update node values with triangle-averaged values
    post_process(mesh);

    // Update the viewer with new node positions
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                     CS207::DefaultColor(), NodePosition(), node_map);
    viewer.center_view();
    viewer.set_label(t);

    // These lines slow down the animation for small meshes.
    // Feel free to remove them or tweak the constants.
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.05);
  }

  return 0;
}
