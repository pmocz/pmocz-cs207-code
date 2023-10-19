/**
 * @file hw4b.cpp
 * Test script for solving the Shallow Water equations with our Mesh class.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Triangles (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "Mesh.hpp"
#include <fstream>

/** data structure for conserved variables of shallow water equations */
struct q_type{
  double h;
  double hu;
  double hv;
  q_type operator+(q_type q) //XXX
  q_type operator-(q_type q) //XXX
  /** constructor */
  q_type() : h(0), hu(0), hv(0) {}
};

/** data held by node */
struct node_data{
  q_type q;
};

/** data held by triangle */
struct triangle_data{
  q_type q;
  q_type flux_sum;
};

/** define Mesh types */
typedef Mesh<node_data,triangle_data> MeshType;
typedef MeshType::Node Node;
typedef MeshType::Face Face;
typedef MeshType::Triangle Triangle;

/** Lax-Friedrich's flux function */
q_type flux_function(q_type qbar_k, q_type qbar_m, Point n_km) {
  // to be supplied by Cris
  // XXX
} 


template<typename MESH>
void shallow_water_step(MESH& m, double t, double dt){
  // loop over all faces to compute and store the fluxes
  q_type f;
  for (auto fit = m.face_begin(); fit != m.face_end(); ++fit) {
    f = flux_function( (*fit).triangle1().value().q, (*fit).triangle2().value().q, (*fit).norm_vector_12() );
    (*fit).triangle1().value().flux_sum += f;  //XXX check signs
    (*fit).triangle2().value().flux_sum -= f;
  }

  // loop over all the triangles to add fluxes to q
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit) {
    (*tit).value().q -= dt/(*tit).area()*(*tit).value().flux_sum;
    (*tit).value().flux_sum = 0.;  // clear flux_sum for next iteration
  }

  // update the node values by averaging over the adjacent triangles
  for (auto nit = m.nodes_begin(); nit != m.nodes_end(); ++nit) {
    double totalarea = 0;
    q_type acc;
    for (auto tit = (*nit).triangle_begin(); tit !=  (*nit).triangle_end(); ++tit) {
      totalarea += (*tit).area();
      acc += (*tit).area()*(*tit).value().q;
    }
    (*nit).value() = acc/totalrea;
  }

}


void main{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE TRIANGLES_FILE\n";
    exit(1);
  }


  MeshType mesh;
  std::vector<typename MeshType::node_type> nodes;
 
  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Mesh
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(mesh.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as three ints which refer to nodes of triangle
  std::array<int,3> t;
  while (CS207::getline_parsed(tets_file, t))
    mesh.add_triangle(nodes[t[0]], nodes[t[1]], nodes[t[2]]);
 
  // Set ICs
  double t_start = 0;
  double t_end = 1.0;//XXX
  double dt = 0.1;//XXX
  // XXX
 
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch(); 
  
  auto node_map = viewer.empty_node_map(graph); /// XXX change
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
  viewer.add_edges(mesh.face_begin(), mesh.face_begin(), node_map);

  viewer.center_view();
 
  // Solve the shallow water equations and output to viewer
  for (double t = t_start; t < t_end; t += dt){
    shallow_water_step(mesh, t, dt);
    // Update viewer with nodes' new positions
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
  }
}