
#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "Mesh.hpp"
#include <fstream>

/** define Mesh types */
typedef Mesh<int,int,int> MeshType;
typedef MeshType::Node Node;

int main(int argc, char* argv[])
{
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
    
  //for(auto ni = mesh.node_begin(); ni != mesh.node_end(); ++ni ) {
  //  std::cout << (*ni).position() << std::endl;
  //} 
  
  for(auto ei = mesh.edge_begin(); ei != mesh.edge_end(); ++ei ) {
    std::cout << (*ei).node1().position() << std::endl;
  } 
  
  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as three ints which refer to nodes of triangle
  std::array<int,3> t;
  while (CS207::getline_parsed(tets_file, t))
    mesh.add_triangle(nodes[t[0]], nodes[t[1]], nodes[t[2]]);
 
  for(auto ti = mesh.triangle_begin(); ti != mesh.triangle_end(); ++ti ) {
    std::cout << (*ti).center_of_mass() << " area " << (*ti).norm_vector(1) <<  std::endl;
  } 
  
 std::cout << "testing incident iterator" <<  std::endl; 
  auto n = *(mesh.node_begin());
  std::cout << n.degree() <<  std::endl;
  for(auto ii = n.triangle_begin(); ii != n.triangle_end(); ++ii ) {
    std::cout << (*ii).index() <<  std::endl;
  }   
 
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch(); 
  
  auto node_map = viewer.empty_node_map(mesh); /// XXX change
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(), node_map);
  viewer.add_edges(mesh.edge_begin(), mesh.edge_end(), node_map);

  viewer.center_view();
 

 return 0;

}