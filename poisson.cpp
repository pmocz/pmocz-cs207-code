/**
 * @file poisson.cpp
 * Test script for treating the Graph as a MTL Matrix
 * and solving a Poisson equation.
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles.
 * Second file: Edges (one per line) defined by 2 indices into the point list
 *              of the first file.
 *
 * Launches an SDLViewer to visualize the solution.
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include <math.h>       // cos, fabs

#include "Graph.hpp"

#include <fstream>

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// data held by a node
struct node_data {
  bool is_bndry;
};

// data held by an edge
struct edge_data {
};

typedef Graph<node_data,edge_data> GraphType;
typedef GraphType::Node Node;
typedef GraphType::Edge Edge;

double poisson_f(Point x) { return 5.0*cos(norm_1(x)); }

double poisson_g(Point x) {
  if( norm_inf(x) == 1 )
    return 0.0;
  if( norm_inf(x-Point(0.6,0.6,0.0))< 0.2)
    return -0.2;
  if( norm_inf(x-Point(0.6,-0.6,0.0))< 0.2)
    return -0.2;
  if( norm_inf(x-Point(-0.6,0.6,0.0))< 0.2)
    return -0.2;
  if( norm_inf(x-Point(-0.6,-0.6,0.0))< 0.2)
    return -0.2;    
  return 1.0;
}


// Define a GraphSymmetricMatrix that maps
// your Graph concept to MTL's Matrix concept. This shouldn't need to copy or
// modify Graph at all!
struct GraphSymmetricMatrix {
  /** Compute the product of a vector with the Poisson matrix represented by our graph
  */
  
  /** constructor **/
  GraphSymmetricMatrix(GraphType& g_in) : g(g_in), s(g.num_nodes()) {}
  
  /** Helper function to perform multiplication . Allows for delayed
   * evaluation of results and various assignment operations such
   * as += and -=.
   * @pre @a size (v) == size (w) 
   */
  template <typename VectorIn , typename VectorOut , typename Assign >
  void mult ( const VectorIn& v, VectorOut& w, Assign ) const {
  
    assert(std::size_t(size(v)) == s);
    assert(size(v) == size(w));
    
    unsigned k;
    double wk;
    for(GraphType::node_iterator ni=g.node_begin(); ni!=g.node_end(); ++ni) {
      k = (*ni).index(); 
      if( (*ni).value().is_bndry ) {  // boundary
        Assign::apply(w[k], v[k]);          
      }
      else {                          // interior
        wk = -v[k]*( (double) ((*ni).degree()) );
        for(GraphType::incident_iterator ii = (*ni).edge_begin(); ii!=(*ni).edge_end(); ++ii)
          if( !((*ii).node2().value().is_bndry) )
            wk += v[(*ii).node2().index()];
        Assign::apply(w[k], wk);
      }
    }
  }

  /** Matrix-vector multiplication forwards to MTL's lazy mat_cvec_multiplier operator */
  template <typename VectorIn >
  mtl::vector::mat_cvec_multiplier<GraphSymmetricMatrix, VectorIn> operator*( const VectorIn& v) const {
    return mtl::vector::mat_cvec_multiplier<GraphSymmetricMatrix, VectorIn>(*this, v);
  }
 
  GraphType& g;
  std::size_t s; // size of matrix
  
  private :
    //Empty!
};


/** The number of elements in the matrix . */
inline std::size_t size( const GraphSymmetricMatrix & A ) { return A.s*A.s; }

/** The number of rows in the matrix . */
inline std::size_t num_rows( const GraphSymmetricMatrix & A ) { return A.s; }

/** The number of columns in the matrix . */
inline std::size_t num_cols( const GraphSymmetricMatrix & A ) { return A.s; }


/** Traits that MTL uses to determine properties of our GraphSymmetricMatrix . */
namespace mtl { namespace ashape {

/** Define GraphSymmetricMatrix to be a non - scalar type . */
template <>
struct ashape_aux < GraphSymmetricMatrix > {
  typedef nonscal type ;
};

} // end namespace ashape

/** GraphSymmetricMatrix implements the Collection concept
* with value_type and size_type */
template <>
struct Collection <GraphSymmetricMatrix> {
  typedef double value_type ;
  typedef unsigned size_type ;
};

} // end namespace mtl




/** A color functor that returns colors for values it receives
 * @param[in] norm normalization value
 */
 template<class Vector>
struct HeatMap {
  template <typename NODE>
  CS207::Color operator()(const NODE& node) {
    return color_.make_heat((  max_ - max(min_,min(max_,(float) u_[node.index()] )) )/(max_-min_));
  }
  HeatMap(const float min, const float max, Vector& u)  : min_(min), max_(max), u_(u) {};
  private:
    float min_;
    float max_;
    CS207::Color color_;
    Vector& u_;
};


/** A position functor that sets z = the solution of the PDE
*/
template<class Vector>
struct MyPosition {
  MyPosition(Vector& u)  : u_(u) {};

  template <typename NODE>
  Point operator()(const NODE& node) {
    return node.position() + Point(0,0,u_[node.index()]);
  }
  
  private:
    Vector& u_;
};


/** visual iterator -- extends itl::cyclic_iterator to draw solution with viewer
*/
namespace itl {

  template <class Real, class Viewer, class Graph, class Vector, class OStream = std::ostream>
  class visual_iteration : public cyclic_iteration<Real> 
  {
      typedef cyclic_iteration<Real,OStream> super;
      typedef visual_iteration self;

    public:
  
      // constructors
      visual_iteration(Viewer& viewer_, Graph& graph_, Vector& current_sol_, const Vector& r0, int max_iter_, Real tol_, Real atol_ = Real(0), int cycle_ = 100,
                       OStream& out = std::cout)
        : super(r0, max_iter_, tol_, atol_, cycle_, out), viewer(viewer_), graph(graph_), current_sol(current_sol_)
      {}
      
      bool finished() { return super::finished(); }

      template <typename T>
      bool finished(const T& r) 
      {
          bool ret = super::finished(r);
          // Construct a Color functor and view with the SDLViewer
          viewer.clear();
          auto node_map = viewer.empty_node_map(graph);
          viewer.add_nodes(graph.node_begin(), graph.node_end(), HeatMap<Vector>(-0.2,1.0,current_sol), MyPosition<Vector>(current_sol), node_map);
          viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);
          viewer.set_label(this->i);
          viewer.center_view();
          //CS207::sleep(0.02);
          return ret;
      }
      
    protected:
      Viewer& viewer;
      Graph& graph;
      Vector& current_sol;
  };

} // namespace itl





int main(int argc, char** argv)
{
  // Check arguments
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NODES_FILE EDGES_FILE\n";
    exit(1);
  }
  
  typedef mtl::dense_vector<double> vec_type;
  
  GraphType graph;
  std::vector<typename GraphType::node_type> nodes;

  
  // Construct the Graph
  
  // Create a nodes_file from the first input argument
  std::ifstream nodes_file(argv[1]);
  // Interpret each line of the nodes_file as a 3D Point and add to the Graph
  Point p;
  while (CS207::getline_parsed(nodes_file, p))
    nodes.push_back(graph.add_node(p));

  // Create a tets_file from the second input argument
  std::ifstream tets_file(argv[2]);
  // Interpret each line of the tets_file as four ints which refer to nodes
  std::array<int,2> t;
  while (CS207::getline_parsed(tets_file, t))
    for (unsigned i = 1; i < t.size(); ++i)
      for (unsigned j = 0; j < i; ++j)
        graph.add_edge(nodes[t[i]], nodes[t[j]]);
  
  // find and label boundary nodes
  for(GraphType::node_iterator ni=graph.node_begin(); ni!=graph.node_end(); ++ni) {
    bool on_bndry = false;
    if( norm_inf((*ni).position()) == 1 || 
        norm_inf((*ni).position()-Point( 0.6, 0.6,0.0))< 0.2 ||
        norm_inf((*ni).position()-Point( 0.6,-0.6,0.0))< 0.2 ||
        norm_inf((*ni).position()-Point(-0.6, 0.6,0.0))< 0.2 ||
        norm_inf((*ni).position()-Point(-0.6,-0.6,0.0))< 0.2 ||
        (fabs((*ni).position().x)<= 0.6 && fabs((*ni).position().y)<= 0.2) ) {
      on_bndry = true;
    }   
    (*ni).value().is_bndry = on_bndry;
  }
  
  // Construct the GraphSymmetricMatrix
  GraphSymmetricMatrix A(graph);
  unsigned N = graph.num_nodes();
  GraphType::edge_iterator ei=graph.edge_begin();
  double h = (*ei).length();
  std::cout << "N nodes: " << graph.num_nodes() << std::endl;
  std::cout << "N edges: " << graph.num_edges() << std::endl;
  std::cout << "h: " << h << std::endl;
  
  // Define b
  vec_type b(N);
  for(GraphType::node_iterator ni=graph.node_begin(); ni!=graph.node_end(); ++ni) {
    if( (*ni).value().is_bndry ) {
      b[(*ni).index()] = poisson_g((*ni).position());
    }
    else {
      b[(*ni).index()] = h*h*poisson_f((*ni).position());
      for(GraphType::incident_iterator ii = (*ni).edge_begin(); ii!=(*ni).edge_end(); ++ii) {
        if( (*ii).node2().value().is_bndry )
          b[(*ni).index()] -= poisson_g((*ii).node2().position());
      }
    }
  }

  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  // Solve Au = b using MTL
    
  // Set initial guess u = 0
  vec_type x(N);
  x = 0.0;
   
  // Termination criterion: r < 1e-10 * b or 1000 iterations, print every 1 cycle
  //itl::cyclic_iteration<double>  iter(b, 1000, 1.e-10, 0.0, 50);
  itl::visual_iteration<double,CS207::SDLViewer,GraphType,vec_type>  iter(viewer, graph, x, b, 1000, 1.e-10, 0.0, 1);
    
  // Solve Ax == b with left preconditioner P
  // Create an ILU(0) preconditioner
  //itl::pc::identity<GraphSymmetricMatrix>  P(A);
  //itl::cg(A, x, b, P, iter);
  
  // Solve Ax == b
  itl::cg(A, x, b, iter);
  
  //auto node_map = viewer.empty_node_map(graph);
  //viewer.add_nodes(graph.node_begin(), graph.node_end(), HeatMap<vec_type>(-0.2,1.0,x), MyPosition<vec_type>(x), node_map);
  //viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  return 0;
}
