/**
 * @file mtl_test.cpp
 * Test script for interfacing with MTL4 and it's linear solvers.
 */

#include <boost/numeric/mtl/mtl.hpp>
#include <boost/numeric/itl/itl.hpp>

// Define a IdentityMatrix that interfaces with MTL
struct IdentityMatrix {
  /** Compute the product of a vector with this identity matrix
  */
  
  /** constructor **/
  IdentityMatrix(std::size_t s_in) : s(s_in) {}
  
  /** Helper function to perform multiplication . Allows for delayed
   * evaluation of results and various assignment operations such
   * as += and -=.
   * @pre @a size (v) == size (w) 
   */
  template <typename VectorIn , typename VectorOut , typename Assign >
  void mult ( const VectorIn& v, VectorOut& w, Assign ) const {
  
    assert(std::size_t(size(v)) == s);
    assert(size(v) == size(w));
    
    for (std::size_t k = 0; k < s; ++k)
      Assign::apply(w[k], v[k]); 

  }

  /** Matrix-vector multiplication forwards to MTL's lazy mat_cvec_multiplier operator */
  template <typename VectorIn >
  mtl::vector::mat_cvec_multiplier<IdentityMatrix, VectorIn> operator*( const VectorIn& v) const {
    return mtl::vector::mat_cvec_multiplier<IdentityMatrix, VectorIn>(*this, v);
  }

  // /** slow,quick * operator */
  // template <typename Vector>
  // Vector operator*( const Vector & x ) const {
  //   return x;
  // }
  
  std::size_t s; // size of matrix
  
  private :
    //Empty!
};


/** The number of elements in the matrix . */
inline std::size_t size( const IdentityMatrix & A ) { return A.s*A.s; }

/** The number of rows in the matrix . */
inline std::size_t num_rows( const IdentityMatrix & A ) { return A.s; }

/** The number of columns in the matrix . */
inline std::size_t num_cols( const IdentityMatrix & A ) { return A.s; }


/** Traits that MTL uses to determine properties of our IdentityMatrix . */
namespace mtl { namespace ashape {

/** Define IdentityMatrix to be a non - scalar type . */
template <>
struct ashape_aux < IdentityMatrix > {
  typedef nonscal type ;
};

} // end namespace ashape

/** IdentityMatrix implements the Collection concept
* with value_type and size_type */
template <>
struct Collection <IdentityMatrix> {
  typedef double value_type ;
  typedef unsigned size_type ;
};

} // end namespace mtl



int main()
{
  // Construct an IdentityMatrix and "solve" Ix = b
  // using MTL's conjugate gradient solver

  using namespace std;
  typedef mtl::dense_vector<double> vec_type;
  typedef IdentityMatrix matrix_type;
  
  const std::size_t N = 20;
  
  cout << endl << "== TEST 1 ==" << endl;
  
  matrix_type A(N);
  
  // Create an ILU(0) preconditioner
  itl::pc::identity<matrix_type>        P(A);
    
  // Set b such that x == 1 is solution; start with x == 0
  vec_type x(N, 1.0), b(N);
  b = A*x;
  x = 0;
    
  // Termination criterion: r < 1e-11 * b or 100 iterations
  itl::cyclic_iteration<double>  iter(b, 100, 1.e-11, 0.0, 5);
    
  // Solve Ax == b with left preconditioner P
  itl::cg(A, x, b, P, iter);

  cout << "b is " << b << endl;
  cout << "x is " << x << endl;  
  
  
  cout << endl << "== TEST 2 ==" << endl;
  
  for(std::size_t i = 0; i < N; ++i)
    x[i] = i;
  b = A*x;
  x = 0;
  itl::cyclic_iteration<double>  iter2(b, 100, 1.e-11, 0.0, 5);
  itl::cg(A, x, b, P, iter2);
  
  cout << "b is " << b << endl;
  cout << "x is " << x << endl;  
  
  
  return 0;
}
