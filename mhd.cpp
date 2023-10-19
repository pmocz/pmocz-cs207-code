/**
 * @file mhd.cpp
 * Implementation of a magnetohydrodynamic (MHD) system using Mesh
 *
 * @brief Reads in two files specified on the command line and an initial condition string.
 * First file: 3D point list (one per line) defined by three doubles
 * Second file: Triangles (one per line) defined by 3 indices into the point list
 */

#include "CS207/SDLViewer.hpp"
#include "CS207/Util.hpp"
#include "CS207/Color.hpp"
#include "Point.hpp"
#include <fstream>
#include <math.h>       // exp,sqrt,M_PI

#include "Mesh.hpp"
#include <array>

// ideal gas gamma
static double ideal_gas_gamma;

// Number of fluid variables (6 for MHD)
static constexpr int n_var = 6;

/** Fluid element characteristics */
struct QVar {
  std::array<double,n_var> q;   // {{rho, rho_vx, rho_vy, rho_e, Bx, By}}
  /** Default constructor. */
  QVar() : q() {};
  /** Construct a given state. */
  QVar(std::array<double,n_var> q_) : q(q_) {};
  // More operators
  QVar& operator+=(QVar qvar2) {
    for (int i=0; i<n_var;++i)
      q[i] += qvar2.q[i];
    return *this;
  }
  QVar& operator-=(QVar qvar2) {
    for (int i=0; i<n_var;++i)
      q[i] -= qvar2.q[i];
    return *this;
  }
  QVar& operator*=(double a) {
    for (int i=0; i<n_var;++i)
      q[i] *= a;
    return *this;
  }
};
// operators
QVar operator*(QVar qvar, double a)    { return qvar *= a; };
QVar operator*(double a, QVar qvar)    { return qvar *= a; };
QVar operator+(QVar qvar1, QVar qvar2) { return qvar1 += qvar2; };
QVar operator/(QVar qvar, double a)    { return qvar *= (1.0/a); };
QVar operator-(QVar qvar1, QVar qvar2) { return qvar1 += -1.0*qvar2; };


/** Mesh data */
struct NodeData {
  QVar q;
  double Ez;         // electric field (for CT)
};
struct EdgeData {
  double phiB;       // magnetic flux thru face -- outward from triangle1() (for CT) 
  double Ez;         // electric field (for CT)
};
struct TriData {
  QVar q;            // conservative fluid variables
  std::array<QVar,4> k; // runge kutta flux
  std::array<QVar,4> y; // runge kutta state
  std::array<std::pair<QVar,QVar>,4> y_grad; // runge kutta gradient
};
// Mesh Type!
typedef Mesh<NodeData,EdgeData,TriData> MeshType;


/** Initial Condition Structures */
std::string ic_options_string = "{Implosion,OrszagTang,BlastWave,HydroInteractive}";
double H(const double& x, const double& y, const double m, const double b) {
  if (y > m*x+b)  // (x + y > 0.15)
    return 1;
  else 
    return 0;
}
double J(const double& x, const double& y) {
  if (sqrt(x*x+y*y) < 0.1)
    return 1;
  else 
    return 0;
}
struct Implosion {
  QVar operator()(const Point& p, const double m, const double b) {
     return QVar({{0.125 + (1-0.125)*H(p.x,p.y,m,b), 0, 0, 0.14/(ideal_gas_gamma-1.0) + (1-0.14)/(ideal_gas_gamma-1.0)*H(p.x,p.y,m,b), 0, 0}});
  }
};
struct OrszagTang {
  QVar operator()(const Point& p) {
     return QVar({{pow(ideal_gas_gamma,2)/(4*M_PI),
                   pow(ideal_gas_gamma,2)/(4*M_PI)*(-sin(2*M_PI*p.y)),
                   pow(ideal_gas_gamma,2)/(4*M_PI)*( sin(2*M_PI*p.x)),
                   ideal_gas_gamma/(4*M_PI*(ideal_gas_gamma-1)) + 0.5*pow(ideal_gas_gamma,2)/(4*M_PI)*(pow(sin(2*M_PI*p.y),2)+pow(sin(2*M_PI*p.x),2)) + 0.5*(pow(sin(2*M_PI*p.y),2)+pow(sin(4*M_PI*p.x),2))/(4*M_PI),
                   -sin(2*M_PI*p.y)/sqrt(4*M_PI),
                   sin(4*M_PI*p.x)/sqrt(4*M_PI)}});    
  }
  double set_bflux(MeshType::edge_type e) {
    double AL, AR;
    Point eL, eR;
    eL = e.node1().position();
    eR = e.node2().position();
    if(cross(e.norm_vector_12(),e.node1().position()-e.triangle1().center_of_mass()).z<0) {
      eL = e.node2().position();
      eR = e.node1().position();
    }
    AL = cos(4*M_PI*eL.x)/(4*M_PI*sqrt(4*M_PI)) + cos(2*M_PI*eL.y)/(2*M_PI*sqrt(4*M_PI));
    AR = cos(4*M_PI*eR.x)/(4*M_PI*sqrt(4*M_PI)) + cos(2*M_PI*eR.y)/(2*M_PI*sqrt(4*M_PI));
    return AL - AR;
  }
};
struct BlastWave {
  QVar operator()(const Point& p) {
     return QVar({{1.0, 0.0, 0.0, (0.1+0.9*J(p.x,p.y))/(ideal_gas_gamma-1.0)+0.5, 1.0/sqrt(2.0), 1.0/sqrt(2.0)}});    
  }
  double set_bflux(MeshType::edge_type e) {
    double AL, AR;
    Point eL, eR;
    eL = e.node1().position();
    eR = e.node2().position();
    if(cross(e.norm_vector_12(),e.node1().position()-e.triangle1().center_of_mass()).z<0) {
      eL = e.node2().position();
      eR = e.node1().position();
    }
    AL = -1.0/sqrt(2)*eL.x + 1.0/sqrt(2)*eL.y;
    AR = -1.0/sqrt(2)*eR.x + 1.0/sqrt(2)*eR.y;
    return AL - AR;
  }
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
  QVar operator()(double nx, double ny,
                  QVar qk, QVar qm,
                  const std::pair<QVar,QVar>& grad_qk, const std::pair<QVar,QVar>& grad_qm,
                  Point dkf, Point dmf
                  ) {
    double e_length = sqrt(nx*nx + ny*ny);
    nx /= e_length;
    ny /= e_length;
    
    // prediction step
    qk += dkf.x*grad_qk.first + dkf.y*grad_qk.second;
    qm += dmf.x*grad_qm.first + dmf.y*grad_qm.second;
    
    // average states
    double q0_star = 0.5*(qk.q[0]+qm.q[0]); 
    double q1_star = 0.5*(qk.q[1]+qm.q[1]); 
    double q2_star = 0.5*(qk.q[2]+qm.q[2]); 
    double q3_star = 0.5*(qk.q[3]+qm.q[3]); 
    double q4_star = 0.5*(qk.q[4]+qm.q[4]); 
    double q5_star = 0.5*(qk.q[5]+qm.q[5]); 
    
    // impose continuous normal B-fields at face
    double qkx_tmp = (q4_star*nx+q5_star*ny)*nx + (qk.q[4]*(-ny)+qk.q[5]*nx)*(-ny);
    double qky_tmp = (q4_star*nx+q5_star*ny)*ny + (qk.q[4]*(-ny)+qk.q[5]*nx)*nx;
    double qmx_tmp = (q4_star*nx+q5_star*ny)*nx + (qm.q[4]*(-ny)+qm.q[5]*nx)*(-ny);
    double qmy_tmp = (q4_star*nx+q5_star*ny)*ny + (qm.q[4]*(-ny)+qm.q[5]*nx)*nx;
    qk.q[4] = qkx_tmp;
    qk.q[5] = qky_tmp;
    qm.q[4] = qmx_tmp;
    qm.q[5] = qmy_tmp;
    
    // calculate total pressures (gas + magnetic)
    double p_k = (ideal_gas_gamma-1)*(qk.q[3]-0.5*(qk.q[1]*qk.q[1]+qk.q[2]*qk.q[2])/qk.q[0]-0.5*(qk.q[4]*qk.q[4]+qk.q[5]*qk.q[5])) + 0.5*(qk.q[4]*qk.q[4]+qk.q[5]*qk.q[5]);
    double p_m = (ideal_gas_gamma-1)*(qm.q[3]-0.5*(qm.q[1]*qm.q[1]+qm.q[2]*qm.q[2])/qm.q[0]-0.5*(qm.q[4]*qm.q[4]+qm.q[5]*qm.q[5])) + 0.5*(qm.q[4]*qm.q[4]+qm.q[5]*qm.q[5]);
    double p_star = (ideal_gas_gamma-1)*(q3_star-0.5*(q1_star*q1_star+q2_star*q2_star)/q0_star-0.5*(q4_star*q4_star+q5_star*q5_star)) + 0.5*(q4_star*q4_star+q5_star*q5_star);
    
    // compute flux
    double f1_x = q1_star;
    double f2_x = q1_star*q1_star/q0_star+p_star-q4_star*q4_star;
    double f3_x = q1_star*q2_star/q0_star-q4_star*q5_star;
    double f4_x = (q3_star+p_star)*q1_star/q0_star-q4_star*(q4_star*q1_star+q5_star*q2_star)/q0_star;
    double f5_x = 0.0;
    double f6_x = q1_star/q0_star*q5_star - q2_star/q0_star*q4_star;
    
    double f1_y = q2_star;
    double f2_y = q2_star*q1_star/q0_star-q5_star*q4_star;
    double f3_y = q2_star*q2_star/q0_star+p_star-q5_star*q5_star;
    double f4_y = (q3_star+p_star)*q2_star/q0_star-q5_star*(q4_star*q1_star+q5_star*q2_star)/q0_star;   
    double f5_y = q2_star/q0_star*q4_star - q1_star/q0_star*q5_star;   
    double f6_y = 0.0;      
    
    // wave speeds
    double c0_k = sqrt(ideal_gas_gamma*(p_k-0.5*(qk.q[4]*qk.q[4]+qk.q[5]*qk.q[5]))/qk.q[0]);
    double c0_m = sqrt(ideal_gas_gamma*(p_m-0.5*(qm.q[4]*qm.q[4]+qm.q[5]*qm.q[5]))/qm.q[0]);
    double can_k = sqrt(pow(qk.q[4]*nx+qk.q[5]*ny,2)/qk.q[0]);
    double can_m = sqrt(pow(qm.q[4]*nx+qm.q[5]*ny,2)/qm.q[0]);
    double ca_k = sqrt((qk.q[4]*qk.q[4]+qk.q[5]*qk.q[5])/qk.q[0]);
    double ca_m = sqrt((qm.q[4]*qm.q[4]+qm.q[5]*qm.q[5])/qm.q[0]);
    double cf_k = sqrt( 0.5*(c0_k*c0_k+ca_k*ca_k) + 0.5*sqrt(pow(c0_k*c0_k+ca_k*ca_k,2)-4*c0_k*c0_k*can_k*can_k) );
    double cf_m = sqrt( 0.5*(c0_m*c0_m+ca_m*ca_m) + 0.5*sqrt(pow(c0_m*c0_m+ca_m*ca_m,2)-4*c0_m*c0_m*can_m*can_m) );
    double vn_k = qk.q[1]/qk.q[0]*nx + qk.q[2]/qk.q[0]*ny;
    double vn_m = qm.q[1]/qm.q[0]*nx + qm.q[2]/qm.q[0]*ny;
    double c = max(max(fabs(vn_k+cf_k),fabs(vn_k-cf_k)),max(fabs(vn_m+cf_m),fabs(vn_m-cf_m)));
    
    // return local Lax Friedrich's flux
    return e_length* QVar({{f1_x*nx+f1_y*ny - 0.5*c*(qm.q[0]-qk.q[0]),
                            f2_x*nx+f2_y*ny - 0.5*c*(qm.q[1]-qk.q[1]),
                            f3_x*nx+f3_y*ny - 0.5*c*(qm.q[2]-qk.q[2]),
                            f4_x*nx+f4_y*ny - 0.5*c*(qm.q[3]-qk.q[3]),
                            f5_x*nx+f5_y*ny - 0.5*c*(qm.q[4]-qk.q[4]),
                            f6_x*nx+f6_y*ny - 0.5*c*(qm.q[5]-qk.q[5])}});

  }
};

/** Node position function object for use in the SDLViewer. */
struct NodePosition {
  template <typename NODE>
  Point operator()(const NODE& n) {
    return Point(n.position().x,n.position().y,n.value().q.q[0]);
  }
};


/** Get gradients of y[i]. */
template <typename MESH>
void get_gradient(int i, MESH& m) {
  // clear gradients
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit) {
    (*tit).value().y_grad[i].first  *=0;
    (*tit).value().y_grad[i].second *=0;
  }
  
  for (auto eit = m.edge_begin(); eit != m.edge_end(); ++eit) { 
    // calculate gradient for boundary triangles
    if ((*eit).num_triangles() == 1) {
      (*eit).triangle1().value().y_grad[i].first += 
        0.5*( (*eit).triangle1().value().y[i] + QVar( {{(*eit).triangle1().value().y[i].q[0],0.0,0.0,(*eit).triangle1().value().y[i].q[3],(*eit).triangle1().value().y[i].q[4],(*eit).triangle1().value().y[i].q[5]}}) )
        *(*eit).norm_vector_12().x;
      (*eit).triangle1().value().y_grad[i].second += 
        0.5*( (*eit).triangle1().value().y[i] + QVar( {{(*eit).triangle1().value().y[i].q[0],0.0,0.0,(*eit).triangle1().value().y[i].q[3],(*eit).triangle1().value().y[i].q[4],(*eit).triangle1().value().y[i].q[5]}}) )
        *(*eit).norm_vector_12().y;
    }
    else {  // calculate gradient for triangles in the inner region
      (*eit).triangle1().value().y_grad[i].first += 
        0.5*( (*eit).triangle1().value().y[i] + (*eit).triangle2().value().y[i] )
        *(*eit).norm_vector_12().x;
      (*eit).triangle1().value().y_grad[i].second += 
        0.5*( (*eit).triangle1().value().y[i] + (*eit).triangle2().value().y[i] )
        *(*eit).norm_vector_12().y;
      (*eit).triangle2().value().y_grad[i].first -= 
        0.5*( (*eit).triangle1().value().y[i] + (*eit).triangle2().value().y[i] )
        *(*eit).norm_vector_12().x;
      (*eit).triangle2().value().y_grad[i].second -= 
        0.5*( (*eit).triangle1().value().y[i] + (*eit).triangle2().value().y[i] )
        *(*eit).norm_vector_12().y;
    }
  }

  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit) {
    (*tit).value().y_grad[i].first  *= 1.0/(*tit).area();  
    (*tit).value().y_grad[i].second *= 1.0/(*tit).area();  
  }
  
  // slope limiter
  QVar delta_y, nbr_diff;
  Point d;
  int j;
  for (auto eit = m.edge_begin(); eit != m.edge_end(); ++eit) { 
    if ((*eit).num_triangles() == 2) {
      // triangle 1
      d = (*eit).center_of_mass() - (*eit).triangle1().center_of_mass();
      nbr_diff = (*eit).triangle2().value().y[i] - (*eit).triangle1().value().y[i];
      delta_y = (*eit).triangle1().value().y_grad[i].first * d.x + (*eit).triangle1().value().y_grad[i].second * d.y;
      for (j=0; j<n_var; ++j) {
        if( (delta_y.q[j] > 0 && nbr_diff.q[j] > 0) || (delta_y.q[j] < 0 && nbr_diff.q[j] < 0) ) {
          (*eit).triangle1().value().y_grad[i].first.q[j] *= min( 1.0, nbr_diff.q[j]/delta_y.q[j] );
          (*eit).triangle1().value().y_grad[i].second.q[j] *= min( 1.0, nbr_diff.q[j]/delta_y.q[j] );
        }
        else {
          (*eit).triangle1().value().y_grad[i].first.q[j] *= 0;
          (*eit).triangle1().value().y_grad[i].second.q[j] *= 0;
        }
      }      
      // triangle 2
      d = (*eit).center_of_mass() - (*eit).triangle2().center_of_mass();
      nbr_diff = (*eit).triangle1().value().y[i] - (*eit).triangle2().value().y[i];
      delta_y = (*eit).triangle2().value().y_grad[i].first * d.x + (*eit).triangle2().value().y_grad[i].second * d.y;
      for (j=0; j<n_var; ++j) {
        if( (delta_y.q[j] > 0 && nbr_diff.q[j] > 0) || (delta_y.q[j] < 0 && nbr_diff.q[j] < 0) ) {
          (*eit).triangle2().value().y_grad[i].first.q[j] *= min( 1.0, nbr_diff.q[j]/delta_y.q[j] );
          (*eit).triangle2().value().y_grad[i].second.q[j] *= min( 1.0, nbr_diff.q[j]/delta_y.q[j] );
        }
        else {
          (*eit).triangle2().value().y_grad[i].first.q[j] *= 0;
          (*eit).triangle2().value().y_grad[i].second.q[j] *= 0;
        }
      }
    }
    else {
      for (j=0; j<n_var;++j) {
        (*eit).triangle1().value().y_grad[i].first.q[j] *= 0;
        (*eit).triangle1().value().y_grad[i].second.q[j] *= 0;
      }
    }
  }

}


/** Calculate volume average magnetic fields from the face fluxes
 */
template <typename MESH>
void correct_b_fld(MESH& m) {
  double bn1,bn2,sign1,sign2,det;
  Point n1,n2;
  for(auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit) {
      sign1 = 1;
      sign2 = 1;     
      if( (*tit) != (*tit).edge(1).triangle1() )
        sign1 = -1;
      if( (*tit) != (*tit).edge(2).triangle1() )
        sign2 = -1;
      
      bn1 = sign1 * (*tit).edge(1).value().phiB / (*tit).edge(1).length();
      bn2 = sign2 * (*tit).edge(2).value().phiB / (*tit).edge(2).length();
    
      n1 = sign1 * (*tit).edge(1).norm_vector_12()/(*tit).edge(1).length();
      n2 = sign2 * (*tit).edge(2).norm_vector_12()/(*tit).edge(2).length();
      
      det = n1.x*n2.y - n1.y*n2.x;   
        
      (*tit).value().q.q[4] = ( n2.y*bn1-n1.y*bn2)/det;
      (*tit).value().q.q[5] = (-n2.x*bn1+n1.x*bn2)/det;
  }
}


/** helper function for Runge-Kutta integrator
 *
 */
template <typename MESH, typename FLUX>
void runge_kutta_step_i(int i, MESH& m, FLUX& flux_function) {
  QVar f;
  for (auto eit = m.edge_begin(); eit != m.edge_end(); ++eit) { 
    // calculate flux for boundary triangles
    if ((*eit).num_triangles() == 1) {
      f = flux_function( (*eit).norm_vector_12().x, 
                         (*eit).norm_vector_12().y, 
                         (*eit).triangle1().value().y[i], 
                         QVar( {{(*eit).triangle1().value().y[i].q[0],0.0,0.0,(*eit).triangle1().value().y[i].q[3],(*eit).triangle1().value().y[i].q[4],(*eit).triangle1().value().y[i].q[5]}}),
                         (*eit).triangle1().value().y_grad[i],
                         std::pair<QVar,QVar>(),
                         (*eit).center_of_mass()-(*eit).triangle1().center_of_mass(),
                         Point(0.,0.,0.)
                         );
    }
    else {  // calculate flux for triangles in the inner region
      f = flux_function( (*eit).norm_vector_12().x, 
                         (*eit).norm_vector_12().y, 
                         (*eit).triangle1().value().y[i], 
                         (*eit).triangle2().value().y[i],
                         (*eit).triangle1().value().y_grad[i],
                         (*eit).triangle2().value().y_grad[i],
                         (*eit).center_of_mass()-(*eit).triangle1().center_of_mass(),
                         (*eit).center_of_mass()-(*eit).triangle2().center_of_mass()
                         );
      (*eit).triangle2().value().k[i] -= f;
   }
    (*eit).triangle1().value().k[i] += f;
    (*eit).value().Ez = -(f.q[5]*(*eit).norm_vector_12().x - f.q[4]*(*eit).norm_vector_12().y)/pow((*eit).length(),2);
  }
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit)
    (*tit).value().k[i] *= -1.0/(*tit).area();
}
 
 
/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time using Runge-Kutta 4th order.
 */
template <typename MESH, typename FLUX>
double runge_kutta4(MESH& m, FLUX& flux_function, double t, double dt) {
  
  // y0
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit)
    (*tit).value().y[0] = (*tit).value().q;
  get_gradient(0, m);
  // k0
  runge_kutta_step_i(0, m, flux_function);
  // y1
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit)
    (*tit).value().y[1] = (*tit).value().q + dt*(*tit).value().k[0]/2.;
  get_gradient(1, m);
  // k1
  runge_kutta_step_i(1, m, flux_function);
  // y2
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit)
    (*tit).value().y[2] = (*tit).value().q + dt*(*tit).value().k[1]/2.;
  get_gradient(2, m);
  // k2
  runge_kutta_step_i(2, m, flux_function);
  // y3
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit)
    (*tit).value().y[3] = (*tit).value().q + dt*(*tit).value().k[2];  
  get_gradient(3, m);
  // k3
  runge_kutta_step_i(3, m, flux_function);
    
  // loop over all the triangles to add fluxes to q
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit) {
    (*tit).value().q += (dt/6.) * ( (*tit).value().k[0] + 2*(*tit).value().k[1] + 2*(*tit).value().k[2] + (*tit).value().k[3] );
    for(int i = 0; i != 4; ++i) {
      (*tit).value().k[i] *= 0;
      (*tit).value().y[i] *= 0;
    }
  }
  return t + dt;
}

/** Integrate a hyperbolic conservation law defined over the mesh m
 * with flux functor f by dt in time using Runge-Kutta 2nd order.
 */
template <typename MESH, typename FLUX>
double runge_kutta2(MESH& m, FLUX& flux_function, double t, double dt) {
  // y0
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit)
    (*tit).value().y[0] = (*tit).value().q;
  get_gradient(0, m);
  // k0
  runge_kutta_step_i(0, m, flux_function);
  // y1
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit)
    (*tit).value().y[1] = (*tit).value().q + dt*(*tit).value().k[0]/2.;
  get_gradient(1, m);
  // k1
  runge_kutta_step_i(1, m, flux_function);
  
  // CT: calculate the node electric fields
  calculate_electric_fields(m);

  // CT: update the face magnetic fluxes
  update_phiB(m,dt);

  // loop over all the triangles to add fluxes to q
  for (auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit) {
    (*tit).value().q += dt * (*tit).value().k[1];
    for(int i = 0; i != 2; ++i) {
      (*tit).value().k[i] *= 0;
      (*tit).value().y[i] *= 0;
    }
  }
  correct_b_fld(m);
  return t + dt;
}


/** Calculate electric fields at triangle nodes from Riemann fluxes for the CT method. */
template <typename MESH>
void calculate_electric_fields(MESH& m) {
  for (auto nit = m.node_begin(); nit != m.node_end(); ++nit) {
    double norm = 0;
    double acc = 0;
    for (auto eit = (*nit).edge_begin(); eit != (*nit).edge_end(); ++eit) {
      norm += 1.0/(*eit).length();
      acc += 1.0/(*eit).length()*(*eit).value().Ez;
    }
    (*nit).value().Ez = acc/norm;
  }
}


/** Calculate electric fields at triangle nodes from Riemann fluxes for the CT method. */
template <typename MESH>
void update_phiB(MESH& m, double dt) {
  double EL, ER;
  for (auto eit = m.edge_begin(); eit != m.edge_end(); ++eit) {
    EL = (*eit).node1().value().Ez;
    ER = (*eit).node2().value().Ez;
    if(cross((*eit).norm_vector_12(),(*eit).node1().position()-(*eit).triangle1().center_of_mass()).z<0) {
      EL = (*eit).node2().value().Ez;
      ER = (*eit).node1().value().Ez;
    }
    (*eit).value().phiB -= dt * (EL - ER);
  }
}


/** Convert the triangle-averaged values to node-averaged values for viewing. */
template <typename MESH>
void post_process(MESH& m) {
  // Translate the triangle-averaged values to node-averaged values (HW4B eqn 8)
  // update the node values by averaging over the adjacent triangles
  for (auto nit = m.node_begin(); nit != m.node_end(); ++nit) {
    double totalarea = 0;
    QVar acc = QVar();
    for (auto tit = (*nit).triangle_begin(); tit != (*nit).triangle_end(); ++tit) {
      totalarea += (*tit).area();
      acc += (*tit).area()*(*tit).value().q;
    }
    (*nit).value().q = acc*(1.0/totalarea);
  }
}

/** Get cell sound speed */
template <typename TRIANGLE>
double get_sound_speed(TRIANGLE t) {
  return sqrt( ideal_gas_gamma * 
              (ideal_gas_gamma-1)*(t.value().q.q[3]-0.5*(pow(t.value().q.q[1],2)+pow(t.value().q.q[2],2))/t.value().q.q[0]-0.5*(pow(t.value().q.q[4],2)+pow(t.value().q.q[5],2))) / 
              t.value().q.q[0] );
}

/** Get cell fast wave speed */
template <typename TRIANGLE>
double get_fast_wave_speed(TRIANGLE t) {
  double c0 = get_sound_speed(t);
  double ca = sqrt((pow(t.value().q.q[4],2)+pow(t.value().q.q[5],2))/t.value().q.q[0]);
  return sqrt( c0*c0+ca*ca );
}

/** Get time step from CFL condition*/
template <typename MESH>
double get_time_step(double min_edge_length, MESH& mesh) {
  auto tit = mesh.triangle_begin();
  double max_vel = get_fast_wave_speed(*tit) + sqrt( pow((*tit).value().q.q[1],2) + pow((*tit).value().q.q[2],2) ) / (*tit).value().q.q[0];
  double cell_vel;
  for (tit = mesh.triangle_begin(); tit != mesh.triangle_end(); ++tit) {
    cell_vel = get_fast_wave_speed(*tit) + sqrt( pow((*tit).value().q.q[1],2) + pow((*tit).value().q.q[2],2) ) / (*tit).value().q.q[0];
    if( cell_vel > max_vel ) 
      max_vel = cell_vel;
  }
  return 0.4*min_edge_length/max_vel;
}


/** print out lots of info about the mesh for debugging purposes*/
template <typename MESH>
void debug_print(MESH& m, double t) {
  for(auto tit = m.triangle_begin(); tit != m.triangle_end(); ++tit ) {
    std::cout << "Triangle " << (*tit).index() << " @" << t << std::endl;
    std::cout << "  Area " << (*tit).area() << std::endl;
    std::cout << "  Node positions (" << (*tit).node(1).position() << ") (" << (*tit).node(2).position() << ") (" << (*tit).node(3).position() << ") ("  << std::endl;
  } 
}


/** A color functor that returns colors for values it receives
 * @param[in] norm normalization value
 */
struct HeatMap {
  template <typename NODE>
  CS207::Color operator()(const NODE& node) {
    return color_.make_heat( (max(cmin_,min(cmax_, (float) node.value().q.q[0] )) - cmin_) /(cmax_ - cmin_) );
  }
  HeatMap(const float cmin, const float cmax) : cmin_(cmin), cmax_(cmax) {};
  private:
    float cmin_;
    float cmax_;
    CS207::Color color_;
};


/** return the distance from a point to line defined by pair of points
 * http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
 */
double dist_pt2line(const Point& pt, const std::pair<Point,Point>& ln) {
  return  norm(cross(pt-ln.first,pt-ln.second))/norm(ln.second-ln.first);
}


/** Basic SDL listener, right click adds energy to triangle
*/
struct add_energy_listener : public CS207::SDL_Listener {
  void handle(SDL_Event e) {   // we are forced to implement this function
    switch (e.type) {
      case SDL_MOUSEBUTTONDOWN: {
        if (e.button.button == SDL_BUTTON_RIGHT ) {
          std::pair<Point,Point> unproj_pts;
          unproj_pts = viewer_.unProject(e.button.x,e.button.y);
          auto tit = mesh_.triangle_begin();
          double min_dist = dist_pt2line((*tit).center_of_mass(),unproj_pts);
          MeshType::triangle_type t;
          double d;
          for (auto tit = mesh_.triangle_begin(); tit != mesh_.triangle_end(); ++tit) { 
            d = dist_pt2line((*tit).center_of_mass(),unproj_pts);
            if( d < min_dist ) {
              min_dist = d;
              t = *tit;
            }
          }
          t.value().q.q[3] *= 2; // double the energy
        }
      } break;
    }
  }
  
  // constructor
  add_energy_listener(CS207::SDLViewer& viewer, MeshType& mesh) : viewer_(viewer), mesh_(mesh) {};
  private:
   CS207::SDLViewer& viewer_;
   MeshType& mesh_;
};


struct set_ic_listener : public CS207::SDL_Listener {
  void handle(SDL_Event e) {   // we are forced to implement this function
    Point p,p1,p2;
    double min_dist1,min_dist2,d;
    MeshType::triangle_type t;
    double m,b;
    m = -1;
    b = 0.15;
    if (counter_<2) {
      if(counter_==0 && !printed_str1_) {
        std::cout << "LEFT-CLICK IC POINT 1\n";
        printed_str1_ = true;
      }
      if(counter_==1 && !printed_str2_) {
        std::cout << "LEFT-CLICK IC POINT 2\n"; 
        printed_str2_ = true;
      }        
      switch (e.type) {
        case SDL_MOUSEBUTTONDOWN: {
          if (e.button.button == SDL_BUTTON_LEFT ) {
            if (counter_==1) {
              unproj_pts2_ = viewer_.unProject(e.button.x,e.button.y);
              auto tit = mesh_.triangle_begin();
              min_dist1 = dist_pt2line((*tit).center_of_mass(),unproj_pts1_);
              min_dist2 = dist_pt2line((*tit).center_of_mass(),unproj_pts2_);
              for (auto tit = mesh_.triangle_begin(); tit != mesh_.triangle_end(); ++tit) { 
                d = dist_pt2line((*tit).center_of_mass(),unproj_pts1_);
                if( d < min_dist1 ) {
                  min_dist1 = d;
                  p1 = (*tit).center_of_mass();
                }
                d = dist_pt2line((*tit).center_of_mass(),unproj_pts2_);
                if( d < min_dist2 ) {
                  min_dist2 = d;
                  p2 = (*tit).center_of_mass();
                }
              }
              if(p2.x!=p1.x) {
                m = (p2.y-p1.y)/(p2.x-p1.x);
                b = p1.y-m*p1.x;
              }
              // ask user for density contrast
              std::cout << "SET DENSITY CONTRAST (input double; e.g. 8):"; 
              double d_constrast, p_constrast, shear_vel;
              bool include_shear_vel;
              bool input_set = false;
              while(!input_set) {
                std::cin >> d_constrast;
                input_set = true;
              }
              std::cout << "SET PRESSURE CONTRAST (input double; e.g. 7):"; 
              input_set = false;
              while(!input_set) {
                std::cin >> p_constrast;
                input_set = true;
              }
              std::cout << "ADD SHEAR VELOCITY? (0 no, 1 yes):"; 
              input_set = false;
              while(!input_set) {
                std::cin >> include_shear_vel;
                shear_vel = 0;
                if(include_shear_vel)
                  shear_vel = 1.0;
                input_set = true;
              }
              
              for (auto tit = mesh_.triangle_begin(); tit != mesh_.triangle_end(); ++tit) {
                p = (*tit).center_of_mass();
                double rho = 1.0/d_constrast + (1-1.0/d_constrast)*H(p.x,p.y,m,b);
                double rho_vx = rho*0.5*shear_vel*(  (1-H(p.x,p.y,m,b))*(sqrt(1.0/(m*m+1))) + H(p.x,p.y,m,b)*(-sqrt(1.0/(m*m+1)))  );
                double rho_vy = rho*0.5*shear_vel*(  (1-H(p.x,p.y,m,b))*(sqrt(m*m/(m*m+1))) + H(p.x,p.y,m,b)*(-sqrt(m*m/(m*m+1)))  );
                double bx = 0;
                double by = 0;
                double rho_e = 1.0/p_constrast/(ideal_gas_gamma-1.0) + (1-1.0/p_constrast)/(ideal_gas_gamma-1.0)*H(p.x,p.y,m,b) + 0.5*(rho_vx*rho_vx+rho_vy*rho_vy)/rho + 0.5*(bx*bx+by*by);
                (*tit).value().q = QVar({{rho, rho_vx, rho_vy, rho_e, bx, by}});
              }
              ++counter_;
              run_flag_ = true;
            }
            
            if (counter_==0) {
              unproj_pts1_ = viewer_.unProject(e.button.x,e.button.y);
              ++counter_;
            }
          }
        } break;
      }
    }
  }
  
  // constructor
  set_ic_listener(CS207::SDLViewer& viewer, MeshType& mesh, bool& run_flag) : 
      viewer_(viewer), 
      mesh_(mesh), 
      counter_(0), 
      run_flag_(run_flag), 
      unproj_pts1_(), 
      unproj_pts2_(), 
      printed_str1_(false), 
      printed_str2_(false) {};
  private:
   CS207::SDLViewer& viewer_;
   MeshType& mesh_;
   unsigned counter_;
   bool& run_flag_;
   std::pair<Point,Point> unproj_pts1_,unproj_pts2_;
   bool printed_str1_,printed_str2_;
};


int main(int argc, char* argv[])
{
  std::cout << "\n== Euler/MHD Solver by Xinyi Guo and Philip Mocz ==\n\n\n";

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

  // Print out stats
  std::cout << mesh.num_nodes() << " "
            << mesh.num_edges() << " "
            << mesh.num_triangles() << std::endl;
      
      
  // Launch the SDLViewer
  CS207::SDLViewer viewer;
  viewer.launch();

  auto node_map = viewer.empty_node_map(mesh);
  viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                   CS207::DefaultColor(), CS207::DefaultPosition(), node_map);
  viewer.add_triangles(mesh.triangle_begin(), mesh.triangle_end(), node_map);
  viewer.xyplane_view();

  // add energy listener
  add_energy_listener* l = new add_energy_listener(viewer,mesh); 
  viewer.add_listener(l); 
            
            
  /** Set the initial conditions */ 
  double cmin = 0;
  double cmax = 1;
  auto init_cond_Implosion = Implosion();
  auto init_cond_OrszagTang = OrszagTang();
  auto init_cond_BlastWave = BlastWave();
  std::string ic_string = argv[3];
  if (ic_string.compare("Implosion")==0) {
    std::cout << "IC = Implosion\n";
    ideal_gas_gamma = 1.4;
    cmin = 0;
    cmax = 1.1;
    for (auto nit = mesh.node_begin(); nit != mesh.node_end(); ++nit)
      (*nit).value().q = init_cond_Implosion((*nit).position(),-1.0,0.15); 
    for (auto tit = mesh.triangle_begin(); tit != mesh.triangle_end(); ++tit)
      (*tit).value().q = init_cond_Implosion((*tit).center_of_mass(),-1.0,0.15);
  } 
  else  if (ic_string.compare("HydroInteractive")==0) {
    std::cout << "IC = HydroInteractive\n";
    ideal_gas_gamma = 1.4;
    cmin = 0;
    cmax = 1.1;
    // add ic listener
    bool run_flag = false;
    set_ic_listener* lic = new set_ic_listener(viewer,mesh,run_flag); 
    viewer.add_listener(lic); 
    while(run_flag == false) {
      CS207::sleep(0.05);
    }
  } 
  else if (ic_string.compare("OrszagTang")==0) {
    std::cout << "IC = OrszagTang\n";
    ideal_gas_gamma = 1.66666667;
    cmin = 0.1;
    cmax = 0.4;
    for (auto nit = mesh.node_begin(); nit != mesh.node_end(); ++nit)
      (*nit).value().q = init_cond_OrszagTang((*nit).position()); 
    for (auto tit = mesh.triangle_begin(); tit != mesh.triangle_end(); ++tit)
      (*tit).value().q = init_cond_OrszagTang((*tit).center_of_mass());
    for (auto eit = mesh.edge_begin(); eit != mesh.edge_end(); ++eit)
      (*eit).value().phiB = init_cond_OrszagTang.set_bflux(*eit);  
    correct_b_fld(mesh);
  }
  else if (ic_string.compare("BlastWave")==0) {
    std::cout << "IC = BlastWave\n";
    ideal_gas_gamma = 1.66666667;
    cmin = 0.5;
    cmax = 1.4;
    for (auto nit = mesh.node_begin(); nit != mesh.node_end(); ++nit)
      (*nit).value().q = init_cond_BlastWave((*nit).position());
    for (auto tit = mesh.triangle_begin(); tit != mesh.triangle_end(); ++tit)
      (*tit).value().q = init_cond_BlastWave((*tit).center_of_mass());
    for (auto eit = mesh.edge_begin(); eit != mesh.edge_end(); ++eit)
      (*eit).value().phiB = init_cond_BlastWave.set_bflux(*eit);  
    correct_b_fld(mesh);
  }
  else {
    std::cerr << "IC_STRING error. \n\n         IC_STRING options " << ic_options_string <<"\n\n";
  }              
  
  
  // Timestep -- CFL stability condition requires dt <= dx / max|velocity|
  double min_edge_length =  (*mesh.edge_begin()).length();
  for (auto eit = mesh.edge_begin(); eit != mesh.edge_end(); ++eit)
    if( (*eit).length() < min_edge_length ) 
      min_edge_length = (*eit).length();
  double dt = get_time_step(min_edge_length,mesh);
  std::cout << "dt = " << dt << std::endl;
  double t_start = 0;
  double t_end = 10; 

  // Preconstruct a Flux functor
  EdgeFluxCalculator flux_calculator;
  
  
  // Begin the time stepping
  for (double t = t_start; t < t_end; t += dt) {
    // debugging
    //debug_print(mesh,t);
    
    // Step forward in time with RK4
    runge_kutta2(mesh, flux_calculator, t, dt);
    
    // Update node values with triangle-averaged values
    post_process(mesh);

    // get next time step
    dt = get_time_step(min_edge_length,mesh);
    
    // Update the viewer with new node positions
    viewer.add_nodes(mesh.node_begin(), mesh.node_end(),
                     HeatMap(cmin,cmax), CS207::DefaultPosition(), node_map);
    //viewer.center_view();
    viewer.set_label(t);

    // slow down the animation for small meshes
    if (mesh.num_nodes() < 100)
      CS207::sleep(0.05);
  }

  // delete listener
  //delete l;
  
  return 0;
}
