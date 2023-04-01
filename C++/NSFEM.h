// Finite Element Method for solving
//   Navier-Stokes equation in two dimensions.

#ifndef __NSFEM_h__
#define __NSFEM_h__

#include "Mat.h"

namespace NSFEM {
    void init(const char *node_file,
              const char *element_file,
              const char *stream_fix=0,
              const char *vortex_fix=0,
              const char *vortex_wall=0,
              const char *stream_free=0,
              const char *vortex_free=0,
              double visc=0., double scale=1., int offset=0);
    // node_file : file name of node data
    //   column 1 = x-coordinate of nodes
    //   column 2 = y-corrdinate of nodes
    //   column 3 = initial values of vortex
    //   if column 3 is absent, initial values are all set to 0
    // element_file : file name of element data
    //   three columns are node indices of three nodes in elements
    //   node indices begin from 0
    // stream_fix :
    //   file name of dirichlet boundary condition for stream function
    //   column 1 = indices of nodes on boundary
    //   column 2 = value of stream function on boundary
    // vortex_fix :
    //   file name of dirichlet boundary condition for vortex
    //   column 1 = indices of nodes on boundary
    //   column 2 = value of vortex on boundary
    // vortex_wall :
    //   file name of wall boundary condition for vortex
    //   column 1 = indices of nodes on boundary
    //   column 2,3 = indices of adjacent nodes
    //   if column 3 < 0, only column 2 is used
    //   else, mid point of 2 and 3 are used as adjacent node
    //   column 4 = wall slipping velocity (rightward from inside)
    // stream_free :
    //   file name of neumann boundary condition for stream function
    //   column 1 = indices of nodes on boundary
    //   column 2 = derivative of stream function on boundary
    // vortex_free :
    //   file name of neumann boundary condition for vortex
    //   column 1 = indices of nodes on boundary
    //   column 2 = derivative of vortex on boundary
    // visc : viscousity coefficient
    //   used only if vortex_free is not 0
    // scale : scale factor to be muliplied to xy-coordinates
    // offset : integer to be subtracted from node indices

    int run(double dt, double visc=0., double eps=0., int maxit=100);
    //  dt : time step
    //  visc : viscousity coefficient (= 1/Re)
    //    if visc <= 0, visc in init() is used
    //  eps : convergence criterion such that
    //    if relative change in max(abs(vortex)) < eps then exit
    //  maxit : maximum number of iterations
    //  return number of iterations

    void elem_center(double& x, double& y, int i);
    // x,y = coordinates of mass center of element i
    
    extern mat_double node; // xy-coordinates of nodes
    extern mat_int element; // node indices of element vertices
    extern vec_double stream; // stream function at nodes
    extern vec_double vortex; // vortex at nodes
    extern mat_double velocity; // velocity at elements
}

#endif// __NSFEM_h__
