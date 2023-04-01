// Finite Element Method for solving
//   Navier-Stokes equation in two dimensions.
// translated from BASIC code in Ohnishi, et al.,
//  "Flow Analysis by a personal computer" chapter 5 (in Japanese)

#include<fstream>
#include<cstdlib>
#include<cmath>
#include "NSFEM.h"
#include "MatBandSym.h"
#include "util.h"
using namespace NSFEM;

mat_double NSFEM::node; // xy-coordinates of nodes
mat_int NSFEM::element; // node indices of elements
vec_double NSFEM::stream; // stream function at nodes
vec_double NSFEM::vortex; // vortex at nodes
mat_double NSFEM::velocity; // velocity at elements

static vec_int ps_bc1_node;// dirichlet boundary nodes for stream function
static vec_int om_bc1_node;// dirichlet boundary nodes for vortex
static mat_int om_bc3_node;// wall boundary nodes for vortex
static vec_int om_bc3_node0;// wall boundary nodes for vortex
static vec_double ps_bc1_val;// dirichlet boundary values for stream function
static vec_double om_bc1_val;// dirichlet boundary values for vortex
static vec_double om_bc3_vs; // slipping velocity
static vec_double om_bc3_dist;// distance between nodes
static vec_double G2;// neumann boundary condition for stream function
static vec_double G4;// neumann boundary condition for vortex
static vec_double det;// twice area of elements

static int nn; // number of nodes
static int ne; // number of elements;
static int nw; // number of wall boundaries
static int bw; // half width of banded matrix
static double visc; // viscosity
static MatBandSym D; // laplacian matrix
static MatBandSym M; // mass matrix

#define roll1(i) (i<2 ? i+1:0)
#define roll2(i) (i>0 ? i-1:2)

static void bce(double *b, double *c, const int *e)
// b = difference of y-coordinates of nodes in element e
// c = difference of x-coordinates of nodes in element e
{
    int i,j,k;
    for(i=0; i<3; i++) {
        j = e[roll1(i)];
        k = e[roll2(i)];
        b[i] = node[j][1] - node[k][1];
        c[i] = node[k][0] - node[j][0];
    }
}

static double dist(int i, int j)
// distance between nodes i and j
{
    return sqrt(pow(node[i][0] - node[j][0], 2)
              + pow(node[i][1] - node[j][1], 2));
}
    
static double dist(int i, int j, int k)
// distance between node i and midpoint of nodes j and k
{
    double d;
    d = pow(node[i][0] - (node[j][0] + node[k][0])/2, 2);
    d+= pow(node[i][1] - (node[j][1] + node[k][1])/2, 2);
    return sqrt(d);
}

static void fix_bc(MatBandSym& A, vec_double& B,
                   const vec_int& node, const vec_double& val)
// modify A and B to take care of
//   dirichlet boundary condition in (node, val)
{
    int i,j,k,l,n,m(B.size());
    for(n=0; n<node.size(); n++) {
        j = node[n];
        k = max(j-bw, 0);
        l = min(j+bw+1, m);
        for(i=k; i<l; i++) B[i] -= A(i,j)*val[n];
        for(i=k; i<l; i++) A(i,j) = 0;
        A(j,j) = 1;
        B[j] = val[n];
    }
}

static void fix_bc(vec_double& B,
                   const vec_int& node, const vec_double& val)
// modify B to take care of
//   dirichlet boundary condition in (node, val)
{
    for(int i=0; i<node.size(); i++)
        B[node[i]] = val[i];
}

void loadtxt(mat_double& A, std::ifstream& f);

void NSFEM::init(const char *node_file,
                 const char *element_file,
                 const char *stream_fix,
                 const char *vortex_fix,
                 const char *vortex_wall,
                 const char *stream_free,
                 const char *vortex_free,
                 double visc_, double scale, int offset)
// initializer to read data files
{
    int i,j,k,n[3],*e;
    double d;
    mat_double A;
    std::ifstream f;

    // read node data
    f.open(node_file);
    loadtxt(A,f);
    f.close();
    nn = A.nrows(); // number of nodes
    node.resize(nn, 2);
    for(i=0; i<nn; i++)
        for(j=0; j<2; j++)
            node[i][j] = A[i][j]*scale;

    vortex.assign(nn, 0.);
    if(A.ncols()>2)// initial values of vortex
        for(i=0; i<nn; i++) vortex[i] = A[i][2];

    // read element data
    f.open(element_file);
    loadtxt(A,f);
    f.close();
    ne = A.nrows(); // number of elements
    element.resize(ne, 3);
    for(i=0; i<ne; i++)
        for(j=0; j<3; j++)
            element[i][j] = int(A[i][j]) - offset;

    velocity.resize(ne, 2);

    // read dirichlet boundary data for stream function
    f.open(stream_fix);
    loadtxt(A,f);
    f.close();
    ps_bc1_node.resize(A.nrows());
    ps_bc1_val.resize(A.nrows());
    for(i=0; i<A.nrows(); i++)
        ps_bc1_node[i] = int(A[i][0]) - offset;
    for(i=0; i<A.nrows(); i++)
        ps_bc1_val[i] = A[i][1];

    // read dirichlet boundary data for vortex
    f.open(vortex_fix);
    loadtxt(A,f);
    f.close();
    om_bc1_node.resize(A.nrows());
    om_bc1_val.resize(A.nrows());
    for(i=0; i<A.nrows(); i++)
        om_bc1_node[i] = int(A[i][0]) - offset;
    for(i=0; i<A.nrows(); i++)
        om_bc1_val[i] = A[i][1];

    // read wall boundary data for vortex
    f.open(vortex_wall);
    loadtxt(A,f);
    f.close();
    nw = A.nrows(); // number of wall boundaries
    om_bc3_node.resize(nw, 3);
    om_bc3_node0.resize(nw);
    om_bc3_dist.resize(nw);
    om_bc3_vs.resize(nw);
    for(i=0; i<nw; i++) {
        for(j=0; j<3; j++)
            n[j] = om_bc3_node[i][j] = int(A[i][j]) - offset;
        om_bc3_node0[i] = n[0];
        om_bc3_vs[i] = A[i][3];
        if(n[2]<0) om_bc3_dist[i] = dist(n[0], n[1]);
        else om_bc3_dist[i] = dist(n[0], n[1], n[2]);
    }

    // width of banded matrix
    for(i=0; i<ne; i++) {
        e = element[i];
        bw = max(bw, abs(e[0]-e[1]));
        bw = max(bw, abs(e[1]-e[2]));
        bw = max(bw, abs(e[2]-e[0]));
    }

    // D = laplacian matrix, M = mass matrix
    double b[3], c[3];
    det.resize(ne);
    D.resize(nn, bw, 0.);
    M.resize(nn, bw, 0.);
    for(k=0; k<ne; k++) {
        e = element[k];
        bce(b,c,e);
        det[k] = b[0]*c[1] - b[1]*c[0];
        for(i=0; i<3; i++) {
            for(j=i; j<3; j++) {
                D(e[i],e[j]) += (b[i]*b[j] + c[i]*c[j])/det[k]/2;
                M(e[i],e[j]) += (i==j ? 2:1)/24.*det[k];
            }
        }
    }

    // read neumann boundary data for stream function
    f.open(stream_free);
    loadtxt(A,f);
    f.close();
    G2.assign(nn, 0.);
    for(i=0; i<A.nrows(); i++) {
        for(j=0; j<2; j++) n[j] = int(A[i][j]) - offset;
        d = dist(n[0],n[1]);
        G2[n[0]] = G2[n[1]] = A[i][2]*d/2;
    }
    
    // read neumann boundary data for vortex
    f.open(vortex_free);
    loadtxt(A,f);
    f.close();
    visc = visc_;
    G4.assign(nn, 0.);
    for(i=0; i<A.nrows(); i++) {
        for(j=0; j<2; j++) n[j] = int(A[i][j]) - offset;
        d = dist(n[0],n[1]);
        G4[n[0]] = G4[n[1]] = visc*A[i][2]*d/2;
    }

    // LDLT decomposition of the matrix D
    fix_bc(D, G2, ps_bc1_node, ps_bc1_val);
    LDLTdecomp(D,D);
}

int NSFEM::run(double dt, double visc_, double eps, int maxit)
// run simulation
{
    int i,j,k,l,m,n,*e;
    double b[3], c[3], u, v, d, ps_n, om_n;
    MatBandSym A, Mdt(M);
    MatBand C;
    vec_double B, om, om_b(nw);

    if(visc_ > 0) visc = visc_;
    for(i=0; i<M.nrows(); i++)
        for(j=0; j<M.ncols(); j++) Mdt[i][j] /= dt;
    for(m=0; m<maxit; m++) {// main loop
        // solve poisson equation for stream function
        B = G2;
        MulAddTo(B, M, vortex);
        fix_bc(B, ps_bc1_node, ps_bc1_val);
        solve(stream, D, B, false); // D is already decomposed
        // solve vortex transport equation
        A = Mdt;
        conv(C,A);
        for(n=0; n<ne; n++) {
            e = element[n];
            bce(b,c,e);
            // compute velocity field
            u = v = 0.;
            for(i=0; i<3; i++) u += c[i]*stream[e[i]];
            for(i=0; i<3; i++) v -= b[i]*stream[e[i]];
            velocity[n][0] = (u /= det[n]);
            velocity[n][1] = (v /= det[n]);
            // numerical viscousity
            u = (visc + u*u*dt/2)/det[n]/2;
            v = (visc + v*v*dt/2)/det[n]/2;
            // compute diffusion terms
            for(i=0; i<3; i++)
                for(j=i; j<3; j++)
                    A(e[i],e[j]) += u*b[i]*b[j] + v*c[i]*c[j];
            // compute convection terms
            for(j=0; j<3; j++) {
                k = e[roll1(j)];
                l = e[roll2(j)];
                d = (stream[k] - stream[l])/6;
                for(i=0; i<3; i++) C(e[i],e[j]) -= d;
            }
        }
        // wall boundary condition for vortex
        for(i=0; i<nw; i++) {
            e = om_bc3_node[i];
            d = om_bc3_dist[i];
            u = om_bc3_vs[i]; // slipping velocity
            if(e[2]<0) {
                ps_n = stream[e[1]];
                om_n = vortex[e[1]];
            }
            else {
                ps_n = (stream[e[1]] + stream[e[2]])/2;
                om_n = (vortex[e[1]] + vortex[e[2]])/2;
            }
            om_b[i] = 3*((stream[e[0]] - ps_n)/d - u)/d - om_n/2;
        }
        B = G4;
        MulAddTo(B, C, vortex);
        fix_bc(A, B, om_bc1_node, om_bc1_val);
        fix_bc(A, B, om_bc3_node0, om_b);
        om = vortex;
        solve(vortex, A, B);
        // check convergence
        u = v = 0.;
        for(i=0; i<nn; i++) {
            u = max(u, fabs(om[i] - vortex[i]));
            v = max(v, fabs(om[i]));
        }
        if(u < eps*v) return m+1;
    }
    return m;
}

void NSFEM::elem_center(double& x, double& y, int i)
// x,y = coordinates of mass center of element i
{
    const int *e(element[i]);
    x = (node[e[0]][0] + node[e[1]][0] + node[e[2]][0])/3;
    y = (node[e[0]][1] + node[e[1]][1] + node[e[2]][1])/3;
}