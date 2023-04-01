# translated from BASIC code in Ohnishi, et al.,
# "Flow Analysis by a personal computer" chapter 5 (in Japanese)

import numpy as np
from scipy.linalg import solveh_banded,norm
import matplotlib.pyplot as plt

def fix_bc(A,B,node,val):
    """ take care of dirichlet boundary condition """
    for n,v in zip(node,val):
        B[:] -= A[:,n]*v; B[n] = v
        A[:,n] = 0; A[n,:] = 0; A[n,n] = 1

def solveh_banded_(A,B,w):
    """ helper function for solving banded matrix eq """
    n = len(B)
    AB = np.empty((w+1,n))
    for j in range(w+1):# construct banded matrix
        AB[w-j,j:] = A[np.r_[j:n], np.r_[:n-j]]
    return solveh_banded(AB,B)

class NSFEM:
    """ Finite Element Method for solving
    Navier-Stokes equation in two dimensions.
    """
    def __init__(self, node, element,
                 stream_fix = None,
                 vortex_fix = None,
                 vortex_wall = None,
                 stream_free = None,
                 vortex_free = None,
                 visc=0, scale=1, offset=0):
        """
        node : string, file name of node data
          column 1 = x-coordinate of nodes
          column 2 = y-corrdinate of nodes
          column 3 = initial values of vortex
          if column3 is absent, initial values are all set to 0
        element : string, file name of element data
          three columns are node indices of three nodes in elements
          node indices begin from 0
        stream_fix : string
          file name of dirichlet boundary condition for stream function
          column 1 = indices of nodes on boundary
          column 2 = value of stream function on boundary
        vortex_fix : string
          file name of dirichlet boundary condition for vortex
          column 1 = indices of nodes on boundary
          column 2 = value of vortex on boundary
        vortex_wall : string
          file name of wall boundary condition for vortex
          column 1 = indices of nodes on boundary
          column 2,3 = indices of adjacent nodes
          if column 3 < 0, only column 2 is used
          else, mid point of 2 and 3 are used as adjacent node
          column 4 = wall slipping velocity (rightward from inside)
        stream_free : string
          file name of neumann boundary condition for stream function
          column 1 = indices of nodes on boundary
          column 2 = derivative of stream function on boundary
        vortex_free : string
          file name of neumann boundary condition for vortex
          column 1 = indices of nodes on boundary
          column 2 = derivative of vortex on boundary
        visc : viscousity coefficient
          used only if vortex_free is not None
        scale : scale factor to be muliplied to xy-coordinates
        offset : integer to be subtracted from node indices
        """
        # read data file
        node = np.loadtxt(node)
        nn = len(node) # number of nodes
        if node.shape[1] > 2:
            self.vortex = node[:,2] # initial values of vortex
            node = node[:,:2]
        else: self.vortex = np.zeros(nn)
        node *= scale

        elem = np.loadtxt(element).astype(np.int) - offset

        if stream_fix is None:
            self.ps_bc1_node, self.ps_bc1_val = [],[]
        else:
            a = np.loadtxt(stream_fix)
            self.ps_bc1_node = a[:,0].astype(np.int) - offset
            self.ps_bc1_val = a[:,1]

        if vortex_fix is None:
            self.om_bc1_node, self.om_bc1_val = [],[]
        else:
            a = np.loadtxt(vortex_fix)
            self.om_bc1_node = a[:,0].astype(np.int) - offset
            self.om_bc1_val = a[:,1]

        if vortex_wall is None:
            self.om_bc3_node, self.om_bc3_vs = [],[]
        else:
            a = np.loadtxt(vortex_wall)
            self.om_bc3_node = a[:,:3].astype(np.int) - offset
            self.om_bc3_vs = a[:,3] # slipping velocity

        if stream_free is None:
            self.ps_bc2_node, self.ps_bc2_val = [],[]
        else:
            a = np.loadtxt(stream_free)
            self.ps_bc2_node = a[:,:2].astype(np.int) - offset
            self.ps_bc3_val = a[:,2]

        if vortex_free is None:
            self.om_bc2_node, self.om_bc2_val = [],[]
        else:
            a = np.loadtxt(vortex_free)
            self.om_bc2_node = a[:,:2].astype(np.int) - offset
            self.om_bc2_val = a[:,2]

        # width of band matrix
        self.w = np.max(np.abs(np.diff(np.c_[elem, elem[:,0]])))

        # det = 2*(area of element) 
        xy = node[elem]
        dxy = np.roll(xy, -1, axis=1) - np.roll(xy, 1, axis=1)
        det = np.einsum('ij,ij->i', xy[:,:,0], dxy[:,:,1])

        # D = laplacian matrix, M = mass matrix
        cb = np.einsum('kij,k->kij', -dxy, 1/det) # c and -b
        bc = np.einsum('kil,kjl,k->kijl', cb, cb, det)/2
        D = np.sum(bc, axis=-1)
        M = (np.ones((3,3)) + np.eye(3))/24
        self.D = np.zeros((nn,nn))
        self.M = np.zeros((nn,nn))
        for k,e in enumerate(elem):
            i,j = np.meshgrid(e,e,indexing='ij')
            self.D[i,j] += D[k]
            self.M[i,j] += det[k]*M

        # distance between nodes for wall boundary condition
        self.om_bc3_dist = np.empty(len(self.om_bc3_node))
        for l,(i,j,k) in enumerate(self.om_bc3_node):
            if k<0: d = norm(node[j] - node[i])
            else: d = norm((node[j] + node[k])/2 - node[i])
            self.om_bc3_dist[l] = d

        # neumann boundary condition for stream function
        self.G2 = np.zeros(nn)
        for n,v in zip(self.ps_bc2_node, self.ps_bc2_val):
            self.G2[n] = v*norm(np.diff(node[n]))/2

        # neumann boundary condition for vortex
        self.G4 = np.zeros(nn)
        for n,v in zip(self.om_bc2_node, self.om_bc2_val):
            self.G4[n] = visc*v*norm(np.diff(node[n]))/2

        self.node = node
        self.element = elem
        self.cb = cb
        self.bc = np.flip(bc, axis=-1) # bb and cc
        self.visc = visc
        
    def run(self, dt, visc=None, eps=0, maxit=100, om=None):
        """
        dt : time step
        visc : viscousity coefficient (= 1/Re)
          if None, visc in __init__() is used
        eps : convergence criterion such that
          if relative change in max(abs(vortex)) < eps then exit
        maxit : maximum number of iterations
        om : initial values of vortex on nodes
          if None, initial values are all set to 0
          if scalar, initial values are all set to om
        return number of iterations
        if "ValueError: array must not contain infs or NaNs"
          occurs, try again
        """
        if visc is None: visc = self.visc
        nn = len(self.node) # number of nodes
        nw = len(self.om_bc3_node) # number of wall boundary conditions
        om_b = np.empty(nw) # vortex at wall boundary
        if om is None: om = self.vortex # vortex
        elif np.isscalar(om): om = np.full(nn, om)
        Mdt = self.M/dt

        for n in range(maxit):
            # poisson equation for stream function
            A = self.D.copy() # LHS
            B = np.dot(self.M, om) + self.G2 # RHS

            # dirichlet boundary condition
            fix_bc(A, B, self.ps_bc1_node, self.ps_bc1_val)

            # stream function
            ps = solveh_banded_(A, B, self.w)

            # compute velocity field (u,v)
            pe = ps[self.element]
            uv = np.einsum('ijk,ij->ik', self.cb, pe)

            # diffusion and convection terms
            A_ = visc + uv**2*dt/2 # numerical viscousity
            A_ = np.einsum('kl,kijl->kij', A_, self.bc)
            B_ = (np.roll(pe, -1, axis=-1) - np.roll(pe, 1, axis=-1))/6

            # vortex transport equation
            A = Mdt.copy()
            B = Mdt.copy()
            for k,e in enumerate(self.element):
                i,j = np.meshgrid(e,e,indexing='ij')
                A[i,j] += A_[k] # LHS
                B[i,j] -= B_[k]
    
            B = np.dot(B, om) + self.G4 # RHS

            # dirichlet boundary condition
            fix_bc(A, B, self.om_bc1_node, self.om_bc1_val)

            # wall boundary condition
            for l,(i,j,k) in enumerate(self.om_bc3_node):
                d = self.om_bc3_dist[l]
                vs = self.om_bc3_vs[l] # slipping velocity
                if k<0: ps_n,om_n = ps[j], om[j]
                else:
                    ps_n = (ps[j] + ps[k])/2
                    om_n = (om[j] + om[k])/2

                om_b[l] = 3*((ps[i] - ps_n)/d - vs)/d - om_n/2

            if nw: fix_bc(A, B, self.om_bc3_node[:,0], om_b)

            # votex
            om_ = om
            om = solveh_banded_(A, B, self.w)
            if np.max(np.abs(om-om_)) < eps*np.max(np.abs(om_)):
                break # converged

        self.stream = ps
        self.vortex = om
        self.velocity = uv
        return n+1

    def plot_element(self, *a, **k):
        x,y = self.node.T
        for e in self.element:
            e = np.r_[e,e[0]]
            plt.plot(x[e], y[e], *a, **k)

    def plot_velocity(self, *a, **k):
        x,y = np.sum(self.node[self.element], axis=1).T/3
        u,v = self.velocity.T
        plt.quiver(x,y,u,v, *a, **k)

    def save(self, node, element,
             stream_fix = None,
             vortex_fix = None,
             vortex_wall = None,
             stream_free = None,
             vortex_free = None,
             save_vortex = True):
        """ if save_vortex is True,
        vortex is saved in column 3 of node file """
        if save_vortex:
            np.savetxt(node, np.c_[self.node,
                                   self.vortex], '%g %g %g')
        else: np.savetxt(node, self.node, '%g %g')

        np.savetxt(element, self.element, '%d %d %d')

        if stream_fix is not None:
            np.savetxt(stream_fix,
                       np.c_[self.ps_bc1_node,
                             self.ps_bc1_val], '%d %g')
        if vortex_fix is not None:
            np.savetxt(vortex_fix,
                       np.c_[self.om_bc1_node,
                             self.om_bc1_val], '%d %g')
        if vortex_wall is not None:
            np.savetxt(vortex_wall,
                       np.c_[self.om_bc3_node,
                             self.om_bc3_vs], '%d %d %d %g')
        if stream_free is not None:
            np.savetxt(stream_free,
                       np.c_[self.ps_bc2_node,
                             self.ps_bc2_val], '%d %g')
        if vortex_free is not None:
            np.savetxt(vortex_free,
                       np.c_[self.om_bc2_node,
                             self.om_bc2_val], '%d %g')
