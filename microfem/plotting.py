import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.tri as mtri
import matplotlib.patches as patches
import scipy.interpolate as interpolate


def plot_topology(cantilever, filename=None):
    """
    Plots a black and white image of the topology of the cantilever where
    black is a solid element and white is a void element. 
    
    Parameters
    ----------
    cantilever : Cantilever object  
        The cantilever object from which contains the topology.
    filename : string    
        The filename for the file to save the plot. If None, no file is saved.
    """
    
    nelx, nely = cantilever.topology.shape
    a = cantilever.a
    b = cantilever.b
    
    x = np.linspace(0, nelx*2*a, nelx+1, endpoint=True)
    y = np.linspace(0, nely*2*b, nely+1, endpoint=True)
    xv, yv = np.meshgrid(x, y)
    
    fig, ax = plt.subplots()
    ax.pcolormesh(xv, yv, cantilever.topology.T, cmap=cm.Greys, vmin=0, vmax=1)
    ax.set_xlim(0, nelx*2*a)
    ax.set_ylim(0, nely*2*b)
    ax.set_aspect('equal')
    plt.show()
    if filename is not None:
        #plt.tight_layout(pad=0.0)
        fig.savefig(filename, dpi=300, bbox_inches='tight')
        
        
def plot_mode(fem, v):
        
    # Create list of points at each node and the deflection DOF.
    n_pts = len(fem.dof.dof_nodes)
    pts = np.zeros((n_pts, 3))
    for n in fem.dof.dof_nodes:
        z = v[n.deflection_dof]
        pts[n.node.index, :] = np.array([n.node.i, n.node.j, z])
    
    # Create two triangles from each element.
    n_tri = round(2 * fem.dof.n_elem)
    tri = np.zeros((n_tri, 3))
    for e in fem.dof.dof_elements:
        index1 = round(2 * e.element.index)
        index2 = round(2 * e.element.index + 1)
        n0, n1, n2, n3 = e.element.nodes
        tri[index1, :] = np.array([n0.index, n1.index, n2.index])
        tri[index2, :] = np.array([n2.index, n3.index, n0.index])

    triang = mtri.Triangulation(pts[:, 0], pts[:, 1], triangles=tri)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_trisurf(triang, pts[:, 2], cmap=cm.viridis)
    ax.plot_trisurf(triang, np.zeros(n_pts), alpha=0.3)
    ax.view_init(30, 45)
    plt.show()
    
    
    
def plot_poisson_solution(fem, uall, show=None):
    
    assert len(fem.dof.dof_nodes) == uall.shape[0]

    a = fem.poisson_domain.a
    b = fem.poisson_domain.b
    nelx, nely = fem.poisson_domain.domain.shape
    
    x = np.arange(nelx + 1) * 2 * a
    y = np.arange(nely + 1) * 2 * b
    mx, my = np.meshgrid(x, y, indexing='ij')
    mask = np.ones(mx.shape, dtype=bool)
    zdata = np.zeros(mx.shape)
    for n in fem.dof.dof_nodes:
        mask[n.node.i, n.node.j] = False
        zdata[n.node.i, n.node.j] = uall[n.dof]
    z_mask = np.ma.array(zdata, mask=mask)

    fig, ax = plt.subplots()
    
    # Add information about the poisson.
    extent = [0, 2*nelx*a, 0, 2*nely*b]
    if show == 'domain':
        img = fem.poisson_domain.domain.T
    elif show == 'conductivity':
        img = fem.poisson_domain.conductivity.T
    elif show == 'source':
        img = fem.poisson_domain.source.T
    ax.imshow(img, extent=extent, origin='lower', cmap=cm.Pastel1)  
    
    # Add the temperature distribution.
    cs = ax.contour(mx, my, z_mask, 10, colors='k', interpolation='bilinear')
    ax.clabel(cs, inline=1, fontsize=10)
    ax.set_title('Temperature Distribution')
    ax.set_xlabel('Width (um)')
    ax.set_ylabel('Length(um)')
    
    plt.show()
