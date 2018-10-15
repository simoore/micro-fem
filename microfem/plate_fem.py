from collections import namedtuple
from matplotlib import patches, cm
from matplotlib.pyplot import figure, show
from mpl_toolkits.mplot3d import Axes3D
from numpy import array, zeros, dot, hstack, meshgrid, arange, full, nan, pi
from numpy import sqrt
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

Options = namedtuple('FemOptions', 'Rho E Nu H A B WID LEN')

class Element:
    def __init__(self, i, j):
        self.i = i
        self.j = j
        self.nodes = [None for _ in range(4)]
        self.void = True


class Node:
    def __init__(self, node_id, i, j):
        self.id = node_id
        self.i = i
        self.j = j
        self.boundary = False





def mesh_add_node(i, j, node, corner, elements, opt):
    """Adds a node to an elements. The coordinates are checked to ensure they
    are valid.
    
    Parameters
    ----------
    i : int
        x-coordinate of the element.
    j : int
        y-coordinate of the element.
    node : Node
        The node to be added to the element.
    corner : int in [0, 1, 2, 3]
        The index of were to add the node to the elements node list.
    elements : list of Element
        The list of elements.
    opt : Options
        The finite element parameters.
    """
    WID, LEN = opt.WID, opt.LEN
    if 0 <= i < WID and 0 <= j < LEN:
        e = elements[i][j]
        e.nodes[corner] = node


def mesh_element(i, j, elements, nodes, bnodes, opt):
    """
    """
    # The key is to understand that their are four nodes in each element.
    # The nodes are index (0,1,2,3). They correspond to the lower-left, 
    # lower-right, upper-right and upper-left corners. The constant arrays
    # are associated with each corner. 'ni' and 'nj' are used to convert 
    # element coordinates to node coordinates. 
    def corner_func(x): return x if x < 4 else x - 4
    ni = [0, 1, 1, 0]   
    nj = [0, 0, 1, 1]
    ii = [[-1, -1, 0], [0, 1, 1], [1, 1, 0], [0, -1, -1]]
    jj = [[0, -1, -1], [-1, -1, 0], [0, 1, 1], [1, 1, 0]]
    
    e = elements[i][j]
    e.void = False
    
    # TODO: Add more flexible boundary conditions.
    new_nodes = [x for x in range(4) if e.nodes[x] is None]
    for k in new_nodes:
        node_id = len(bnodes) if (j == 0) and (k in [0, 1]) else len(nodes)
        node = Node(node_id, i + ni[k], j + nj[k])
        if node.j == 0:
            node.boundary = True
            bnodes.append(node)
        else:
            nodes.append(node)
        e.nodes[k] = node
        
        # Add the new node to the three neighbouring elements.
        for neighbour in range(3):
            iii = i + ii[k][neighbour]
            jjj = j + jj[k][neighbour]
            corner = corner_func(k + neighbour + 1)
            mesh_add_node(iii, jjj, node, corner, elements, opt)


def mesh(topology, opt):
    """Generates a finite element mesh from the cantilever topology.
    
    Parameters
    ----------
    topololgy : list of list of int
        The binary matrix indicating whether an element is present or void.
    opt : Options
        The finite element parameters.
    """
    WID, LEN = opt.WID, opt.LEN
    elements = [[Element(i, j) for j in range(LEN)] for i in range(WID)]
    nodes, bnodes = [], []
    for i in range(WID):
        for j in range(LEN):
            if topology[i][j] == 1:
                mesh_element(i, j, elements, nodes, bnodes, opt)
    elements = [x for flat in elements for x in flat if x.void is False]
    
    # Find tip node.
    for n in nodes:
        if n.i == int(WID/2) and n.j == LEN:
            tip_id = n.id
            
    return elements, nodes, bnodes, tip_id


def assemble(topology, opt):
    """
    """
    elements, nodes, _, tip_id = mesh(topology, opt)
    ke, me = element(opt)
    ndof = 3
    ntriplet = 0
    
    # 'num' is the maximum possible number of triplets needed to create the 
    # system matrices.
    num = 144 * len(elements)
    row, col, kk, mm = zeros(num), zeros(num), zeros(num), zeros(num)
    
    for e in elements:
        dof = [ndof * n.id + y for n in e.nodes for y in range(ndof)]
        boundary = [n.boundary for n in e.nodes for _ in range(ndof)]
        for ii in range(12):
            for jj in range(12):
                if boundary[ii] is False and boundary[jj] is False:
                    row[ntriplet] = dof[ii]
                    col[ntriplet] = dof[jj]
                    mm[ntriplet] = me[ii, jj]
                    kk[ntriplet] = ke[ii, jj]
                    ntriplet += 1
                    
    # Create the sparse mass and stiffness matrices.
    dimen = ndof * len(nodes)
    kks = coo_matrix((kk, (row, col)), shape=(dimen, dimen)).tocsr()
    mms = coo_matrix((mm, (row, col)), shape=(dimen, dimen)).tocsr()
    tip_dof = ndof*tip_id
    return kks, mms, tip_dof


def modal(num_modes, topology, opt):
    
    kks, mms, tip_dof = assemble(topology, opt)
    w, v = eigsh(kks, k=num_modes, M=mms, sigma=0, which='LM')
    return w, v, kks, mms, tip_dof


def plot_mesh(topology, opt):
    elements, _, _, _ = mesh(topology, opt)
    fig = figure()
    subplot = fig.add_subplot(111, aspect='equal')
    for e in elements:
        node = e.nodes[0]
        x = node.i
        y = node.j
        width, height = 1, 1
        subplot.add_patch(patches.Rectangle((y, x), width, height))
    subplot.autoscale_view(True, True, True)
    show()


def plot_topology(topology):
    for r in topology:
        for e in r:
            print(int(e), end="")
        print()


def plot_mode(topology, v, mode_num, opt):
    _, nodes, bnodes, _ = mesh(topology, opt)
    xlim = max([x.i for x in nodes]) + 1
    ylim = max([x.j for x in nodes]) + 1
    ze = full((xlim, ylim), nan)
    zz = full((xlim, ylim), nan)
    w = v[0::3, mode_num]
    x = arange(0, xlim, 1)
    y = arange(0, ylim, 1)
    x, y = meshgrid(y, x)
    for n in nodes:
        ze[n.i, n.j] = 0
        zz[n.i, n.j] = w[n.id]
    for n in bnodes:
        ze[n.i, n.j] = zz[n.i, n.j] = 0
    fig = figure()
    axis = fig.add_subplot(111, projection='3d')
    axis.plot_surface(x, y, zz, rstride=1, cstride=1, linewidth=0, cmap=cm.jet,
                      vmin=min(w), vmax=max(w))
    axis.plot_surface(x, y, ze, rstride=1, cstride=1, linewidth=0, alpha=0.3)
    show()


def identification(num_modes, topology, opt):
    l, _, _, _, _ = modal(num_modes, topology, opt)
    freq = sqrt(l) / 2 / pi
    for i in range(num_modes):
        print('mode: %g' % freq[i])
    for i in range(1, num_modes):
        print('ratio %g' % (freq[i] / freq[0]))
