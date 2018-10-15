import numpy as np


class UniformMesh(object):
    """The optimization problem is defined for a fixed rectangular mesh. The 
    nodes along the x-axis (y=0) are on a clampled boundary. This class stores
    the list of elements that make up the domain and associates the nodes with 
    their corresponding elements. Since several different finite element 
    analyses are generated from the same domain, other classes associate the
    degrees of freedom (dof) with each node using the geometric information
    from this class. This class does manage the density property utilized in
    topology optimization routines.
    
    Public Attributes
    -----------------    
    self.elements : list of objects
    The list of elements that make up the cantilever domain.
    
    self.nodes : list of objects
    The list of nodes that comprise the canitlever domain.
    """
    
    def __init__(self, cantilever):
        """
        :param cantilever: An object describing a the topology of the 
            cantilever. This requires a binary matrix called `topology` to 
            describe which element in the mesh exists.
        """
        
        # Create all elements and nodes on the rectangular domain.
        nelx, nely = cantilever.topology.shape
        self._nodes_2D = [[Node(i, j) for j in range(nely + 1)]
                          for i in range(nelx + 1)]
        self._elements_2D = [[Element(i, j, cantilever, self._nodes_2D) 
                             for j in range(nely)] for i in range(nelx)]

        # Create list of valid elements and nodes.
        gen_elements = (e for row in self._elements_2D for e in row)
        gen_nodes = (n for row in self._nodes_2D for n in row)
        self.elements = [e for e in gen_elements if e.void == False]
        self.nodes = [n for n in gen_nodes if n.void == False]
        
        
        # Set index on elements.
        for i, e in enumerate(self.elements):
            e.index = i
            
        for i, n in enumerate(self.nodes):
            n.index = i
        
        self.a = cantilever.a
        self.b = cantilever.b


    @property
    def n_node(self):
        return len(self.nodes)
    
    
    @property
    def n_elem(self):
        return len(self.elements)
        
        
    def get_densities(self):
        
        return np.array([e.density for e in self.elements])
    
    
    def set_densities(self, densities):
        
        for x, e in zip(densities, self.elements):
            e.density = x
    
    
    def to_console(self):
        
        print('-- Elements --')
        for e in self.elements:
            if e.void == False:
                print(e)
                
        print('\n-- Nodes --')
        for n in self.nodes:
            if n.void == False:
                print(n)
                    
                    
class Element(object):
    """
    Public Attributes
    -----------------
    self.i     : The x-index of the element.
    self.j     : The y-index of the element.
    self.nodes : A four element tuple that stores the nodes of the element. The 
                 element is rectangular with a node in each corner. The nodes, 
                 denoted by their position on a compass, are stored in order 
                 (sw, se, ne, nw).
    self.void  : False if a member of the domain, else True.
    self.index : The index in the list of non-void elements.
    """
    def __init__(self, i, j, cantilever, nodes_2D):
        
        nsw = nodes_2D[i][j]
        nse = nodes_2D[i + 1][j]
        nne = nodes_2D[i + 1][j + 1]
        nnw = nodes_2D[i][j + 1]
        
        # Public Attributes.
        self.i = i
        self.j = j
        self.nodes = (nsw, nse, nne, nnw)
        self.void = False if cantilever.topology[i][j] == 1 else True
        self.index = 0 # set later if not void
        #self.density = cantilever.densities[i][j]
        
        # Set node.void to False if element is non-void.
        if self.void is False:
            for n in self.nodes:
                n.void = False
            
        
    def __repr__(self):
        fields = tuple([self.index] + [n.index for n in self.nodes] + 
                       [self.i, self.j])
        repr_ = 'Element %d: %d %d %d %d (%g, %g)' % fields
        return repr_
    

class Node(object):
    """
    Public Attributes
    -----------------
    self.i     : The x-index of the node in the 2D list of nodes.
    self.j     : The y-index of the node in the 2D list of nodes.
    self.index : The index in the list of non-void nodes.
    self.void  : Indicates whether the node is part of the models domain.
    self.boundary : bool
    Indicates if the node lies on the boundary of the cantilever, that is 
    along the x-axis (y=0 or j=0).
    """
    def __init__(self, i, j):

        self.i = i
        self.j = j
        self.index = 0 # set later if void is false
        self.void = True
        self.boundary = True if j == 0 else False


    def __repr__(self):
        
        s1 = 'Node {0:d}: {1:d} {2:d}'.format(self.index, self.i, self.j)
        return s1
