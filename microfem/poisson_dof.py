import numpy as np


class PoissonDOF(object):
    
    def __init__(self, mesh):

        self.mesh = mesh
        self.dof_nodes = [PoissonNode(n) for n in mesh.nodes]
        self.dof_elements = [PoissonElement(e, self.dof_nodes) 
                             for e in mesh.elements]
        
        self.all_dofs = np.arange(mesh.n_node)
        self.fixed_dofs = [n.dof for n in self.dof_nodes if n.node.j == 0]
        self.free_dofs = np.setdiff1d(self.all_dofs, self.fixed_dofs)
        self.n_dof = len(self.all_dofs)
        
        
class PoissonElement(object):
    
    def __init__(self, element, dof_nodes):
        
        self.element = element
        self.dof_nodes = [dof_nodes[n.index] for n in element.nodes]
        self.dofs = np.array([n.dof for n in self.dof_nodes])
    
    
class PoissonNode(object):
    
    def __init__(self, node):
        
        self.node = node
        self.dof = node.index