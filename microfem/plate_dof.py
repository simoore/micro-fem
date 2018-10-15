import numpy as np


class PlateDOF(object):
    
    def __init__(self, mesh):
        
        self.mesh = mesh
        self.dof_nodes = [PlateNode(n) for n in mesh.nodes]
        self.dof_elements = [PlateElement(e, self.dof_nodes) for e in mesh.elements]
        
        ads = [n.mechanical_dof for n in self.dof_nodes]
        fds = [n.mechanical_dof for n in self.dof_nodes if n.boundary == True]
        
        self.all_dofs = np.concatenate(ads)
        self.fixed_dofs = np.concatenate(fds)
        self.free_dofs = np.setdiff1d(self.all_dofs, self.fixed_dofs)

        self.n_mdof = len(self.all_dofs)
        self.n_elem = len(self.dof_elements)
        
        
class PlateElement(object):
    
    def __init__(self, element, dof_nodes):

        self.element = element
        self.dof_nodes = [dof_nodes[n.index] for n in element.nodes]
        
        mds = [n.mechanical_dof for n in self.dof_nodes]
        self.mechanical_dof = np.concatenate(mds)

        
    def get_displacement(self, u):
        """
        Displacement (u) is the system wide displacement field and this
        function extracts the components related to this element.
        """
        ue = np.zeros((len(self.mechanical_dof), 1))
        ue[:] = u[self.mechanical_dof]
        return ue


class PlateNode(object):
    """
    Plates have three mechanical DOFs per node.
    """
        
    def __init__(self, node):
        
        self.node = node
        self.mechanical_dof = np.arange(3) + 3 * node.index
        self.deflection_dof = self.mechanical_dof[0]
        self.boundary = True if node.j == 0 else False
