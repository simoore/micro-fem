import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from .mesh import UniformMesh
from .laminate_model import LaminateModel
from .laminate_dof import LaminateDOF


class LaminateFEM(object):

    def __init__(self, material, cantilever):
        
        self.cantilever = cantilever
        self.mesh = UniformMesh(cantilever.topology)
        self.dof = LaminateDOF(self.mesh)
        self.model = LaminateModel(material, cantilever.a, cantilever.b)
        self.a = cantilever.a
        self.b = cantilever.b
        self.assemble()
        
    
    def get_mass_matrix(self, free=False):
        
        muu = self.muu.tocsr()
        if free is False:
            return muu
        return muu[self.dof.free_dofs, :][:, self.dof.free_dofs]

    
    def get_stiffness_matrix(self, free=False):
        
        kuu = self.kuu.tocsr()
        if free is False:
            return kuu
        return kuu[self.dof.free_dofs, :][:, self.dof.free_dofs]
    
    
    def get_piezoelectric_matrix(self, free=False):
        
        kuv = self.kuv.tocsr()
        if free is False:
            return kuv
        return kuv[self.dof.free_dofs, :]
    
    
    def get_capacitance_matrix(self):
        return self.kvv
        
        
    def modal_analysis(self, n_modes):
        """The return value (w) are the eigenvalues and the return value (v) 
        are the eigenvectors.
        """
        m = self.muu.tocsc()[self.dof.free_dofs, :][:, self.dof.free_dofs]
        k = self.kuu.tocsc()[self.dof.free_dofs, :][:, self.dof.free_dofs]
        w, v = linalg.eigsh(k, k=n_modes, M=m, sigma=0, which='LM')
        vall = np.zeros((self.dof.n_mdof, n_modes))
        vall[self.dof.free_dofs, :] = v
        return w, v, vall
        
    
    def assemble(self):
        """The mass, stiffness, piezoelectric, and capacitance matricies are 
        assembled in this function.
        """
        muue = self.model.get_mass_element()
        kuue = self.model.get_stiffness_element()
        kuve = self.model.get_piezoelectric_element()
        kvve = self.model.get_capacitance_element()
        
        nm, ne = kuve.shape
        
        k_num = nm * nm * self.mesh.n_elem
        p_num = nm * ne * self.mesh.n_elem
        c_num = ne * ne * self.mesh.n_elem
        
        k_index = list(np.ndindex(nm, nm))
        p_index = list(np.ndindex(nm, ne))
        c_index = list(np.ndindex(ne, ne))
        
        k_row = np.zeros(k_num)
        k_col = np.zeros(k_num)
        k_val = np.zeros(k_num)
        m_val = np.zeros(k_num)
        p_row = np.zeros(p_num)
        p_col = np.zeros(p_num)
        p_val = np.zeros(p_num)
        c_row = np.zeros(c_num)
        c_col = np.zeros(c_num)
        c_val = np.zeros(c_num)
        
        k_ntriplet = 0
        p_ntriplet = 0
        c_ntriplet = 0
        
        for ni, e in enumerate(self.dof.dof_elements):
            
            for ii, jj in k_index:
                k_row[k_ntriplet] = e.mechanical_dof[ii]
                k_col[k_ntriplet] = e.mechanical_dof[jj]
                k_val[k_ntriplet] = kuue[ii, jj]
                m_val[k_ntriplet] = muue[ii, jj]
                k_ntriplet += 1
            
            for ii, jj in p_index:
                p_row[p_ntriplet] = e.mechanical_dof[ii]
                p_col[p_ntriplet] = e.electrical_dof[jj]
                p_val[p_ntriplet] = kuve[ii, jj]
                p_ntriplet += 1
            
            for ii, jj in c_index:
                c_row[c_ntriplet] = e.electrical_dof[ii]
                c_col[c_ntriplet] = e.electrical_dof[jj]
                c_val[c_ntriplet] = kvve[ii, jj]
                c_ntriplet += 1
        
        muu_shape = (self.dof.n_mdof, self.dof.n_mdof)
        kuu_shape = (self.dof.n_mdof, self.dof.n_mdof)
        kuv_shape = (self.dof.n_mdof, self.dof.n_edof)
        kvv_shape = (self.dof.n_edof, self.dof.n_edof)
        
        self.muu = sparse.coo_matrix((m_val, (k_row, k_col)), shape=muu_shape)
        self.kuu = sparse.coo_matrix((k_val, (k_row, k_col)), shape=kuu_shape)
        self.kuv = sparse.coo_matrix((p_val, (p_row, p_col)), shape=kuv_shape)
        self.kvv = sparse.coo_matrix((c_val, (c_row, c_col)), shape=kvv_shape)
        