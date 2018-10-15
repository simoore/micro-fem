import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from mesh_v2 import UniformMesh
from laminate_model import LaminateModel
from laminate_dof import LaminateDOF


class LaminateFEM(object):
    """
    Public Attributes
    -----------------
    :self.cantilever:
    :self.mesh: The collections of elements and their associated nodes that 
                define the domain over which the finite element analysis is 
                applied.
    """
    #def __init__(self, cantilever, material, to_connect=False, pmu=0.0):
    def __init__(self, cantilever, material):
        
        #self.to_connect = to_connect
        #self.pmu = pmu
        self.cantilever = cantilever
        self.mesh = UniformMesh(cantilever)
        self.dof = LaminateDOF(self.mesh)
        self.model = LaminateModel(material, cantilever.a, cantilever.b)
        self.assemble(self.mesh.get_densities())
        
    
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
        return w, vall
        
    
#    def set_penalty(self, xs, mu=None):
#        """self.piezo_temp_grad is the coefficient for df1/dmu
#        """
#        self.density = xs
#                
#        self.density_penalty = xs ** 3 
#        self.elastic_penalty = xs
#        self.piezo_penalty = xs ** 5
#        self.cap_penalty = xs
#        
#        self.density_grad = 3 * xs ** 2
#        self.elastic_grad = np.ones(xs.shape)
#        self.piezo_grad = 5 * xs ** 4
#        self.cap_grad = np.ones(xs.shape)
#
#        if self.to_connect is True:
#            mu = np.zeros_like(xs) if mu is None else mu
#            mu_scale = np.exp(-self.pmu * mu)
#            self.piezo_penalty = self.piezo_penalty * mu_scale
#            self.piezo_grad = self.piezo_grad * mu_scale
#            self.piezo_temp_grad = -self.pmu * self.piezo_penalty
#        else:
#            self.piezo_temp_grad = np.zeros_like(self.piezo_penalty)
        
    
    #def assemble(self, xs=None, mu=None):
    def assemble(self):
        """The mass, stiffness, piezoelectric, and capacitance matricies are 
        assembled in this function.
        """
        #if xs is None:
        #    return
        
        #self.set_penalty(xs, mu)
        
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
                #k_val[k_ntriplet] = self.elastic_penalty[ni] * kuue[ii, jj]
                #m_val[k_ntriplet] = self.density_penalty[ni] * muue[ii, jj]
                k_val[k_ntriplet] = kuue[ii, jj]
                m_val[k_ntriplet] = muue[ii, jj]
                k_ntriplet += 1
            
            for ii, jj in p_index:
                p_row[p_ntriplet] = e.mechanical_dof[ii]
                p_col[p_ntriplet] = e.electrical_dof[jj]
                #p_val[p_ntriplet] = self.piezo_penalty[ni] * kuve[ii, jj]
                p_val[p_ntriplet] = kuve[ii, jj]
                p_ntriplet += 1
            
            for ii, jj in c_index:
                c_row[c_ntriplet] = e.electrical_dof[ii]
                c_col[c_ntriplet] = e.electrical_dof[jj]
                #c_val[c_ntriplet] = self.cap_penalty[ni] * kvve[ii, jj]
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
        