import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg

from .plate_model import PlateModel
from .plate_dof import PlateDOF
from .mesh import UniformMesh


class PlateFEM(object):

    def __init__(self, material, cantilever):
        """
        The initialization rountine creates the element models. The mesh and 
        penalization are updated seperately.
        
        Parameters
        ----------
        material : microfem.PlateMaterial 
            An object containing the material properties of the plate.
        """
        
        self._model = PlateModel(material, cantilever.a, cantilever.b)
        self.a = cantilever.a
        self.b = cantilever.b
        self._mesh = UniformMesh(cantilever.topology)
        self.dof = PlateDOF(self._mesh)
        self._assemble()


    def get_mass_matrix(self, free=False):
        
        muu = self._muu.tocsr()
        if free is False:
            return muu
        return muu[self.dof.free_dofs, :][:, self.dof.free_dofs]

    
    def get_stiffness_matrix(self, free=False):
        
        kuu = self._kuu.tocsr()
        if free is False:
            return kuu
        return kuu[self.dof.free_dofs, :][:, self.dof.free_dofs]
    
    
    def modal_analysis(self, n_modes):
        """
        The return value (w) are the eigenvalues and the return value (v) 
        are the eigenvectors.
        """
        
        m = self._muu.tocsc()[self.dof.free_dofs, :][:, self.dof.free_dofs]
        k = self._kuu.tocsc()[self.dof.free_dofs, :][:, self.dof.free_dofs]
        w, v = linalg.eigsh(k, k=n_modes, M=m, sigma=0, which='LM')
        vall = np.zeros((self.dof.n_mdof, n_modes))
        vall[self.dof.free_dofs, :] = v
        return w, v, vall
    
    
    def _assemble(self):
        """
        Assembles the mass and stiffness matrix of the finite element model of 
        the plate.
        """
        
        muue = self._model.me
        kuue = self._model.ke
        nm = kuue.shape[0]
        k_num = nm * nm * self._mesh.n_elem
        k_index = list(np.ndindex(nm, nm))
        k_row = np.zeros(k_num)
        k_col = np.zeros(k_num)
        k_val = np.zeros(k_num)
        m_val = np.zeros(k_num)
        k_ntriplet = 0

        for ni, e in enumerate(self.dof.dof_elements):
            
            for ii, jj in k_index:
                k_row[k_ntriplet] = e.mechanical_dof[ii]
                k_col[k_ntriplet] = e.mechanical_dof[jj]
                k_val[k_ntriplet] = kuue[ii, jj]
                m_val[k_ntriplet] = muue[ii, jj]
                k_ntriplet += 1
            
        muu_shape = (self.dof.n_mdof, self.dof.n_mdof)
        kuu_shape = (self.dof.n_mdof, self.dof.n_mdof)

        self._muu = sparse.coo_matrix((m_val, (k_row, k_col)), shape=muu_shape)
        self._kuu = sparse.coo_matrix((k_val, (k_row, k_col)), shape=kuu_shape)
