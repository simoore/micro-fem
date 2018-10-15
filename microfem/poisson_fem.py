import numpy as np
import scipy.sparse as sparse
from .poisson_dof import PoissonDOF
from .poisson_model import PoissonModel
from .mesh import UniformMesh


class PoissonFEM(object):
    """
    Public Attributes
    -----------------
    self.poisson_domain : microfem.PoissonDomain
        The object that describes the shape of the domain, the conductivty of
        the elements, and the energy sources on each element.
    self.dof : microfem.PoissonDOF
        The object provides access to the degrees-of-freedom and mesh 
        parameters.
    """
    def __init__(self, poisson_domain):
        
        mesh = UniformMesh(poisson_domain.domain)
        model = PoissonModel(poisson_domain.a, poisson_domain.b)
        self.poisson_domain = poisson_domain
        self.dof = PoissonDOF(mesh)
        self._ke = model.ke
        self._fe = model.fe
        self._k = mesh.domain2array(poisson_domain.conductivity)
        self._q = mesh.domain2array(poisson_domain.source)
        
        
        # Create row and col vectors for sparce matrices 
        def kr_term(e): return np.hstack([e.dofs for _ in range(4)])
        def kc_term(e): return np.hstack([d*np.ones(4) for d in e.dofs]) 
        kr = np.hstack([kr_term(e) for e in self.dof.dof_elements])
        kc = np.hstack([kc_term(e) for e in self.dof.dof_elements])
        fr = np.hstack([e.dofs for e in self.dof.dof_elements])
        fc = np.zeros(fr.shape)
        
        # Precomputed indexes of the sparse matrices.
        self._k_index = (kr, kc)
        self._f_index = (fr, fc)
        self._k_shape = (self.dof.n_dof, self.dof.n_dof)
        self._f_shape = (self.dof.n_dof, 1) 
    
        # Assemble the conductivity and source matrices.
        self._ktau, self._ftau = self._assemble()
        
    
    def get_conduction_matrix(self, free=False):

        if free is False:
            return self._ktau
        return self._ktau[self.dof.free_dofs, :][:, self.dof.free_dofs]
    
    
    def get_heating_matrix(self, free=False):
        
        if free is False:
            return self._ftau
        return self._ftau[self.dof.free_dofs, :]
    
    
    def solve(self):
        """Applies the boundary conditions onto the system matrices (ktau) 
        and (ftau). Computes the solution of the equation Ku=f. Reinserts the
        boundary DOFs into the solution then returns.
        """
        sysk = self.get_conduction_matrix(free=True)
        sysf = self.get_heating_matrix(free=True)
        ufree = sparse.linalg.spsolve(sysk, sysf)
        uall = np.zeros(self.dof.all_dofs.shape)
        uall[self.dof.free_dofs] = ufree
        return uall, ufree
    
        
    def _assemble(self):
        """
        Determine conductance and heat source for each element from the
        puesdo density (xs). Scales the conductance and heat source element
        matrices (self._ke, self._fe) for each element and assembles them into
        sparse matrices (self.ktau, self.ftau).
        """
        kexpand = np.expand_dims(self._k, axis=1)
        fexpand = np.expand_dims(self._q, axis=1)
        kv = np.ravel(kexpand @ np.expand_dims(self._ke.ravel(), axis=0))
        fv = np.ravel(fexpand @ self._fe.T)

        sysk = sparse.coo_matrix((kv, self._k_index), shape=self._k_shape)
        sysf = sparse.coo_matrix((fv, self._f_index), shape=self._f_shape)
        ktau = sysk.tocsr()
        ftau = sysf.tocsr()
        return ktau, ftau
    