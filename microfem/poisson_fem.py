import numpy as np
import scipy.sparse as sparse
from poisson_dof import PoissonDOF
from poisson_model import PoissonModel


class PoissonFEM(object):
    
    def __init__(self, mesh):
        
        model = PoissonModel(mesh.a, mesh.b)
        
        self._dof = PoissonDOF(mesh)
        self._ke = model.ke
        self._fe = model.fe
        
        # Create row and col vectors for sparce matrices 
        def kr_term(e): return np.hstack([e.dofs for _ in range(4)])
        def kc_term(e): return np.hstack([d*np.ones(4) for d in e.dofs]) 
        kr = np.hstack([kr_term(e) for e in self._dof.dof_elements])
        kc = np.hstack([kc_term(e) for e in self._dof.dof_elements])
        fr = np.hstack([e.dofs for e in self._dof.dof_elements])
        fc = np.zeros(fr.shape)
        
        self._k_index = (kr, kc)
        self._f_index = (fr, fc)
        self._k_shape = (self._dof.n_dof, self._dof.n_dof)
        self._f_shape = (self._dof.n_dof, 1) 
        self.assemble(mesh.get_densities())
        
        
    @property
    def dof(self):
        return self._dof
    
    
    def get_conduction_matrix(self, free=False):

        if free is False:
            return self._ktau
        return self._ktau[self._dof.free_dofs, :][:, self._dof.free_dofs]
    
    
    def get_heating_matrix(self, free=False):
        
        if free is False:
            return self._ftau
        return self._ftau[self._dof.free_dofs, :]
    
    
    def thermal_analysis(self):
        """Applies the boundary conditions onto the system matrices (ktau) 
        and (ftau). Computes the solution of the equation Ku=f. Reinserts the
        boundary DOFs into the solution then returns.
        """
        sysk = self.get_conduction_matrix(free=True)
        sysf = self.get_heating_matrix(free=True)
        ufree = sparse.linalg.spsolve(sysk, sysf)
        uall = np.zeros(self._dof.all_dofs.shape)
        uall[self._dof.free_dofs] = ufree
        return uall, ufree
    
        
    def assemble(self, xs):
        """
        Determine conductance and heat source for each element from the
        puesdo density (xs). Scales the conductance and heat source element
        matrices (self._ke, self._fe) for each element and assembles them into
        sparse matrices (self.ktau, self.ftau).
        """
        if xs is None:
            return
        
        k, q = self._set_penalty(xs)
        kexpand = np.expand_dims(k, axis=1)
        fexpand = np.expand_dims(q, axis=1)
        kv = np.ravel(kexpand @ np.expand_dims(self._ke.ravel(), axis=0))
        fv = np.ravel(fexpand @ self._fe.T)
        
        sysk = sparse.coo_matrix((kv, self._k_index), shape=self._k_shape)
        sysf = sparse.coo_matrix((fv, self._f_index), shape=self._f_shape)
        self._ktau = sysk.tocsr()
        self._ftau = sysf.tocsr()
    
    
    def _set_penalty(self, xs):
        
        q0, k0, eps = 1, 1e4, 1e-4
        k = k0*xs + (1 - xs)*eps*k0
        q = q0*xs
        self.thermal_grad = (k0 - eps*k0)*np.ones_like(xs)
        self.heat_grad = q0*np.ones_like(xs)
        return k, q
