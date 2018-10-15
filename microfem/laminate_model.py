import numpy as np

class LaminateModel(object):
    
    _points = [[-0.577350269189626, -0.577350269189626],
               [0.577350269189626, -0.577350269189626],
               [0.577350269189626, 0.577350269189626],
               [-0.577350269189626, 0.577350269189626]]
    
    def __init__(self, material, a, b):
        """
        Parameters
        ----------
        a : float
            Units are in um.
        b : float
            Units are in um.
        """
        self._a = a * 1e-6
        self._b = b * 1e-6
        self._jacobian = a * b
        self._material = material
        
        elements = self._generate_element_matrices()
        self._muue, self._kuue, self._kuve, self._kvve = elements
        
    
    def get_mass_element(self):
        return self._muue
    
    
    def get_stiffness_element(self):
        return self._kuue
    
    
    def get_piezoelectric_element(self):
        return self._kuve
    
    
    def get_capacitance_element(self):
        return self._kvve
    
        
    def _generate_element_matrices(self):  
        
        cs1, cs2, cs3, ce1, ce2, cc, cm = self._material.get_fem_parameters()
                
        muue = np.zeros((20, 20))
        kuue = np.zeros((20, 20))
        kuve = np.zeros((20, 1))
        
        # Four point integration for bending stiffness and piezoelectric effect.
        be = self._dofs_to_electric_field_matrix()
        for p in self._points:
            bs1, bs2, bs3 = self._dofs_to_strain_matrix(p)
            nu = self._dofs_to_displacement_matrix(p)
            muue += self._jacobian * nu.T @ cm @ nu
            kuue += self._jacobian * (bs1.T @ cs1 @ bs1 + bs2.T @ cs2 @ bs1)
            kuue += self._jacobian * (bs1.T @ cs2 @ bs2 + bs2.T @ cs3 @ bs2)
            kuve += self._jacobian * ((bs1.T + bs3.T) @ ce1 * be)
            kuve += self._jacobian * (bs2.T @ ce2 * be)
        
        
        # Single point integration for shear stiffness or capacitance.
        point, weight = (0, 0), 4
        _, _, bs3 = self._dofs_to_strain_matrix(point)
        kuue += weight * self._jacobian * (bs3.T @ cs1 @ bs3)
        kvve = np.array([[weight * self._jacobian * cc]])
        
        # Enforce symmetry.
        muue = 0.5 * (muue + muue.T)
        kuue = 0.5 * (kuue + kuue.T)
        return muue, kuue, kuve, kvve
    
    
    def _dofs_to_strain_matrix(self, point):
        bs1 = [None for _ in range(4)]
        bs2 = [None for _ in range(4)]
        bs3 = [None for _ in range(4)]
        for i in range(4):
            n, dndx, dndy = self._shapes(point, i)
            bs1[i] = np.array([[dndx, 0, 0, 0, 0], 
                               [0, dndy, 0, 0, 0], 
                               [0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0],
                               [dndy, dndx, 0, 0, 0]])
            bs2[i] = np.array([[0, 0, 0, 0, dndx], 
                               [0, 0, 0, -dndy, 0], 
                               [0, 0, 0, 0, 0],
                               [0, 0, 0, 0, 0],
                               [0, 0, 0, -dndx, dndy]])
            bs3[i] = np.array([[0, 0, 0, 0, 0], 
                               [0, 0, 0, 0, 0], 
                               [0, 0, dndy, -n, 0],
                               [0, 0, dndx, 0, n],
                               [0, 0, 0, 0, 0]])
        bs1 = np.hstack(bs1)
        bs2 = np.hstack(bs2)
        bs3 = np.hstack(bs3)
        return bs1, bs2, bs3
    
    
    def _dofs_to_displacement_matrix(self, point):
        nu = [None for _ in range(4)]
        for i in range(4):
            n, _, _ = self._shapes(point, i)
            nu[i] = np.diag((n, n, n, n, n))
        nu = np.hstack(nu)
        return nu
    
    
    def _dofs_to_electric_field_matrix(self):
        be = 1/self._material.he
        return be
    
    
    def _shapes(self, point, index):
        """The index refers to a node in the normalized element.
        index = 0 : node sw
        index = 1 : node se
        index = 2 : node ne
        index = 3 : node nw
        """
        xi, eta = point
        xi_sign = [-1, 1, 1, -1]
        eta_sign = [-1, -1, 1, 1]
        xs, es = xi_sign[index], eta_sign[index]
        n = 0.25 * (1 + xs * xi) * (1 + es * eta)
        dndx = xs * 0.25 * (1 + es * eta) / self._a
        dndy = es * 0.25 * (1 + xs * xi) / self._b
        return n, dndx, dndy
    