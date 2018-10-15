import numpy as np


class PlateModel(object):
    
    
    _points = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) / np.sqrt(3)

    
    def __init__(self, material, a, b):
        """
        Parameters
        ----------
        a : float
            Units in um.
        b : float 
            Units in um.
        """
        self.ke = self._calculate_ke(a * 1e-6, b * 1e-6, material)
        self.me = self._calculate_me(a * 1e-6, b * 1e-6, material)
         
    
    @staticmethod
    def _shapes(point, a, b):
        
        xi, eta = point
        xi_sign = np.array([-1, 1, 1, -1])
        eta_sign = np.array([-1, -1, 1, 1])
        n = 0.25 * (1 + xi_sign * xi) * (1 + eta_sign * eta)
        dndx = xi_sign * 0.25 * (1 + eta_sign * eta) / a
        dndy = eta_sign * 0.25 * (1 + xi_sign * xi) / b
        return n, dndx, dndy
    
    
    def _calculate_me(self, a, b, material):
        
        rho = material.rho
        h = material.h
        iw = np.array([[1, 0, 0], [0, 0, 0], [0, 0, 0]])
        it = np.array([[0, 0, 0], [0, 1, 0], [0, 0, 1]])
        jacobian = a * b
        
        # Loops execute four point numerical integration.
        mw = np.zeros((12, 12))
        for i in range(4):
            n, _, _ = self._shapes(PlateModel._points[i], a, b)
            n_full = np.hstack([x * iw for x in n])
            mw = mw + (n_full.T @ n_full)
        mw *= jacobian * rho * h
    
        mt = np.zeros((12, 12))
        for i in range(4):
            n, _, _ = self._shapes(PlateModel._points[i], a, b)
            n_full = np.hstack([x * it for x in n])
            mt = mt + (n_full.T @ n_full)
        mt *= jacobian * rho * h * h * h / 12
        
        me = mw + mt
        return me
        

    def _calculate_ke(self, a, b, material):

        e = material.elastic 
        nu = material.nu
        h = material.h
        kappa = np.pi ** 2 / 12     # shear corretion factor
        jacobian = a * b
        g = 0.5 * e / (1 + nu)      # shear modulus
        ci = np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])

        ki = np.zeros((12, 12))
        for i in range(4):
            f = lambda x, y: np.array([[0, 0, -x], [0, y, 0], [0, x, -y]])
            _, dndx, dndy = self._shapes(PlateModel._points[i], a, b)
            bi = np.hstack([f(x, y) for x, y in zip(dndx, dndy)])
            ki = ki + (bi.T @ ci @ bi)
        ki *= jacobian * h * h * h / 12 * e / (1 - nu * nu)
    
        # One point numerical integration.
        f = lambda x, y, z: np.array([[y, 0, x], [z, -x, 0]])
        n, dndx, dndy = self._shapes([0, 0], a, b)
        bo = np.hstack([f(x, y, z) for x, y, z in zip(n, dndx, dndy)])
        ko = jacobian * kappa * h * g * 4 * (bo.T @ bo)
    
        ke = ki + ko
        ke = 0.5 * (ke + ke.T)  # enforce symmetry
        return ke
