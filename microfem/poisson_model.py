import numpy as np

class PoissonModel(object):
    
    def __init__(self, a, b):
        """
        Parameters
        ----------
        a : float 
            Half the width of an element in the x-direction. Units in um.
        b : float
            Half the width of an element in the y-direction. Units in um.
        """
        a = a * 1e-6
        b = b * 1e-6
        points = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) / np.sqrt(3)
        self._ke = np.zeros((4, 4))
        self._fe = np.zeros((4, 1))
        for p in points:
            n, dndxi, dndeta = self.shapes(p)
            self._ke += a*b * (dndxi.T @ dndxi / a*a + dndeta.T @ dndeta / b*b)
            self._fe += a*b * n.T
            
    
    @property
    def ke(self):
        return self._ke
    
    
    @property
    def fe(self):
        return self._fe
        
          
    @staticmethod
    def shapes(point):
        """
        Parameters
        ----------
        point : tuple 
            The point (xi, eta) to evaluate the shape functions.
        """
        xi, eta = point
        xs = np.array([[-1, 1, 1, -1]])
        es = np.array([[-1, -1, 1, 1]])
        n = 0.25 * (1 + xs * xi) * (1 + es * eta)
        dndxi = xs * 0.25 * (1 + es * eta)
        dndeta = es * 0.25 * (1 + xs * xi) 
        return n, dndxi, dndeta