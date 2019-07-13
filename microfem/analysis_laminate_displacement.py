import numpy as np
import scipy.sparse as sparse


class LaminateDisplacement(object):

    def __init__(self, fem, coords):

        self._elements = fem.dof.dof_elements
        self._xtip = coords[0]
        self._ytip = coords[1]
        self._a = fem.a  # distance in um
        self._b = fem.b  # distance in um
        self._dimension = fem.dof.n_mdof
    
        # Initial values of the operator.
        self._assemble(self._operator_element)

        
    def get_operator(self):
        
        return self._opr
    
    
    def _assemble(self, element_func):
        
        num = 20 * len(self._elements)
        row = np.zeros(num)
        col = np.zeros(num) 
        val = np.zeros(num)
        ntriplet = 0
        
        for e in self._elements:
            dof = e.mechanical_dof
            ge = element_func(e)
            if ge is not None:
                for ii in range(20):
                    row[ntriplet] = 0
                    col[ntriplet] = dof[ii]
                    val[ntriplet] = ge[0, ii]
                    ntriplet += 1

        shape = (1, self._dimension)
        self._opr = sparse.coo_matrix((val, (row, col)), shape=shape).tocsr()
    
    
    def _operator_element(self, element):
        
        x0 = element.element.i + 0.5
        y0 = element.element.j + 0.5
        xi = (self._xtip / self._a - 2 * x0) 
        eta = (self._ytip / self._b - 2 * y0) 
        
        #if (x0 == 29.5 or x0 == 30.5) and (y0 == 59.5):
        #    print(xi, eta)
        #    print(self._xtip, self._ytip)
        #    print(self._xtip - 2 * self._a * x0)
        
        if -1 < xi <= 1 and -1 < eta <= 1:
            xi_sign = np.array([-1, 1, 1, -1])
            eta_sign = np.array([-1, -1, 1, 1])
            n = 0.25 * (1 + xi_sign * xi) * (1 + eta_sign * eta)
            ge = np.array([[0, 0, n[0], 0, 0, 0, 0, n[1], 0, 0, 0, 0, n[2], 0, 0, 0, 0, n[3], 0, 0]])
            return ge
        
        return None
