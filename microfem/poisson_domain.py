class PoissonDomain(object):
    """Cantilevers are defined by three parameters. (topology) a binary matrix
    indicating the existence of elements in a rectangular mesh, (a) half the 
    width of the element, (b) have the length of the element. The first index
    is the x-coordinate and the second index is the y-coordinate. Elements
    along the x-axis (y == 0) are on the boundary.
    """
    
    def __init__(self, domain, conductivity, source, a, b):
        """
        Parameters
        ----------
        a : float
            Half the element width in x-direction. The units are um.
        b : float
            Half the element width in y-direction. The units are um.
        """
        
        self.domain = domain
        self.conductivity = conductivity
        self.source = source
        self.a = a
        self.b = b
        

    def to_console(self):
        
        nelx, nely = self.domain.shape
        dimensions = (2 * self.a, 2 * self.b)
        mesh = (2 * self.a * nelx, 2 * self.b * nely)
        name = 'Poisson Equation Domain'
        
        print(''.join(('--- ', name, ' ---\n')))
        print('Each element is %g x %g um' % dimensions)
        print('The design area is %g x %g um' % mesh)
      