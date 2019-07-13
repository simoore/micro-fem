class Cantilever(object):
    
    def __init__(self, topology, a, b, xtip, ytip):
        """
        Parameters
        ----------
        topology : binary rank 2 ndarray
            A binary matrix indicating the existence of elements in a 
            rectangular mesh.
        a : float
            Half the element width in x-direction. The units are um.
        b : float
            Half the element width in y-direction. The units are um.
        xtip : float
            The x-coodinate of the tip in the design space. The units are um.
        ytip : float
            The y-coordiate of the tip in the design space. The units are um.
        """
        
        nelx, nely = topology.shape
        self.topology = topology
        self.a = a
        self.b = b
        self.xtip = xtip
        self.ytip = ytip


    def to_console(self):
        
        nelx, nely = self.topology.shape
        dimensions = (2 * self.a, 2 * self.b)
        mesh = (2 * self.a * nelx, 2 * self.b * nely)
        tip_location = (self.xtip, self.ytip)
        name = 'MicroFEM Cantilever'
        
        print(''.join(('--- ', name, ' ---\n')))
        print('Each element is %g x %g um' % dimensions)
        print('The design area is %g x %g um' % mesh)
        print('(xtip, ytip) = (%g, %g) um\n' % tip_location)
      