class Cantilever(object):
    """Cantilevers are defined by three parameters. (topology) a binary matrix
    indicating the existence of elements in a rectangular mesh, (a) half the 
    width of the element, (b) have the length of the element. The first index
    is the x-coordinate and the second index is the y-coordinate. Elements
    along the x-axis (y == 0) are on the boundary.
    """
    
    def __init__(self, topology, a, b, xtip, ytip):
        """
        Parameters
        ----------
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
        #self.xtip = 1e6 * (a * nelx)
        #self.ytip = 1e6 * (2 * nely * b - b)
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
      