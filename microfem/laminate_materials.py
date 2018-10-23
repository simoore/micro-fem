import numpy as np


class LaminateMaterial(object):
    """
    Attributes
    ----------
    self.he : float
    The thickness of the piezoelectric layer in m.
    """
    
    def __init__(self, layers, piezo_layer):

        # Set z coordinate of layers.
        thickness = sum((l.h for l in layers))
        layers[0].zl = -thickness / 2
        layers[0].zu = layers[0].zl + layers[0].h
        for i in range(1, len(layers)):
            layers[i].zl = layers[i - 1].zu
            layers[i].zu = layers[i].zl + layers[i].h

        # Piezoelectric layer properties.
        self.he = piezo_layer.h
        self._ce1 = piezo_layer.e.T * (piezo_layer.zu - piezo_layer.zl)
        self._ce2 = piezo_layer.e.T * (piezo_layer.zu**2 - piezo_layer.zl**2)/2
        self._cc = piezo_layer.epsilon / piezo_layer.h

        # Laminate mass properties.
        self._cm0, self._cm1, self._cm2 = self._generate_cm(layers)

        # Set elastic matrices.
        self._cs1, self._cs2, self._cs3 = self._generate_cs(layers)


    @staticmethod
    def _generate_cs(layers):
        cs1 = np.zeros((5, 5))
        cs2 = np.zeros((5, 5))
        cs3 = np.zeros((5, 5))
        for l in layers:
            cs1 += l.c * (l.zu - l.zl)
            cs2 += l.c * (l.zu ** 2 - l.zl ** 2) / 2
            cs3 += l.c * (l.zu ** 3 - l.zl ** 3) / 3
        return cs1, cs2, cs3


    @staticmethod
    def _generate_cm(layers):
        p0, p1, p2 = 0, 0, 0
        for l in layers:
            p0 += l.rho * (l.zu - l.zl)
            p1 += l.rho * (l.zu ** 2 - l.zl ** 2) / 2
            p2 += l.rho * (l.zu ** 3 - l.zl ** 3) / 3
        return p0, p1, p2
    

    def get_fem_parameters(self):
        return (self._cs1, self._cs2, self._cs3, self._ce1, self._ce2,
                self._cc, self._cm0, self._cm1, self._cm2)


class LaminateLayer(object):
    """
    Attributes
    ----------
    self.zl : float
    The z-coordinate of the bottom of the layer.
    
    self.zh : float
    The z-coordinate of the top of the layer.
    
    self.h : float
    The thickness of the layer in m.
    
    self.rho : float
    The density of the layer in kg/m^3.
    
    self.c : ndarray
    The elastic modulus of the layer.
    
    self.e : ndarray
    The piezoelectric constants of the layer.
    
    self.epsilon : float
    The permitivity of the layer.
    """
    
    
    def __init__(self, h, rho, elastic, nu, epsilon=0.0, e_piezo=0.0):
        """
        Parameters
        ----------
        e_piezo : float
        This is the strain piezoelectric coefficient e=e_{31}=e_{32} in C/m^2.
        """
        
        # Layer coordinates, thickness, and density.
        self.zl = 0.0
        self.zu = 0.0
        self.h = h
        self.rho = rho

        # Elastic, piezoelectricv and capacitance coefficents.
        self.c = self._generate_elastic_matrix(elastic, nu)
        self.e = self._generate_piezoelectric_matrix(e_piezo)
        self.epsilon = epsilon
        
    
    @staticmethod
    def _generate_elastic_matrix(elastic, nu):
        kappa = np.pi ** 2 / 12
        shear_modulus = 0.5 * elastic / (1 + nu)
        c11 = elastic / (1 - nu ** 2)
        c12 = c11 * nu
        c44 = kappa * shear_modulus
        c66 = elastic / (2 * (1 + nu))
        ce = np.array([[c11, c12, 0, 0, 0],
                       [c12, c11, 0, 0, 0],
                       [0, 0, c44, 0, 0],
                       [0, 0, 0, c44, 0],
                       [0, 0, 0, 0, c66]])
        return ce
   
    
    @staticmethod
    def _generate_piezoelectric_matrix(e_piezo):
        e = np.array([[e_piezo, e_piezo, 0, 0, 0]])
        return e
    

class PiezoMumpsMaterial(LaminateMaterial):
    """PiezoMUMPS is a micofabrication process provided by the company 
    MEMSCAP. The devices produced using this process consist of three 
    primary layers. A 10um layer of doped single crystal silicon, a 0.5um
    layer of aluminium nitride, and a 1um layer of aluminium.
    
    The silicon material properties are from:
        
    Matthew A. Hopcroft, William D. Nix, and Thomas W. Kenny; What is the 
    Young's Modulus of Silicon?; Journal of Microelectromechanical systems 
    19(2) 2010 pp 229-238.
    
    Material properties of aluminium nitride are published in the following 
    article. The material properties are expressed as an orthotropic material
    and are converted to isotropic plate properties.
    
    Kazuo Tsubouch, and Nobuo Mikoshiba; Zero-Temperature-Coefficient SAW 
    Devices on AlN Epitaxial Films; IEEE Transactions on Sonics and Ultrasonics 
    SU-32(5) 1985 pp 634-644.
    
    This article is referenced anymore.
    
    Joseph C Doll and Beth L Pruitt; Design of piezoresistive versus 
    piezoelectric contact mode scanning probes; J. Micromech. Microeng.
    20 (2010) 095023
    
    Need reference for all poisson's ratios. Need reference for aluminium 
    material properties.
    """

    def __init__(self):
        
        perm_free = 8.85418782e-12

        si_h = 10e-6
        si_rho = 2330
        si_elastic = 130e9
        si_nu = 0.28
        
        si = LaminateLayer(si_h, si_rho, si_elastic, si_nu)

        aln_h = 0.5e-6
        aln_rho = 3260
        aln_e = 300e9
        aln_nu = 0.36
        aln_perm = perm_free * 10.2
        aln_piezo = 0.58 # This is e_31=e_32.
        aln = LaminateLayer(aln_h, aln_rho, aln_e, aln_nu, aln_perm, aln_piezo)

        al_h = 1e-6
        al_rho = 2700
        al_elastic = 70e9
        al_nu = 0.33
        
        al = LaminateLayer(al_h, al_rho, al_elastic, al_nu)

        layers = (si, aln, al)
        super().__init__(layers, aln)
        
        
class SoiMumpsModel(LaminateMaterial):
    
    def __init__(self):
        perm_free = 8.85418782e-12
        
        si_h = 10e-6
        si_rho = 2330
        si_elastic = 169e9
        si_nu = 0.064
        # si_nu = 0.29
        si_perm = perm_free * 0
        si_piezo = 0
        
        si = LaminateLayer(si_h, si_rho, si_elastic, si_nu, si_perm, si_piezo)
        layers = (si,)
        super().__init__(layers, si)
        