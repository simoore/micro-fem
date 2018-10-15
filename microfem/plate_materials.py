class PlateMaterial(object):
    
    def __init__(self, h, rho, elastic, nu):
        
        self.h = h
        self.rho = rho
        self.elastic = elastic
        self.nu = nu
    
    
class SoiMumpsMaterial(PlateMaterial):
    
    def __init__(self):
        
        si_h = 10e-6
        si_rho = 2330
        si_elastic = 130e9
        si_nu = 0.29
        
        super().__init__(h=si_h, rho=si_rho, elastic=si_elastic, nu=si_nu)
