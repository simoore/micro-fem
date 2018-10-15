import numpy as np
from .analysis_plate_displacement import PlateDisplacement


class ModeIdentification(object):
    
    def __init__(self, fem, cantilever):
        
        x1 = cantilever.xtip - 1
        y0 = cantilever.ytip
        x2 = cantilever.xtip + 1
        self._optr1 = PlateDisplacement(fem, (x1, y0)).get_operator()
        self._optr2 = PlateDisplacement(fem, (x2, y0)).get_operator()

    
    def is_mode_flexural(self, mode):
        
        disp1 = self._optr1 @ mode
        disp2 = self._optr2 @ mode
        if np.sign(disp1) == np.sign(disp2):
            return True
        return False
