import numpy as np
from .analysis_plate_displacement import PlateDisplacement
from .analysis_laminate_displacement import LaminateDisplacement


class ModeIdentification(object):
    
    def __init__(self, fem, cantilever, type_='laminate'):
        
        x1 = cantilever.xtip - 1
        y0 = cantilever.ytip
        x2 = cantilever.xtip + 1
        
        if type_ == 'plate':
            self._optr1 = PlateDisplacement(fem, (x1, y0)).get_operator()
            self._optr2 = PlateDisplacement(fem, (x2, y0)).get_operator()
        if type_ == 'laminate':
            self._optr1 = LaminateDisplacement(fem, (x1, y0)).get_operator()
            self._optr2 = LaminateDisplacement(fem, (x2, y0)).get_operator()

    
    def is_mode_flexural(self, mode):
        
        disp1 = self._optr1 @ mode
        disp2 = self._optr2 @ mode
        if np.sign(disp1) == np.sign(disp2):
            return True
        return False
