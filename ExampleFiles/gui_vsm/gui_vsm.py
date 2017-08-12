from PyQt5.QtWidgets import QApplication
from slider_fit_app import SliderFitApp, cPlotAndFit

import numpy as np
import sys, lmfit
from VSMmodels import vsm_magnetization_functions as vsmmodel

class GuiApp(cPlotAndFit):
    def init_data(self):
        self.x = np.linspace(-10, 10, 300)
        self.p = lmfit.Parameters()
        self.p.add("Ms", 200, min=0, max=500)
        self.p.add("xi", 5,min=1e-10, max=100)
        self.p.add("chi", 0,min=-10, max=10)
        self.ymodel = self.get_model(self.p, self.x)

    def get_model(self, p, x):
        return vsmmodel.langevin_function(x, p["Ms"], p["xi"]) + p["chi"]*x
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    aw = SliderFitApp(GuiApp)
    aw.show()
    app.exec_()
