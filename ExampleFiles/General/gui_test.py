from PyQt5.QtWidgets import QApplication
from slider_fit_app import SliderFitApp, cPlotAndFit

import numpy as np
import sys, lmfit

class GuiApp(cPlotAndFit):
    def init_data(self):
        self.x = np.linspace(-10, 10, 300)
        self.p = lmfit.Parameters()
        self.p.add("min", -5, min=-10, max=10)
        self.p.add("max", 5,min=-10, max=10)

        self.ymodel = self.get_model(self.p, self.x)

    def get_model(self, p, x):
        return np.arcsin(2*(x-p["min"])/(p["max"]-p["min"]) - 1)
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    aw = SliderFitApp(GuiApp)
    aw.show()
    app.exec_()
