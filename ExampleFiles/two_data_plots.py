from PyQt5.QtWidgets import QApplication
from slider_fit_app import SliderFitApp, cPlotAndFit

import numpy as np
import sys, lmfit


def load_data(datafile):
    rawdata = np.loadtxt(datafile)
    x = rawdata[:,0]
    y = rawdata[:,1]
    sy = rawdata[:,2]
    return x, y, sy
    
class GuiApp(cPlotAndFit):
    def init_data(self):
        self.data_path1 = "datafile1.dat"
        self.data_path2 = "datafile2.dat"
        self.data_path3 = "datafile3.dat"
        self.x_1, self.y_1, self.sy_1 = load_data(self.data_path1)
        self.x_2, self.y_2, self.sy_2 = load_data(self.data_path2)
        self.x_3, self.y_3, self.sy_3 = load_data(self.data_path3)
        
        self.p = lmfit.Parameters()
        self.p.add("n", -5, min=-10, max=10)

        self.ymodel = self.get_model(self.p)

    def get_model(self, p):
        ymodel = p["n"]*self.y_1 + (1-p["n"])*self.y_2
        return ymodel

    def figure_of_merit(self, p):
        self.ymodel = self.get_model(p)
        return (self.ymodel-self.y_3)/self.sy_3
            
    def update_plot(self):
        self.ymodel = self.get_model(self.p)
        self.model_plot.set_ydata(self.ymodel)
        self.draw()
        
    def define_plot_canvas(self):
        self.ax1 = self.fig.add_subplot(111)
        
        self.ax1.set_xlabel("$\mathit{x}$")
        self.ax1.set_ylabel("$\mathit{y}$")
        if self.x_1 is not None and self.y_1 is not None and self.sy_1 is not None:
            self.ax1.errorbar(self.x_1, self.y_1, self.sy_1, marker='.',\
                    linestyle='None', color='blue', label=self.data_path1)
        if self.x_2 is not None and self.y_2 is not None and self.sy_2 is not None:
            self.ax1.errorbar(self.x_2, self.y_2, self.sy_2, marker='.',\
                    linestyle='None', color='red', label=self.data_path2)
        if self.x_3 is not None and self.y_3 is not None and self.sy_3 is not None:
            self.ax1.errorbar(self.x_3, self.y_3, self.sy_3, marker='.',\
                    linestyle='None', color='green', label=self.data_path3)

        self.model_plot, = self.ax1.plot(self.x, self.ymodel, marker='None',\
                linestyle='-', color='#ca0020', lw=1, label="Model")
if __name__ == "__main__":
    app = QApplication(sys.argv)
    aw = SliderFitApp(GuiApp)
    aw.show()
    app.exec_()
