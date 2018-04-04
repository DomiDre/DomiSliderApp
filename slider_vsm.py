try:
    from .slider_fit_app import cPlotAndFit
except:
    from slider_fit_app import cPlotAndFit

import PyQt5.QtWidgets as pyqt5widget
import lmfit
import numpy as np

class cPlotAndFitVSM(cPlotAndFit):
    def __init__(self, parent=None):
        self.modelfile = "vsm_modelfile.dat"
        self.sldmodelfile = "vsm_sldfile.dat"
        super().__init__(parent)
        
    def define_plot_canvas(self):
        self.ax1 = self.fig.add_subplot(111)
        
        if self.x is not None and self.y is not None and self.sy is not None:
            self.ax1.errorbar(self.x, self.y, self.sy, marker='.',\
                    linestyle='None', color='#0571b0', label=self.data_path, zorder=0)

        self.model_plot, = self.ax1.plot(self.x, self.ymodel, marker='None',\
                linestyle='-', color='#ca0020', lw=1, label="Model", zorder=1)
        
        self.ax1.set_xlim([min(self.x), max(self.x)])
        if self.y is not None:
            self.ax1.set_ylim([min(self.y)*1.05, max(self.y)*1.05])
        else:
            self.ax1.set_ylim([min(self.ymodel)*0.95, max(self.ymodel)*1.05])
        self.ax1.set_xlabel("$\mathit{\mu_0 H} \, / \, T$")
        self.ax1.set_ylabel("$\mathit{M} \, / \, kAm^{-1}$")
        
        self.fig.subplots_adjust(wspace=0.5)
        
    def update_plot(self):
        self.ymodel = self.get_model(self.p, self.x)
        self.model_plot.set_ydata(self.ymodel)
        
        if self.y is not None:
            fom = self.figure_of_merit(self.p)
            self.chi2 = sum(fom**2)/self.dof
        else:
            fom = 0
            self.chi2 = 0
        self.update_chi2()
        
        self.draw()
        
    def figure_of_merit(self, p):
        self.ymodel = self.get_model(p, self.x)
        return (self.ymodel-self.y)/self.sy
        
    def export_model(self):
        savefile = open(self.modelfile, "w")
        if self.fit_result is not None:
            savefile.write("#"+\
                       lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        savefile.write("\n#B / T \t M / kAm-1 \t sM / kAm-1"+\
            " \t Mmodel / kAm-1\n")

        for iB, Bval in enumerate(self.x):
            savefile.write(str(Bval) +"\t"+ str(self.y[iB])+"\t"+\
                    str(self.sy[iB]) + "\t" + str(self.ymodel[iB])+"\n")
        print("Wrote results to " + self.modelfile)
        savefile.close()