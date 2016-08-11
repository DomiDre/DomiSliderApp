from slider_fit_app import cPlotAndFit
import PyQt5.QtWidgets as pyqt5widget
import lmfit
import numpy as np
class cPlotAndFitRefl(cPlotAndFit):
    def __init__(self, parent=None):
        self.modelfile = "refl_modelfile.dat"
        self.sldmodelfile = "sld_modelfile.dat"
        self.data_path = None
        super().__init__(parent)
        
    def define_plot_canvas(self):
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
        
        if self.x is not None and self.y is not None and self.sy is not None:
            self.ax1.errorbar(self.x, self.y, self.sy, marker='.',\
                    linestyle='None', color='#4dac26', label=self.data_path)

        self.model_plot, = self.ax1.plot(self.x, self.ymodel, marker='None',\
                linestyle='-', color='#ca0020', lw=1, label="Model")
        
        self.ax1.set_yscale('log')
        self.ax1.set_xlim([min(self.x), max(self.x)])
        if self.y is not None:
            self.ax1.set_ylim([min(self.y)*0.8, max(self.y)*1.2])
        else:
            self.ax1.set_ylim([min(self.ymodel)*0.8, max(self.ymodel)*1.2])
        self.ax1.set_xlabel("$\mathit{q_z} \, / \, \AA^{-1}$")
        self.ax1.set_ylabel("$\mathit{I} \, / \, a.u.$")
        
        
        self.sld_plot, = self.ax2.plot(self.xsld/10., self.ysld*1e6,\
                marker='None', ls='-', color='#e66101')
        self.ax2.set_xlim(min(self.xsld/10.), max(self.xsld/10.))
        self.ax2.set_xlabel("$\mathit{z} \, / \, nm$")
        self.ax2.set_ylabel("$\mathit{SLD} \, / \, 10^{-6} \AA^{-2}$")
        self.fig.subplots_adjust(wspace=0.5)
        
    def update_plot(self):
        self.get_model(self.p)

        self.model_plot.set_ydata(self.ymodel)
        self.sld_plot.set_ydata(self.ysld*1e6)
        self.draw()
        
    def figure_of_merit(self, p):
        self.get_model(p)
        return (np.log(self.ymodel)-np.log(self.y))/self.sy*self.y
        
    def export_model(self):
        savefile = open(self.modelfile, "w")
        if self.fit_result is not None:
            savefile.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        if self.y is not None:
            savefile.write("#q / A-1 \t I / a.u. \t sI / a.u."+\
                " \t Imodel / a.u.\n")
            for iq, qval in enumerate(self.x):
                savefile.write(str(qval) +"\t"+ str(self.y[iq])+"\t"+\
                        str(self.sy[iq]) + "\t" + str(self.ymodel[iq])+"\n")
        else:
            savefile.write("#q / A-1 \t Imodel / a.u.\n")
            for iq, qval in enumerate(self.x):
                savefile.write(str(qval) + "\t" + str(self.ymodel[iq])+"\n")
        print("Wrote results to " + self.modelfile)
        savefile.close()
        
        savefile = open(self.sldmodelfile, "w")
        if self.fit_result is not None:
            savefile.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        savefile.write("#z / nm \t SLD / 1e-6 A-2\n")
        for iq, xval in enumerate(self.xsld):
            savefile.write(str(xval/10.) +"\t"+ str(self.ysld[iq]*1e6)+"\n")
        print("Wrote results to " + self.sldmodelfile)
        savefile.close()
