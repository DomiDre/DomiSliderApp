try:
    from .slider_fit_app import cPlotAndFit
except:
    from slider_fit_app import cPlotAndFit
import PyQt5.QtWidgets as pyqt5widget
import lmfit
import numpy as np
class cPlotAndFitRefl(cPlotAndFit):
    def __init__(self, parent=None):
        self.reflsavefile = "gui_refl.dat"
        self.sldsavefile = "gui_sld.dat"
        self.data_path = None
        super().__init__(parent)
        
    def define_plot_canvas(self):
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
        
        if self.x is not None and self.y is not None and self.sy is not None:
            self.ax1.errorbar(self.x, self.y, self.sy, marker='.',\
                    linestyle='None', color='#2c7bb6', label=self.data_path,\
                    zorder=0)

        self.model_plot, = self.ax1.plot(self.x, self.ymodel, marker='None',\
                linestyle='-', color='#ca0020', lw=1, label="Model",\
                zorder=1)
        
        self.ax1.set_yscale('log')
        self.ax1.set_xlim([min(self.x), max(self.x)])
        if self.y is not None:
            self.ax1.set_ylim([min(self.y)*0.8, max(self.y)*1.2])
        else:
            self.ax1.set_ylim([min(self.ymodel)*0.8, max(self.ymodel)*1.2])
        self.ax1.set_xlabel("$\mathit{q_z} \, / \, \AA^{-1}$")
        self.ax1.set_ylabel("$\mathit{I} \, / \, a.u.$")
        
        self.sld_plot, = self.ax2.plot(self.xsld/10., np.real(self.ysld)*1e6,\
                marker='None', ls='-', color='#e66101')
        self.sld_plot_imag, = self.ax2.plot(self.xsld/10., np.imag(self.ysld)*1e6,\
                marker='None', ls='-', color='red')
        self.ax2.set_xlim(min(self.xsld/10.), max(self.xsld/10.))
        self.ax2.set_ylim(0, max(np.real(self.ysld)*1e6)*1.2)
        self.ax2.set_xlabel("$\mathit{z} \, / \, nm$")
        self.ax2.set_ylabel("$\mathit{SLD} \, / \, 10^{-6} \AA^{-2}$")
        self.fig.subplots_adjust(wspace=0.5)
        
    def update_plot(self):
        self.get_model(self.p)
        self.model_plot.set_ydata(self.ymodel)
        self.sld_plot.set_ydata(np.real(self.ysld)*1e6)
        self.sld_plot_imag.set_ydata(np.imag(self.ysld)*1e6)
        
        if self.y is not None:
            fom = self.figure_of_merit(self.p)
            self.chi2 = sum(fom**2)/self.dof
        else:
            fom = 0
            self.chi2 = 0
        self.update_chi2()

        self.draw()
        
    def figure_of_merit(self, p):
        self.get_model(p)
        return (np.log(self.ymodel)-np.log(self.y))/self.sy*self.y
        
    def export_model(self):
        savefile = open(self.reflsavefile, "w")
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
        print("Wrote results to " + self.reflsavefile)
        savefile.close()
        
        savefile = open(self.sldsavefile, "w")
        if self.fit_result is not None:
            savefile.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        savefile.write("\n#z / nm \t real(SLD) / 1e-6 A-2\t imag(SLD) / 1e-6 A-2\n")
        for iq, xval in enumerate(self.xsld):
            savefile.write(str(xval/10.) +"\t"+ str(np.real(self.ysld[iq])*1e6)+\
            "\t"+ str(np.imag(self.ysld[iq])*1e6)+"\n")
        print("Wrote results to " + self.sldsavefile)
        savefile.close()
        
    def save_plot(self):
        self.fig.savefig(self.reflsavefile.rsplit(".",1)[0]+"_plot.png")

