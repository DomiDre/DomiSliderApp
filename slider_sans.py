from slider_fit_app import cPlotAndFit
import PyQt5.QtWidgets as pyqt5widget
import lmfit
import numpy as np
class cPlotAndFitSANS(cPlotAndFit):
    def __init__(self, parent=None):
        self.modelfile = "sans_modelfile.dat"
        self.sldmodelfile = "sans_sld_modelfile.dat"
        self.data_path_sa=None
        self.data_path_la=None
        super().__init__(parent)
        
    def get_sld(self, p, x):
        sys.exit("Define sld=get_sld(p, x) in cPlotAndFitSAXS")
    
    def get_model(self, p, x_sa, x_la):
        sys.exit("Define y_sa, y_la=get_model(p, x_sa, x_la) in cPlotAndFitSANS")
        
    def define_plot_canvas(self):
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
        
        if self.x_sa is not None and self.y_sa is not None and self.sy_sa is not None:
            self.ax1.errorbar(self.x_sa, self.y_sa, self.sy_sa, marker='.',\
                    linestyle='None', color='#4dac26', label=self.data_path_sa)

        if self.x_la is not None and self.y_la is not None and self.sy_la is not None:
            self.ax1.errorbar(self.x_la, self.y_la, self.sy_la, marker='.',\
                    linestyle='None', color='#b8e186', label=self.data_path_la)

        self.model_plot_sa, = self.ax1.plot(self.x_sa, self.ymodel_sa, marker='None',\
                linestyle='-', color='#ca0020', lw=1, label="Model")
        self.model_plot_la, = self.ax1.plot(self.x_la, self.ymodel_la, marker='None',\
                linestyle='-', color='#f4a582', lw=1, label="Model")
        
        self.ax1.set_xscale('log')
        self.ax1.set_yscale('log')
        self.ax1.set_xlim([min(self.x_sa), max(self.x_la)])
        self.ax1.set_ylim([min(self.y_la[self.y_la>0])*0.8, max(self.y_sa)*1.2])
        self.ax1.set_xlabel("$\mathit{q_z} \, / \, \AA^{-1}$")
        self.ax1.set_ylabel("$\mathit{I} \, / \, cm^{-1}$")
        
        
        self.sld_plot, = self.ax2.plot(self.xsld/10., self.ysld*1e6,\
                marker='None', ls='-', color='#e66101')
        self.ax2.set_xlim(min(self.xsld/10.), max(self.xsld/10.))
        self.ax2.set_ylim(0, max(self.ysld)*1.5e6)
        self.ax2.set_xlabel("$\mathit{z} \, / \, nm$")
        self.ax2.set_ylabel("$\mathit{SLD} \, / \, 10^{-6} \AA^{-2}$")
        self.fig.subplots_adjust(wspace=0.5)
        
    def update_plot(self):
        self.ymodel_sa, self.ymodel_la = self.get_model(self.p, \
                                                        self.x_sa, self.x_la)
        self.model_plot_sa.set_ydata(self.ymodel_sa)
        self.model_plot_la.set_ydata(self.ymodel_la)
        
        self.ysld = self.get_sld(self.p, self.xsld)
        self.sld_plot.set_ydata(self.ysld*1e6)
        self.update_chi2()
        self.draw()
        
    def figure_of_merit(self, p):
        self.ymodel_sa, self.ymodel_la = self.get_model(p, \
                                                        self.x_sa, self.x_la)
        resi_sa = (np.log(self.ymodel_sa)-np.log(self.y_sa))/self.sy_sa*\
                  self.y_sa
        resi_la = (np.log(self.ymodel_la)-np.log(self.y_la))/self.sy_la*\
                  self.y_la
        resi = np.concatenate([resi_sa, resi_la])
        return resi
        
    def export_model(self):
        savefile = open(self.modelfile, "w")
        if self.fit_result is not None:
            savefile.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        savefile.write("#q / A-1 \t I / cm-1 \t sI / cm-1"+\
                " \t Imodel / cm-1\n")
        for iq, qval in enumerate(self.x_sa):
            savefile.write(str(qval) +"\t"+ str(self.y_sa[iq])+"\t"+\
                    str(self.sy_sa[iq]) + "\t" + str(self.ymodel_sa[iq])+"\n")
        savefile.write("#\n")
        for iq, qval in enumerate(self.x_la):
            savefile.write(str(qval) +"\t"+ str(self.y_la[iq])+"\t"+\
                    str(self.sy_la[iq]) + "\t" + str(self.ymodel_la[iq])+"\n")
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
