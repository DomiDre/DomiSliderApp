from SliderApp.slider_fit_app import cPlotAndFit
import PyQt5.QtWidgets as pyqt5widget
import lmfit
import numpy as np
class cPlotAndFitSAXSSANS(cPlotAndFit):
    def __init__(self, parent=None):
        self.modelfile = "sanspol_modelfile.dat"
        self.sldmodelfile = "sanspol_sldfile.dat"
        
        self.data_path_saxs=None
        self.x_saxs = None
        self.y_saxs = None
        self.sy_saxs = None
        
        self.data_path_sans_sa=None
        self.data_path_sans_la=None
        
        self.x_sans_sa = None
        self.x_sans_la = None
        self.y_sans_sa = None
        self.y_sans_la = None
        self.sy_sans_sa = None
        self.sy_sans_la = None
        
        super().__init__(parent)
        
    def get_sld(self, p):
        sys.exit("Define ysaxssld, ynucsld = get_sld(p)")
    
    def get_model(self, p):
        sys.exit("Define y_saxs, y_sans_sa, y_sans_la = get_model(p)")

    def get_dof(self):
        self.dof =  len(self.x_saxs) +\
                    len(self.x_sans_la) +\
                    len(self.x_sans_sa) 
        for param in self.p:
            if self.p[param].vary:
                self.dof -= 1

    def define_plot_canvas(self):
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)

        # Plot data points (if loaded)
        if self.x_saxs is not None and self.y_saxs is not None and\
           self.sy_saxs is not None:
            self.ax1.errorbar(self.x_saxs, self.y_saxs, self.sy_saxs, marker='.',\
                    linestyle='None', color='#a6dba0',\
                    label=self.data_path_saxs, zorder=0)
        
        if self.x_sans_sa is not None and self.y_sans_sa is not None and\
           self.sy_sans_sa is not None:
            self.ax1.errorbar(self.x_sans_sa, self.y_sans_sa, self.sy_sans_sa, marker='.',\
                    linestyle='None', color='#4dac26',\
                    label=self.data_path_sans_sa, zorder=0)

        if self.x_sans_la is not None and self.y_sans_la is not None and\
           self.sy_sans_la is not None:
            self.ax1.errorbar(self.x_sans_la, self.y_sans_la, self.sy_sans_la, marker='.',\
                    linestyle='None', color='#b8e186',\
                    label=self.data_path_sans_la, zorder=0)

        self.model_plot_saxs, = self.ax1.plot(self.x_saxs, self.ymodel_saxs,\
                            marker='None', linestyle='-', color='#008837',\
                            lw=1, label="Model", zorder=1)
        self.model_plot_sans_sa, = self.ax1.plot(self.x_sans_sa, self.ymodel_sans_sa,\
                            marker='None', linestyle='-', color='#0571b0',\
                            lw=1, label="Model", zorder=1)
        self.model_plot_sans_la, = self.ax1.plot(self.x_sans_la, self.ymodel_sans_la,\
                            marker='None', linestyle='-', color='#92c5de',\
                            lw=1, label="Model", zorder=1)
        
        self.ax1.set_xscale('log')
        self.ax1.set_yscale('log')
        self.ax1.set_xlim([min(min(self.x_saxs), min(self.x_sans_sa)),\
                           max(max(self.x_saxs), max(self.x_sans_la))])
        if self.y_sans_la is not None:
            self.ax1.set_ylim([min(min(self.y_saxs), min(self.y_sans_la))*0.8,\
                               max(max(self.y_saxs), max(self.y_sans_sa))*1.2])
        else:
            self.ax1.set_ylim([min(min(self.ymodel_saxs), min(self.ymodel_sans_la))*0.8,\
                               max(max(self.ymodel_saxs), max(self.ymodel_sans_sa))*1.2])
        self.ax1.set_xlabel("$\mathit{q} \, / \, \AA^{-1}$")
        self.ax1.set_ylabel("$\mathit{I} \, / \, cm^{-1}$")
        
        
        self.sldsaxs_plot, = self.ax2.plot(self.xsld/10., self.ysaxssld*1e6,\
                marker='None', ls='-', color='#008837')
        self.sldnuc_plot, = self.ax2.plot(self.xsld/10., self.ynucsld*1e6,\
                marker='None', ls='-', color='#0571b0')

        self.ax2.set_xlim(min(self.xsld/10.), max(self.xsld/10.))
        self.ax2.set_ylim(0, max(max(self.ynucsld), max(self.ysaxssld))*1.1e6)
        self.ax2.set_xlabel("$\mathit{z} \, / \, nm$")
        self.ax2.set_ylabel("$\mathit{SLD} \, / \, 10^{-6} \AA^{-2}$")
        self.fig.subplots_adjust(wspace=0.5)
        
    def update_plot(self):
        self.get_model(self.p)
        
        self.model_plot_saxs.set_ydata(self.ymodel_saxs)
        self.model_plotp_sa.set_ydata(self.ymodel_sans_sa)
        self.model_plotp_la.set_ydata(self.ymodel_sans_la)
        
        self.get_sld(self.p)
        
        self.sldsaxs_plot.set_ydata(self.ysaxssld*1e6)
        self.sldnuc_plot.set_ydata(self.ynucsld*1e6)
        
        if self.y_sans_sa is not None:
            fom = self.figure_of_merit(self.p)
            self.chi2 = sum(fom**2)/self.dof
            self.update_chi2()
        
        self.draw()
        
    def figure_of_merit(self, p):
        self.get_model(p)

        resi_saxs = (np.log(self.ymodel_saxs)-np.log(self.y_saxs))/self.sy_saxs*\
                  self.y_saxs

        resi_sans_sa = (np.log(self.ymodel_sans_sa)-np.log(self.y_sans_sa))/self.sy_sans_sa*\
                  self.y_sans_sa
        resi_sans_la = (np.log(self.ymodel_sans_la)-np.log(self.y_sans_la))/self.sy_sans_la*\
                  self.y_sans_la
        if "weightSAXS" in p:
            resi_saxs = resi_saxs * p["weightSAXS"].value
        resi = np.concatenate([resi_saxs, resi_sans_sa, resi_sans_la])
        
        return resi
        
    def ex_sansort_model(self):
        savefile = open(self.modelfile, "w")
        if self.fit_result is not None:
            savefile.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        savefile.write("\n#q / A-1 \t I / cm-1 \t sI / cm-1"+\
                " \t Imodel / cm-1\n")
        savefile.write("#I saxs\n")
        for iq, qval in enumerate(self.x_saxs):
            savefile.write(str(qval) +"\t"+ str(self.y_saxs[iq])+"\t"+\
                    str(self.sy_saxs[iq]) + "\t" + str(self.ymodel_saxs[iq])+"\n")
        savefile.write("#I sans sa\n")
        for iq, qval in enumerate(self.x_sans_sa):
            savefile.write(str(qval) +"\t"+ str(self.y_sans_sa[iq])+"\t"+\
                    str(self.sy_sans_sa[iq]) + "\t" + str(self.ymodel_sans_sa[iq])+"\n")
        savefile.write("#I sans la\n")
        for iq, qval in enumerate(self.x_sans_la):
            savefile.write(str(qval) +"\t"+ str(self.y_sans_la[iq])+"\t"+\
                    str(self.sy_sans_la[iq]) + "\t" + str(self.ymodel_sans_la[iq])+"\n")
        print("Wrote results to " + self.modelfile)
        savefile.close()
        
        savefile = open(self.sldmodelfile, "w")
        if self.fit_result is not None:
            savefile.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        savefile.write("\n#z / nm \t SLDsaxs / 1e-6 A-2 \t SLDnuc / 1e-6 A-2 \t SLDmag / 1e-6 A-2\n")
        for iq, xval in enumerate(self.xsld):
            savefile.write(str(xval/10.) +\
                            "\t"+ str(self.ysaxssld[iq]*1e6)+\
                            "\t"+ str(self.ynucsld[iq]*1e6)+"\n")
        print("Wrote results to " + self.sldmodelfile)
        savefile.close()
