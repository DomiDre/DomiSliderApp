from SliderApp.slider_fit_app import cPlotAndFit
import PyQt5.QtWidgets as pyqt5widget
import lmfit
import numpy as np
class cPlotAndFitSAXSSANSPOL(cPlotAndFit):
    def __init__(self, parent=None):
        self.modelfile = "sanspol_modelfile.dat"
        self.sldmodelfile = "sanspol_sldfile.dat"
        
        self.data_path_saxs=None
        self.x_saxs = None
        self.y_saxs = None
        self.sy_saxs = None
        
        self.data_path_sa=None
        self.data_path_la=None
        
        self.xp_sa = None
        self.xp_la = None
        self.xm_sa = None
        self.xm_la = None
        self.yp_sa = None
        self.yp_la = None
        self.ym_sa = None
        self.ym_la = None
        self.syp_sa = None
        self.syp_la = None
        self.sym_sa = None
        self.sym_la = None

        super().__init__(parent)
        
    def get_sld(self, p, x):
        sys.exit("Define ysaxssld, ynucsld, ymagsld=get_sld(p, x) in cPlotAndFitSAXS")
    
    def get_model(self, p, x_saxs, x_sa, x_la):
        sys.exit("Define yp_sa, yp_la, ym_sa, ym_la=get_model(p, xp_sa, "+\
                 "xp_la, xm_sa, xm_la) in cPlotAndFitSANSPOL")

    def get_dof(self):
        self.dof =  len(self.x_saxs) +\
                    len(self.xp_la) + len(self.xm_la) +\
                    len(self.xp_sa) + len(self.xm_sa) 
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
                    linestyle='None', color='#a6dba0', label=self.data_path_saxs)
        
        if self.xp_sa is not None and self.yp_sa is not None and\
           self.syp_sa is not None:
            self.ax1.errorbar(self.xp_sa, self.yp_sa, self.syp_sa, marker='.',\
                    linestyle='None', color='#4dac26', label=self.data_pathp_sa)

        if self.xp_la is not None and self.yp_la is not None and\
           self.syp_la is not None:
            self.ax1.errorbar(self.xp_la, self.yp_la, self.syp_la, marker='.',\
                    linestyle='None', color='#b8e186', label=self.data_pathp_la)

        if self.xm_sa is not None and self.ym_sa is not None and\
           self.sym_sa is not None:
            self.ax1.errorbar(self.xm_sa, self.ym_sa, self.sym_sa, marker='.',\
                    linestyle='None', color='#d01c8b', label=self.data_pathm_sa)

        if self.xm_la is not None and self.ym_la is not None and\
           self.sym_la is not None:
            self.ax1.errorbar(self.xm_la, self.ym_la, self.sym_la, marker='.',\
                    linestyle='None', color='#f1b6da', label=self.data_pathm_la)

        self.model_plot_saxs, = self.ax1.plot(self.x_saxs, self.ymodel_saxs,\
                            marker='None', linestyle='-', color='#008837',\
                            lw=1, label="Model")
        self.model_plotp_sa, = self.ax1.plot(self.xp_sa, self.ymodelp_sa,\
                            marker='None', linestyle='-', color='#0571b0',\
                            lw=1, label="Model")
        self.model_plotp_la, = self.ax1.plot(self.xp_la, self.ymodelp_la,\
                            marker='None', linestyle='-', color='#92c5de',\
                            lw=1, label="Model")
        self.model_plotm_sa, = self.ax1.plot(self.xm_sa, self.ymodelm_sa,\
                            marker='None', linestyle='-', color='#ca0020',\
                            lw=1, label="Model")
        self.model_plotm_la, = self.ax1.plot(self.xm_la, self.ymodelm_la,\
                            marker='None', linestyle='-', color='#f4a582',\
                            lw=1, label="Model")
        
        self.ax1.set_xscale('log')
        self.ax1.set_yscale('log')
        self.ax1.set_xlim([min(min(self.x_saxs), min(self.xp_sa), min(self.xm_sa)),\
                           max(max(self.x_saxs), max(self.xp_la), max(self.xm_la))])
        if self.yp_la is not None:
            self.ax1.set_ylim([min(min(self.y_saxs), min(self.yp_la), min(self.ym_la))*0.8,\
                               max(max(self.y_saxs), max(self.yp_sa), max(self.ym_sa))*1.2])
        else:
            self.ax1.set_ylim([min(min(self.ymodel_saxs), min(self.ymodelp_la), min(self.ymodelm_la))*0.8,\
                               max(max(self.ymodel_saxs), max(self.ymodelp_sa), max(self.ymodelm_sa))*1.2])
        self.ax1.set_xlabel("$\mathit{q} \, / \, \AA^{-1}$")
        self.ax1.set_ylabel("$\mathit{I} \, / \, cm^{-1}$")
        
        
        self.sldsaxs_plot, = self.ax2.plot(self.xsld/10., self.ysaxssld*1e6,\
                marker='None', ls='-', color='#008837')
        self.sldnuc_plot, = self.ax2.plot(self.xsld/10., self.ynucsld*1e6,\
                marker='None', ls='-', color='#0571b0')
        self.sldmag_plot, = self.ax2.plot(self.xsld/10., self.ymagsld*1e6,\
                marker='None', ls='-', color='#ca0020')

        self.ax2.set_xlim(min(self.xsld/10.), max(self.xsld/10.))
        self.ax2.set_ylim(0, max(self.ynucsld)*1.5e6)
        self.ax2.set_xlabel("$\mathit{z} \, / \, nm$")
        self.ax2.set_ylabel("$\mathit{SLD} \, / \, 10^{-6} \AA^{-2}$")
        self.fig.subplots_adjust(wspace=0.5)
        
    def update_plot(self):
        self.ymodel_saxs,\
            self.ymodelp_sa, self.ymodelp_la,\
            self.ymodelm_sa, self.ymodelm_la =\
                self.get_model(self.p,\
                                self.x_saxs,\
                                self.xp_sa, self.xp_la,\
                                self.xm_sa, self.xm_la)
        
        self.model_plot_saxs.set_ydata(self.ymodel_saxs)
        self.model_plotp_sa.set_ydata(self.ymodelp_sa)
        self.model_plotp_la.set_ydata(self.ymodelp_la)
        self.model_plotm_sa.set_ydata(self.ymodelm_sa)
        self.model_plotm_la.set_ydata(self.ymodelm_la)
        
        self.ysaxssld, self.ynucsld, self.ymagsld = self.get_sld(self.p, self.xsld)
        
        self.sldsaxs_plot.set_ydata(self.ysaxssld*1e6)
        self.sldnuc_plot.set_ydata(self.ynucsld*1e6)
        self.sldmag_plot.set_ydata(self.ymagsld*1e6)
        
        if self.yp_sa is not None:
            fom = self.figure_of_merit(self.p)
            self.chi2 = sum(fom**2)/self.dof
            self.update_chi2()
        
        self.draw()
        
    def figure_of_merit(self, p):
        self.ymodel_saxs,\
            self.ymodelp_sa, self.ymodelp_la,\
            self.ymodelm_sa, self.ymodelm_la =\
                self.get_model(self.p,\
                                self.x_saxs,\
                                self.xp_sa, self.xp_la,\
                                self.xm_sa, self.xm_la)
        resi_saxs = (np.log(self.ymodel_saxs)-np.log(self.y_saxs))/self.sy_saxs*\
                  self.y_saxs
        resip_sa = (np.log(self.ymodelp_sa)-np.log(self.yp_sa))/self.syp_sa*\
                  self.yp_sa
        resip_la = (np.log(self.ymodelp_la)-np.log(self.yp_la))/self.syp_la*\
                  self.yp_la
        resim_sa = (np.log(self.ymodelm_sa)-np.log(self.ym_sa))/self.sym_sa*\
                  self.ym_sa
        resim_la = (np.log(self.ymodelm_la)-np.log(self.ym_la))/self.sym_la*\
                  self.ym_la
        resi = np.concatenate([resi_saxs, resip_sa, resip_la, resim_sa, resim_la])
        
        return resi
        
    def export_model(self):
        savefile = open(self.modelfile, "w")
        if self.fit_result is not None:
            savefile.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        savefile.write("\n#q / A-1 \t I / cm-1 \t sI / cm-1"+\
                " \t Imodel / cm-1\n")
        savefile.write("#I saxs\n")
        for iq, qval in enumerate(self.x_saxs):
            savefile.write(str(qval) +"\t"+ str(self.y_saxs[iq])+"\t"+\
                    str(self.sy_saxs[iq]) + "\t" + str(self.ymodel_saxs[iq])+"\n")
        savefile.write("#I+ sa\n")
        for iq, qval in enumerate(self.xp_sa):
            savefile.write(str(qval) +"\t"+ str(self.yp_sa[iq])+"\t"+\
                    str(self.syp_sa[iq]) + "\t" + str(self.ymodelp_sa[iq])+"\n")
        savefile.write("#I+ la\n")
        for iq, qval in enumerate(self.xp_la):
            savefile.write(str(qval) +"\t"+ str(self.yp_la[iq])+"\t"+\
                    str(self.syp_la[iq]) + "\t" + str(self.ymodelp_la[iq])+"\n")
        savefile.write("#I- sa\n")
        for iq, qval in enumerate(self.xm_sa):
            savefile.write(str(qval) +"\t"+ str(self.ym_sa[iq])+"\t"+\
                    str(self.sym_sa[iq]) + "\t" + str(self.ymodelm_sa[iq])+"\n")
        savefile.write("#I- la\n")
        for iq, qval in enumerate(self.xm_la):
            savefile.write(str(qval) +"\t"+ str(self.ym_la[iq])+"\t"+\
                    str(self.sym_la[iq]) + "\t" + str(self.ymodelm_la[iq])+"\n")
        print("Wrote results to " + self.modelfile)
        savefile.close()
        
        savefile = open(self.sldmodelfile, "w")
        if self.fit_result is not None:
            savefile.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        savefile.write("\n#z / nm \t SLDsaxs / 1e-6 A-2 \t SLDnuc / 1e-6 A-2 \t SLDmag / 1e-6 A-2\n")
        for iq, xval in enumerate(self.xsld):
            savefile.write(str(xval/10.) +\
                            "\t"+ str(self.ysaxssld[iq]*1e6)+\
                            "\t"+ str(self.ynucsld[iq]*1e6)+\
                            "\t"+ str(self.ymagsld[iq]*1e6)+"\n")
        print("Wrote results to " + self.sldmodelfile)
        savefile.close()
