try:
    from .slider_fit_app import cPlotAndFit
except:
    from slider_fit_app import cPlotAndFit

import PyQt5.QtWidgets as pyqt5widget
import lmfit, datetime
import numpy as np
class cPlotAndFitXRPNRefl(cPlotAndFit):
    def __init__(self, parent=None):
        self.reflsavefile = "gui_refl.dat"
        self.sldsavefile = "gui_sld.dat"
        self.data_path_xrr = None
        self.data_path_pnr_p = None
        self.data_path_pnr_m = None

        self.x_xrr = None
        self.y_xrr = None
        self.sy_xrr = None
        self.ymodel_xrr = None
        
        self.x_p = None
        self.y_p = None
        self.sy_p = None
        self.ymodel_p = None
        
        self.x_m = None
        self.y_m = None
        self.sy_m = None
        self.ymodel_m = None
        super().__init__(parent)
        
    def get_dof(self):
        self.dof = len(self.x_xrr) + len(self.x_p) + len(self.x_m) 
        for param in self.p:
            if self.p[param].vary:
                self.dof -= 1
                
    def define_plot_canvas(self):
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
        
        if self.x_xrr is not None and self.y_xrr is not None and self.sy_xrr is not None:
            self.ax1.errorbar(\
                self.x_xrr, self.y_xrr, self.sy_xrr, marker='.',\
                markersize=2, linestyle='None', color='#a6dba0',\
                label=self.data_path_xrr, zorder=0)
        if self.x_p is not None and self.y_p is not None and self.sy_p is not None:
            self.ax1.errorbar(\
                self.x_p, self.y_p, self.sy_p, marker='.',\
                markersize=2, linestyle='None', color='#92c5de',\
                label=self.data_path_pnr_p, zorder=0)
        if self.x_m is not None and self.y_m is not None and self.sy_m is not None:
            self.ax1.errorbar(\
                self.x_m, self.y_m, self.sy_m, marker='.',\
                markersize=2, linestyle='None', color='#f4a582',\
                label=self.data_path_pnr_m, zorder=0)

        self.model_plot_xrr, = self.ax1.plot(self.x_xrr, self.ymodel_xrr,\
             marker='None', linestyle='-', color='#008837', lw=1,\
            label="XRR", zorder=1, alpha=0.4)
        self.model_plot_plus, = self.ax1.plot(self.x_p, self.ymodel_p,\
            marker='None', linestyle='-', color='#0571b0', lw=1,\
            label="PNR I+", zorder=1, alpha=0.4)
        self.model_plot_minus, = self.ax1.plot(self.x_m, self.ymodel_m,\
            marker='None', linestyle='-', color='#ca0020', lw=1,\
            label="PNR I-", zorder=1, alpha=0.4)
        
        self.ax1.set_yscale('log')
        self.ax1.set_xlim([min(min(self.x_xrr), min(self.x_p), min(self.x_m)),\
                           max(max(self.x_xrr), max(self.x_p), max(self.x_m))])

        if self.y is not None:
            self.ax1.set_ylim(\
                [min(min(self.y_xrr), min(self.y_p), min(self.y_m))*0.8,\
                max(max(self.y_xrr), max(self.y_p), max(self.y_m))*1.2])
        else:
            self.ax1.set_ylim(\
                [min(min(self.ymodel_xrr),\
                    min(self.ymodel_p), min(self.ymodel_m),\
                    min(self.y_p), min(self.y_m))*0.8,\
                max(max(self.ymodel_xrr),\
                    max(self.ymodel_p), max(self.ymodel_m),\
                    max(self.y_p), max(self.y_m))*1.2])
        self.ax1.set_xlabel("$\mathit{q_z} \, / \, \AA^{-1}$")
        self.ax1.set_ylabel("$\mathit{I} \, / \, a.u.$")
        
        
        self.sld_xrr_plot, =\
            self.ax2.plot(self.xsld/10., np.real(self.ysld_xrr)*1e6,\
                marker='None', ls='-', color='#008837', label="$SLD_{xrr}$")
        self.sld_nuc_plot, =\
            self.ax2.plot(self.xsld/10., np.real(self.ysld_nuc)*1e6,\
                marker='None', ls='-', color='#0571b0', label="$SLD_{nuc}$")
        self.sld_mag_plot, =\
            self.ax2.plot(self.xsld/10., np.real(self.ysld_mag)*1e6,\
                marker='None', ls='-', color='#ca0020', label="$SLD_{mag}$")
        self.ax2.legend(fontsize=8).draw_frame(True)
        self.ax2.set_xlim(min(self.xsld/10.), max(self.xsld/10.))
        self.ax2.set_xlabel("$\mathit{z} \, / \, nm$")
        self.ax2.set_ylabel("$\mathit{SLD} \, / \, 10^{-6} \AA^{-2}$")
        self.fig.subplots_adjust(wspace=0.5)
        
    def update_plot(self):
        self.get_model(self.p)
        self.model_plot_xrr.set_ydata(self.ymodel_xrr)
        self.model_plot_plus.set_ydata(self.ymodel_p)
        self.model_plot_minus.set_ydata(self.ymodel_m)

        self.sld_xrr_plot.set_ydata(np.real(self.ysld_xrr)*1e6)
        self.sld_nuc_plot.set_ydata(np.real(self.ysld_nuc)*1e6)
        self.sld_mag_plot.set_ydata(np.real(self.ysld_mag)*1e6)

        fom = self.figure_of_merit(self.p)
        self.chi2 = sum(fom**2)/self.dof
        self.update_chi2()

        self.draw()
        
    def figure_of_merit(self, p):
        self.get_model(p)
        resi_xrr = (np.log(self.ymodel_xrr)-np.log(self.y_xrr))/self.sy_xrr*self.y_xrr
        resi_p = (np.log(self.ymodel_p)-np.log(self.y_p))/self.sy_p*self.y_p
        resi_m = (np.log(self.ymodel_m)-np.log(self.y_m))/self.sy_m*self.y_m
        if "weightXRR" in p:
            resi_xrr = resi_xrr * p["weightXRR"].value
        return np.concatenate([resi_xrr, resi_p, resi_m])
        
    def export_model(self):
        savefile = open(self.reflsavefile, "w")
        if self.fit_result is not None:
            savefile.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        if self.y_p is not None:
            savefile.write("#qz / A-1 \t I / a.u. \t sI / a.u."+\
                " \t Imodel / a.u.\n")
            savefile.write("#XRR\n")
            for iq, qval in enumerate(self.x_xrr):
                savefile.write(str(qval) +"\t"+ str(self.y_xrr[iq])+"\t"+\
                        str(self.sy_xrr[iq]) + "\t" + str(self.ymodel_xrr[iq])+"\n")
            savefile.write("#PNR I+\n")
            for iq, qval in enumerate(self.x_p):
                savefile.write(str(qval) +"\t"+ str(self.y_p[iq])+"\t"+\
                        str(self.sy_p[iq]) + "\t" + str(self.ymodel_p[iq])+"\n")
            savefile.write("#PNR I-\n")
            for iq, qval in enumerate(self.x_m):
                savefile.write(str(qval) +"\t"+ str(self.y_m[iq])+"\t"+\
                        str(self.sy_m[iq]) + "\t" + str(self.ymodel_m[iq])+"\n")
        else:
            savefile.write("#q / A-1 \t Imodel / a.u.\n")
            savefile.write("#XRR\n")
            for iq, qval in enumerate(self.x_xrr):
                savefile.write(str(qval) + "\t" + str(self.ymodel_xrr[iq])+"\n")
            savefile.write("#PNR I+\n")
            for iq, qval in enumerate(self.x_p):
                savefile.write(str(qval) + "\t" + str(self.ymodel_p[iq])+"\n")
            savefile.write("#PNR I-\n")
            for iq, qval in enumerate(self.x_m):
                savefile.write(str(qval) + "\t" + str(self.ymodel_m[iq])+"\n")
        print("Wrote results to " + self.reflsavefile)
        savefile.close()
        
        savefile = open(self.sldsavefile, "w")
        if self.fit_result is not None:
            savefile.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        savefile.write("\n#z / nm \t SLD_XRR / 1e-6 A-2 \t SLD_n / 1e-6 A-2 \t SLD_m / 1e-6 A-2\n")
        for iq, xval in enumerate(self.xsld):
            savefile.write(str(xval/10.) +\
                           "\t"+ str(self.ysld_xrr[iq]*1e6)+\
                           "\t"+ str(self.ysld_nuc[iq]*1e6)+\
                           "\t"+ str(self.ysld_mag[iq]*1e6)+"\n")
        print("Wrote results to " + self.sldsavefile)
        savefile.close()
        
    def save_plot(self):
        self.fig.savefig(self.reflsavefile.rsplit(".",1)[0]+"_plot.png")
