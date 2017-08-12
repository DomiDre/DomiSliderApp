from PyQt5.QtWidgets import QApplication
from slider_fit_app import SliderFitApp
from slider_sans import cPlotAndFitSANS

import sas_models
import sas_methods
import numpy as np
import sys, lmfit

class SANSGuiApp(cPlotAndFitSANS):
    def init_data(self):
        self.qmin = 1e-2
        self.qmax = np.inf

        self.data_path_sa='DD67_nuclear20_sa.dat'
        self.data_path_la='DD67_nuclear20_la_scaled.dat'

        self.x_sa, self.y_sa, self.sy_sa, q_exp_excl, I_exp_excl, sI_exp_excl =\
                sas_methods.load_sas_data(self.data_path_sa, self.qmin, self.qmax)

        self.x_la, self.y_la, self.sy_la, q_exp_excl, I_exp_excl, sI_exp_excl =\
                sas_methods.load_sas_data(self.data_path_la, self.qmin, self.qmax)

        self.p = lmfit.Parameters()
        self.p.add("R", 62.5306, min=0, max=100)
        self.p.add("SLDsphere", 42.146e-6,min=30e-6, max=60e-6,vary=0) # CoFe2O4
        self.p.add("SLDmatrix", 7.694e-6, min=1e-6, max=20e-6,vary=False) # c-Hexan
        self.p.add("sigR", 0.0959, min=0, max=0.2)
        self.p.add("i0", 0.02733847, min=0, max= 3e-2)
        self.p.add("dth_sa", 2e-3, min=0, max=5e-3, vary=1)
        self.p.add("dth_la", 2e-3, min=0, max=5e-3, vary=1)
        self.p.add("bg", 0.00147369, min=0, max=1e-1, vary=0)
        # D33 resolution
        self.wavelength = 6.0 #A
        self.dlam_ov_lam = 4.247e-3

        self.ymodel_sa, self.ymodel_la = self.get_model(self.p,\
                                                        self.x_sa, self.x_la)

        self.xsld = np.linspace(0, 150, 150)
        self.ysld = self.get_sld(self.p, self.xsld)
        
    def get_model(self, p, x_sa, x_la):
        R = p["R"]
        SLDsphere = p["SLDsphere"]
        SLDmatrix = p["SLDmatrix"]
        sigR = p["sigR"]
        dth_sa = p["dth_sa"]
        dth_la = p["dth_la"]
        i0 = p["i0"]
        bg = p["bg"]
        Imodel_sa = sas_models.sphere.formfactor(x_sa, R,\
                                                 SLDsphere, SLDmatrix, sigR) 
        Imodel_la = sas_models.sphere.formfactor(x_la, R,\
                                                 SLDsphere, SLDmatrix, sigR) 

        sigQ_sa = np.sqrt((self.dlam_ov_lam * x_sa)**2 +\
                                         (4.*np.pi/self.wavelength * dth_sa)**2)
        sigQ_la = np.sqrt((self.dlam_ov_lam * x_la)**2 +\
                                         (4.*np.pi/self.wavelength * dth_la)**2)

        Imodel_sa = sas_models.math.resolution_smear(x_sa, Imodel_sa,\
                                                        sigQ_sa)
        Imodel_la = sas_models.math.resolution_smear(x_la, Imodel_la,\
                                                        sigQ_la)
                                                        
        Imodel_sa = i0*Imodel_sa + bg
        Imodel_la = i0*Imodel_la + bg
        
        return Imodel_sa, Imodel_la
        
    def get_sld(self, p, x):
        R = p["R"]
        SLDsphere = p["SLDsphere"]
        SLDmatrix = p["SLDmatrix"]
        sigR = p["sigR"]
        sld = sas_models.sphere.sld(x, R, SLDsphere, SLDmatrix, sigR)
        return sld

if __name__ == "__main__":
    app = QApplication(sys.argv)
    aw = SliderFitApp(SANSGuiApp)
    aw.show()
    app.exec_()
