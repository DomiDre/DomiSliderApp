from PyQt5.QtWidgets import QApplication
from slider_fit_app import SliderFitApp
from slider_saxs import cPlotAndFitSAXS

import sas_models
import sas_methods
import numpy as np
import sys, lmfit

class SAXSGuiApp(cPlotAndFitSAXS):
    def init_data(self):
        self.qmin = 1e-2
        self.qmax = np.inf

        self.data_path = "finalDD67.xy"

        self.x, self.y, self.sy, q_exp_excl, I_exp_excl, sI_exp_excl =\
                sas_methods.load_saxs_data(self.data_path, self.qmin, self.qmax)

        self.p = lmfit.Parameters()
        self.p.add("R", 63.3522145282, min = 0, max = 633.0, vary = True)
        self.p.add("SLDsphere", 4.2146e-05, min = 3e-05, max = 6e-05, vary = False) # CoFe2O4
        self.p.add("SLDmatrix", 7.694e-06, min = 1e-06, max = 2e-05, vary = False) # c-Hexan
        self.p.add("sigR", 0.0883982402697, min = 0, max = 0.2, vary = True)
        self.p.add("i0", 0.0260111896747, min = 0, max = 0.03, vary = True)
        self.p.add("bg", 0.00147369, min = 0, max = 0.1, vary = False)

        # GALAXI resolution
        self.dth = 0.3e-3 # incident angle divergence from paper
        self.dlam_ov_lam = 0.005 # Ga Ka1 and Ka2 not resolved, estimate energy dispersion 
        self.wavelength = 1.3414

        self.ymodel = self.get_model(self.p, self.x)

        self.xsld = np.linspace(0, 150, 150)
        self.ysld = self.get_sld(self.p, self.xsld)
        
        
    def get_model(self, p, x):
        R = p["R"]
        SLDsphere = p["SLDsphere"]
        SLDmatrix = p["SLDmatrix"]
        sigR = p["sigR"]
        i0 = p["i0"]
        bg = p["bg"]
        Imodel = sas_models.sphere.formfactor(x, R, SLDsphere, SLDmatrix, sigR) 
        sigQ = np.sqrt((self.dlam_ov_lam * x)**2 + (4.*np.pi/self.wavelength * self.dth)**2)
        Imodel = sas_models.math.resolution_smear(x, Imodel, sigQ)
        Imodel = i0*Imodel + bg
        return Imodel
        
    def get_sld(self, p, x):
        R = p["R"]
        SLDsphere = p["SLDsphere"]
        SLDmatrix = p["SLDmatrix"]
        sigR = p["sigR"]
        sld = sas_models.sphere.sld(x, R, SLDsphere, SLDmatrix, sigR)
        return sld

if __name__ == "__main__":
    app = QApplication(sys.argv)
    aw = SliderFitApp(SAXSGuiApp)
    aw.show()
    app.exec_()
