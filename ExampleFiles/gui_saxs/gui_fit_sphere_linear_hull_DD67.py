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
        self.p.add("R2", 63.0131115636, min = 0, max = 90.0, vary = True)
        self.p.add("dR", 72.17999999999999, min = 0, max = 90.0, vary = True)
        self.p.add("SLDcore", 4.1490000000000004e-05, min = 3e-05, max = 6e-05, vary = True) # CoFe2O4
        self.p.add("SLDhull", 4.59871113608e-05, min = 3e-05, max = 6e-05, vary = True) # CoFe2O4
        self.p.add("SLDmatrix", 7.694e-06, min = 1e-06, max = 2e-05, vary = False) # c-Hexan
        self.p.add("sigR2", 0.08900000000000001, min = 0, max = 0.2, vary = True)
        self.p.add("sigdR", 0.0, min = 0, max = 0.2, vary = False)
        self.p.add("i0", 0.018209999999999997, min = 0, max = 0.03, vary = True)
        self.p.add("bg", 0.0011, min = 0, max = 0.1, vary = True)

        # GALAXI resolution
        self.dth = 0.3e-3 # incident angle divergence from paper
        self.dlam_ov_lam = 0.005 # Ga Ka1 and Ka2 not resolved, estimate energy dispersion 
        self.wavelength = 1.3414

        self.ymodel = self.get_model(self.p, self.x)

        self.xsld = np.linspace(0, 150, 150)
        self.ysld = self.get_sld(self.p, self.xsld)
        
        
    def get_model(self, p, x):
        R2 = p["R2"]
        dR = p["dR"]
        SLDcore = p["SLDcore"]
        SLDhull = p["SLDhull"]
        SLDmatrix = p["SLDmatrix"]
        sigR2 = p["sigR2"]
        sigdR = p["sigdR"]
        i0 = p["i0"]
        bg = p["bg"]
        Imodel = sas_models.sphere_linear_hull.formfactor(x, R2, dR, SLDcore, SLDhull,\
                                                         SLDmatrix, sigR2, sigdR) 
        sigQ = np.sqrt((self.dlam_ov_lam * x)**2 + (4.*np.pi/self.wavelength * self.dth)**2)
        Imodel = sas_models.math.resolution_smear(x, Imodel, sigQ)
        Imodel = i0*Imodel + bg
        return Imodel
        
    def get_sld(self, p, x):
        R2 = p["R2"]
        dR = p["dR"]
        SLDcore = p["SLDcore"]
        SLDhull = p["SLDhull"]
        SLDmatrix = p["SLDmatrix"]
        sigR2 = p["sigR2"]
        sigdR = p["sigdR"]
        sld = sas_models.sphere_linear_hull.sld(x, R2, dR, SLDcore, SLDhull,\
                                                         SLDmatrix, sigR2, sigdR)
        return sld

if __name__ == "__main__":
    app = QApplication(sys.argv)
    aw = SliderFitApp(SAXSGuiApp)
    aw.show()
    app.exec_()
