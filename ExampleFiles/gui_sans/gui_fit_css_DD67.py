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
        self.p.add("R", 62.5306, min=0, max=100, vary=0)
        self.p.add("dshell", 0., min=0, max=100, vary=1)
        self.p.add("dsurfactant", 8.36927, min=0, max = 30, vary=1)
        self.p.add("SLDcore", 7.519e-6,min=0, max=9e-6,vary=0) # CoFe2O4
        self.p.add("SLDshell", 6.069e-6, min=0,  max=9e-6, vary=0) # CoFe2O4, dens = 5.3 g/mL
        self.p.add("SLDsurfactant", 0.078e-6, min=0, max=5e-6, vary=0) # Oleic Acid
        self.p.add("SLDmatrix", 5.664e-6, min=1e-6, vary=False) # c-Hexan
        self.p.add("sigR", 0.0959, min=0, max=0.2)
        self.p.add("sigdshell", 0., min=0, max=0.2, vary=0)
        self.p.add("sigdsurfactant", 0., min=0, max=0.2, vary=0)
        self.p.add("dth_sa", 2e-3, min=0, max=5e-3, vary=1)
        self.p.add("dth_la", 2e-3, min=0, max=5e-3, vary=1)
        self.p.add("i0", 0.9423, min=0, max=2, vary=1)
        self.p.add("bg", 0.0176, min=0, max=0.03, vary=False)

        # D33 resolution
        self.wavelength = 6.0 #A
        self.dlam_ov_lam = 4.247e-3


        self.ymodel_sa, self.ymodel_la = self.get_model(self.p,\
                                                        self.x_sa, self.x_la)

        self.xsld = np.linspace(0, 150, 150)
        self.ysld = self.get_sld(self.p, self.xsld)
        
        
        
    def get_model(self, p, x_sa, x_la):
        R = p["R"]
        dshell = p["dshell"]
        dsurfactant = p["dsurfactant"]
        SLDcore = p["SLDcore"]
        SLDshell = p["SLDshell"]
        SLDsurfactant = p["SLDsurfactant"]
        SLDmatrix = p["SLDmatrix"]
        sigR = p["sigR"]
        sigdshell = p["sigdshell"]
        sigdsurfactant = p["sigdsurfactant"]
        dth_sa = p["dth_sa"]
        dth_la = p["dth_la"]
        i0 = p["i0"]
        bg = p["bg"]
        Imodel_sa = sas_models.coreshellsurfactant.formfactor(x_sa, R, dshell,\
                    dsurfactant, SLDcore, SLDshell, SLDsurfactant, SLDmatrix,\
                    sigR, sigdshell, sigdsurfactant) 
        Imodel_la = sas_models.coreshellsurfactant.formfactor(x_la,R, dshell,\
                    dsurfactant, SLDcore, SLDshell, SLDsurfactant, SLDmatrix,\
                    sigR, sigdshell, sigdsurfactant) 
                    
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
        dshell = p["dshell"]
        dsurfactant = p["dsurfactant"]
        SLDcore = p["SLDcore"]
        SLDshell = p["SLDshell"]
        SLDsurfactant = p["SLDsurfactant"]
        SLDmatrix = p["SLDmatrix"]
        sigR = p["sigR"]
        sigdshell = p["sigdshell"]
        sigdsurfactant = p["sigdsurfactant"]
        sld = sas_models.coreshellsurfactant.sld(x, R, dshell,\
            dsurfactant, SLDcore, SLDshell, SLDsurfactant, SLDmatrix,\
            sigR, sigdshell, sigdsurfactant) 
        return sld

if __name__ == "__main__":
    app = QApplication(sys.argv)
    aw = SliderFitApp(SANSGuiApp)
    aw.show()
    app.exec_()
