from PyQt5.QtWidgets import QApplication
from slider_fit_app import SliderFitApp
from slider_sanspol import cPlotAndFitSANSPOL

import sas_models
import sas_methods
import numpy as np
import sys, lmfit

class SANSPOLGuiApp(cPlotAndFitSANSPOL):
    def init_data(self):
        self.qmin = 1e-3
        self.qmax = 1

        self.xp_sa = np.linspace(self.qmin, self.qmax, 300)
        self.xp_la = np.linspace(self.qmin, self.qmax, 300)
        self.xm_sa = np.linspace(self.qmin, self.qmax, 300)
        self.xm_la = np.linspace(self.qmin, self.qmax, 300)

        self.p = lmfit.Parameters()
        self.p.add("R1", 62.5306, min=0, max=70, vary=0)
        self.p.add("dR", 15.32, min=0, max=30, vary=1)
        self.p.add("SLDcore", 6.069e-6,min=5e-6, max=8e-6,vary=0) # CoFe2O4
        self.p.add("SLDshell", 2.2259e-6, min=0,  max=8e-6, vary=1) # CoFe2O4, dens = 5.3 g/mL
        self.p.add("SLDmatrix", 5.664e-6, min=5e-6, max=7e-6,vary=0) # c-Hexan
        self.p.add("sigR", 0.0942, min=0, max=0.2, vary=0)
        self.p.add("sigdR", 0., min=0, max=0.2, vary=0)
        self.p.add("magSLDcore", 1e-6, min=0, max=2e-6, vary=1)
        self.p.add("magSLDshell", 1e-6, min=0, max=2e-6, vary=1)
        self.p.add("i0", 0.9423, min=0, max=2, vary=1)
        self.p.add("bg", 0.0176, min=0, max=0.03, vary=False)
        self.p.add("dth_sa", 1e-3, min=0, max=4e-3, vary=False)
        self.p.add("dth_la", 1e-3, min=0, max=4e-3, vary=False)

        # D33 options
        self.xi = 1
        self.sin2alpha = 0.98990776678064341
        self.wavelength = 6.0 #A
        self.dlam_ov_lam = 4.247e-3

        self.ymodelp_sa, self.ymodelp_la, self.ymodelm_sa, self.ymodelm_la =\
                self.get_model(self.p, self.xp_sa, self.xp_la, self.xm_sa,\
                self.xm_la)

        self.xsld = np.linspace(0, 150, 150)
        self.ynucsld, self.ymagsld = self.get_sld(self.p, self.xsld)
        
        
        
    def get_model(self, p, xp_sa, xp_la, xm_sa, xm_la):
        R1 = p["R1"]
        dR = p["dR"]
        SLDcore = p["SLDcore"]
        SLDshell = p["SLDshell"]
        SLDmatrix = p["SLDmatrix"]
        sigR = p["sigR"]
        sigdR = p["sigdR"]
        magSLDcore = p["magSLDcore"]
        magSLDshell = p["magSLDshell"]
        dth_sa = p["dth_sa"]
        dth_la = p["dth_la"]
        i0 = p["i0"]
        bg = p["bg"]
        
        Imodelp_sa = sas_models.sphere_linear_hull.magnetic_formfactor(\
                        xp_sa, R1, dR,\
                        SLDcore, SLDshell, SLDmatrix,\
                        sigR, sigdR, magSLDcore, magSLDshell, 0, self.xi,\
                        self.sin2alpha, 1)
        Imodelp_la = sas_models.sphere_linear_hull.magnetic_formfactor(\
                        xp_la, R1, dR,\
                        SLDcore, SLDshell, SLDmatrix,\
                        sigR, sigdR, magSLDcore, magSLDshell, 0, self.xi,\
                        self.sin2alpha, 1)
        Imodelm_sa = sas_models.sphere_linear_hull.magnetic_formfactor(\
                        xm_sa, R1, dR,\
                        SLDcore, SLDshell, SLDmatrix,\
                        sigR, sigdR, magSLDcore, magSLDshell, 0, self.xi,\
                        self.sin2alpha, -1) 
        Imodelm_la = sas_models.sphere_linear_hull.magnetic_formfactor(\
                        xm_la, R1, dR,\
                        SLDcore, SLDshell, SLDmatrix,\
                        sigR, sigdR, magSLDcore, magSLDshell, 0, self.xi,\
                        self.sin2alpha, -1) 
                        
        self.sigQp_sa = np.sqrt((self.dlam_ov_lam * self.xp_sa)**2 +\
                                         (4.*np.pi/self.wavelength * dth_sa)**2)
        self.sigQp_la = np.sqrt((self.dlam_ov_lam * self.xp_la)**2 +\
                                         (4.*np.pi/self.wavelength * dth_la)**2)
        self.sigQm_sa = np.sqrt((self.dlam_ov_lam * self.xm_sa)**2 +\
                                         (4.*np.pi/self.wavelength * dth_sa)**2)
        self.sigQm_la = np.sqrt((self.dlam_ov_lam * self.xm_la)**2 +\
                                         (4.*np.pi/self.wavelength * dth_la)**2)

        Imodelp_sa = sas_models.math.resolution_smear_interpolating(xp_sa,\
                                        Imodelp_sa, self.sigQp_sa)
        Imodelp_la = sas_models.math.resolution_smear_interpolating(xp_la,\
                                        Imodelp_la, self.sigQp_la)
        Imodelm_sa = sas_models.math.resolution_smear_interpolating(xm_sa,\
                                        Imodelm_sa, self.sigQm_sa)
        Imodelm_la = sas_models.math.resolution_smear_interpolating(xm_la,\
                                        Imodelm_la, self.sigQm_la)

        Imodelp_sa = i0*Imodelp_sa + bg
        Imodelp_la = i0*Imodelp_la + bg
        Imodelm_sa = i0*Imodelm_sa + bg
        Imodelm_la = i0*Imodelm_la + bg
        
        return Imodelp_sa, Imodelp_la, Imodelm_sa, Imodelm_la
        
    def get_sld(self, p, x):
        R1 = p["R1"]
        dR = p["dR"]
        SLDcore = p["SLDcore"]
        SLDshell = p["SLDshell"]
        SLDmatrix = p["SLDmatrix"]
        magSLDcore = p["magSLDcore"]
        magSLDshell = p["magSLDshell"]
        sigR = p["sigR"]
        sigdR = p["sigdR"]
        sldnuc = sas_models.sphere_linear_hull.sld(x, R1, dR,\
                                                SLDcore, SLDshell, SLDmatrix,\
                                                sigR, sigdR) 
        sldmag = sas_models.sphere_linear_hull.sld(x, R1, dR,\
                                                magSLDcore, magSLDshell, 0,\
                                                sigR, sigdR) 
        return sldnuc, sldmag

if __name__ == "__main__":
    app = QApplication(sys.argv)
    aw = SliderFitApp(SANSPOLGuiApp)
    aw.show()
    app.exec_()
