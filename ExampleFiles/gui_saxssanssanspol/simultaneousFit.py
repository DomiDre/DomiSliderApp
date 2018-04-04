from PyQt5.QtWidgets import QApplication
from SliderApp.slider_fit_app import SliderFitApp
from SliderApp.slider_saxssanssanspol import cPlotAndFitSAXSSANSSANSPOL

from SASModels import sas_models
from SASModels import sas_methods
import numpy as np
import sys, lmfit

sld = {
    'CoFe2O4_n': 6.069e-6,
    'OA_n': 0.078e-6,
    'd8Toluene_n': 5.664e-06,
    'CoFe2O4_xr': 40.825e-6,
    'OA_xr': 8.472e-6,
    'cHexane_xr': 6.461e-6
}

class SANSPOLGuiApp(cPlotAndFitSAXSSANSSANSPOL):
    def init_data(self):
        self.qmin = 1e-2
        self.qmax = np.inf

        self.data_path_saxs='./1SAXS/0rawdata/finalDD67.xy'
        self.data_path_sans_sa='./2SANS/SANS_nuclear/0rawdata/DD67_nuclear20_sa.dat'
        self.data_path_sans_la='./2SANS/SANS_nuclear/0rawdata/DD67_nuclear20_la_scaled.dat'
        self.data_pathp_sanspol_sa='./2SANS/SANSPOL_magnetic/0rawdata/dd67_rfm_sa_scaled.dat'
        self.data_pathp_sanspol_la='./2SANS/SANSPOL_magnetic/0rawdata/dd67_rfm_la.dat'
        self.data_pathm_sanspol_sa='./2SANS/SANSPOL_magnetic/0rawdata/dd67_rfp_sa_scaled.dat'
        self.data_pathm_sanspol_la='./2SANS/SANSPOL_magnetic/0rawdata/dd67_rfp_la.dat'

        self.modelfile = "SAXS_SANSPOL_DD67_fit.dat"
        self.sldmodelfile = "SAXS_SANSPOL_DD67_sld.dat"

        self.x_saxs, self.y_saxs, self.sy_saxs, q_exp_excl, I_exp_excl, sI_exp_excl =\
                sas_methods.load_sas_data(self.data_path_saxs, self.qmin, self.qmax)

        self.x_sans_sa, self.y_sans_sa, self.sy_sans_sa, q_exp_excl, I_exp_excl, sI_exp_excl =\
                sas_methods.load_sas_data(self.data_path_sans_sa, self.qmin, self.qmax)

        self.x_sans_la, self.y_sans_la, self.sy_sans_la, q_exp_excl, I_exp_excl, sI_exp_excl =\
                sas_methods.load_sas_data(self.data_path_sans_la, self.qmin, self.qmax)

        self.xp_sanspol_sa, self.yp_sanspol_sa, self.syp_sanspol_sa, q_exp_excl, I_exp_excl, sI_exp_excl =\
                sas_methods.load_sas_data(self.data_pathp_sanspol_sa, self.qmin, self.qmax)

        self.xp_sanspol_la, self.yp_sanspol_la, self.syp_sanspol_la, q_exp_excl, I_exp_excl, sI_exp_excl =\
                sas_methods.load_sas_data(self.data_pathp_sanspol_la, self.qmin, self.qmax)

        self.xm_sanspol_sa, self.ym_sanspol_sa, self.sym_sanspol_sa, q_exp_excl, I_exp_excl, sI_exp_excl =\
                sas_methods.load_sas_data(self.data_pathm_sanspol_sa, self.qmin, self.qmax)

        self.xm_sanspol_la, self.ym_sanspol_la, self.sym_sanspol_la, q_exp_excl, I_exp_excl, sI_exp_excl =\
                sas_methods.load_sas_data(self.data_pathm_sanspol_la, self.qmin, self.qmax)

        self.p = lmfit.Parameters()
        self.p.add("weightSAXS", 0.1, min = 0, max = 1, vary = False)
        self.p.add("i0_xr", 0.025480498736948187, min = 0, max = 0.05, vary = True)
        self.p.add("bg_xr", 0.00132, min = 0, max = 0.03, vary = False)
        self.p.add("i0_n", 0.214, min = 0, max = 1, vary = True)
        self.p.add("bg_n", 0.016139999999999998, min = 0, max = 0.03, vary = False)
        self.p.add("dth_sa", 0.001297, min = 0.001, max = 0.004, vary = False)
        self.p.add("dth_la", 0.001372, min = 0.001, max = 0.004, vary = False)
        self.p.add("R", 63.09230380086919, min = 0, max = 70, vary = True)
        self.p.add("dshell", 15.511037377491995, min = 0, max = 30, vary = True)
        self.p.add("sigR", 0.0906845441066217, min = 0, max = 0.2, vary = True)
        self.p.add("sigdshell", 0.0, min = 0, max = 0.2, vary = False)
        self.p.add("magSLDcore", 2.7399999999999994e-07, min = 0, max = 2e-06, vary = True)

        # D33 options
        self.xi = 1
        self.sin2alpha = 0.98990776678064341
        self.wavelength = 6.0 #A
        self.dlam_ov_lam = 4.247e-3

        self.get_model(self.p)

        self.xsld = np.linspace(0, 150, 150)
        self.get_sld(self.p)
        
        
        
    def get_model(self, p):
        R = p["R"]
        dshell = p["dshell"]
        sigR = p["sigR"]
        sigdshell = p["sigdshell"]
        magSLDcore = p["magSLDcore"]
        
        i0_xr = p["i0_xr"]
        bg_xr = p["bg_xr"]
        i0_n = p["i0_n"]
        bg_n = p["bg_n"]
        dth_sa = p["dth_sa"]
        dth_la = p["dth_la"]

        SLDcore_n = sld['CoFe2O4_n']
        SLDshell_n = sld['OA_n']
        SLDmatrix_n = sld['d8Toluene_n']
        SLDcore_xr = sld['CoFe2O4_xr']
        SLDshell_xr = sld['OA_xr']
        SLDmatrix_xr = sld['cHexane_xr']
        
        Imodel_saxs = sas_models.coreshell.formfactor(\
                        self.x_saxs, R, dshell,\
                        SLDcore_xr, SLDshell_xr, SLDmatrix_xr,\
                        sigR, sigdshell)
        Imodel_sans_sa = sas_models.coreshell.formfactor(\
                        self.x_sans_sa, R, dshell,\
                        SLDcore_n, SLDshell_n, SLDmatrix_n,\
                        sigR, sigdshell)
        Imodel_sans_la = sas_models.coreshell.formfactor(\
                        self.x_sans_la, R, dshell,\
                        SLDcore_n, SLDshell_n, SLDmatrix_n,\
                        sigR, sigdshell)
        Imodelp_sanspol_sa = sas_models.coreshell.magnetic_formfactor(\
                        self.xp_sanspol_sa, R, dshell,\
                        SLDcore_n, SLDshell_n, SLDmatrix_n,\
                        sigR, sigdshell, 0, magSLDcore, 0, 0, self.xi,\
                        self.sin2alpha, 1)
        Imodelp_sanspol_la = sas_models.coreshell.magnetic_formfactor(\
                        self.xp_sanspol_la, R, dshell,\
                        SLDcore_n, SLDshell_n, SLDmatrix_n,\
                        sigR, sigdshell, 0, magSLDcore, 0, 0, self.xi,\
                        self.sin2alpha, 1)
        Imodelm_sanspol_sa = sas_models.coreshell.magnetic_formfactor(\
                        self.xm_sanspol_sa, R, dshell,\
                        SLDcore_n, SLDshell_n, SLDmatrix_n,\
                        sigR, sigdshell, 0, magSLDcore, 0, 0, self.xi,\
                        self.sin2alpha, -1) 
        Imodelm_sanspol_la = sas_models.coreshell.magnetic_formfactor(\
                        self.xm_sanspol_la, R, dshell,\
                        SLDcore_n, SLDshell_n, SLDmatrix_n,\
                        sigR, sigdshell, 0, magSLDcore, 0, 0, self.xi,\
                        self.sin2alpha, -1) 
                        
        self.sigQ_sans_sa = np.sqrt((self.dlam_ov_lam * self.x_sans_sa)**2 +\
                                         (4.*np.pi/self.wavelength * dth_sa)**2)
        self.sigQ_sans_la = np.sqrt((self.dlam_ov_lam * self.x_sans_la)**2 +\
                                         (4.*np.pi/self.wavelength * dth_la)**2)
        self.sigQp_sanspol_sa = np.sqrt((self.dlam_ov_lam * self.xp_sanspol_sa)**2 +\
                                         (4.*np.pi/self.wavelength * dth_sa)**2)
        self.sigQp_sanspol_la = np.sqrt((self.dlam_ov_lam * self.xp_sanspol_la)**2 +\
                                         (4.*np.pi/self.wavelength * dth_la)**2)
        self.sigQm_sanspol_sa = np.sqrt((self.dlam_ov_lam * self.xm_sanspol_sa)**2 +\
                                         (4.*np.pi/self.wavelength * dth_sa)**2)
        self.sigQm_sanspol_la = np.sqrt((self.dlam_ov_lam * self.xm_sanspol_la)**2 +\
                                         (4.*np.pi/self.wavelength * dth_la)**2)

        Imodel_sans_sa = sas_models.math.resolution_smear_interpolating(self.x_sans_sa,\
                                        Imodel_sans_sa, self.sigQ_sans_sa)
        Imodel_sans_la = sas_models.math.resolution_smear_interpolating(self.x_sans_la,\
                                        Imodel_sans_la, self.sigQ_sans_la)
        Imodelp_sanspol_sa = sas_models.math.resolution_smear_interpolating(self.xp_sanspol_sa,\
                                        Imodelp_sanspol_sa, self.sigQp_sanspol_sa)
        Imodelp_sanspol_la = sas_models.math.resolution_smear_interpolating(self.xp_sanspol_la,\
                                        Imodelp_sanspol_la, self.sigQp_sanspol_la)
        Imodelm_sanspol_sa = sas_models.math.resolution_smear_interpolating(self.xm_sanspol_sa,\
                                        Imodelm_sanspol_sa, self.sigQm_sanspol_sa)
        Imodelm_sanspol_la = sas_models.math.resolution_smear_interpolating(self.xm_sanspol_la,\
                                        Imodelm_sanspol_la, self.sigQm_sanspol_la)

        self.ymodel_saxs = i0_xr*Imodel_saxs + bg_xr
        self.ymodel_sans_sa = i0_n*Imodel_sans_sa + bg_n
        self.ymodel_sans_la = i0_n*Imodel_sans_la + bg_n
        self.ymodelp_sanspol_sa = i0_n*Imodelp_sanspol_sa + bg_n
        self.ymodelp_sanspol_la = i0_n*Imodelp_sanspol_la + bg_n
        self.ymodelm_sanspol_sa = i0_n*Imodelm_sanspol_sa + bg_n
        self.ymodelm_sanspol_la = i0_n*Imodelm_sanspol_la + bg_n
        
        
    def get_sld(self, p):
        R = p["R"]
        dshell = p["dshell"]
        SLDcore_n = sld['CoFe2O4_n']
        SLDshell_n = sld['OA_n']
        SLDmatrix_n = sld['d8Toluene_n']
        magSLDcore = p["magSLDcore"]
        SLDcore_xr = sld['CoFe2O4_xr']
        SLDshell_xr = sld['OA_xr']
        SLDmatrix_xr = sld['cHexane_xr']

        self.ysaxssld = sas_models.coreshell.sld(self.xsld, R, dshell,\
                                          SLDcore_xr, SLDshell_xr, SLDmatrix_xr) 
        self.ynucsld = sas_models.coreshell.sld(self.xsld, R, dshell,\
                                          SLDcore_n, SLDshell_n, SLDmatrix_n) 
        self.ymagsld = sas_models.coreshell.sld(self.xsld, R, dshell,\
                                          magSLDcore, 0, 0) 

if __name__ == "__main__":
    app = QApplication(sys.argv)
    aw = SliderFitApp(SANSPOLGuiApp)
    aw.show()
    app.exec_()
