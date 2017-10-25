from PyQt5.QtWidgets import QApplication
from DomiSliderApp.slider_fit_app import SliderFitApp
from DomiSliderApp.slider_refl import cPlotAndFitRefl

import ReflNanoparticles.reflectivity as reflectivity
import numpy as np
import sys, lmfit

class ReflGuiApp(cPlotAndFitRefl):
    def init_data(self):
        self.qmin = 1e-2
        self.qmax = np.inf

        self.data_path = "DD101_7_bruker.xy"

        self.x = np.linspace(0.001, 0.1, 100)
        self.y = np.linspace(0.001, 0.1, 100)*1
        self.sy = np.linspace(0.001, 0.1, 100)*1
        
        self.p = lmfit.Parameters()
        self.p.add('dai',              1e-4,          min=0.,   max=1e-3, vary=0)
        self.p.add('I0',               4.585e8,       min=1e7,  max=1e9,  vary=1)
        self.p.add('Ibg',              70,            min=0.,   max=100,  vary=0)
        self.p.add('th_offset',        0.0,           min=-0.1,  max=0.1, vary=0)
        self.p.add('roughness',        5.82017266161, min=0.0,  max=20,   vary=0)
        self.p.add('oleic_sld',        8.467*1e-6,    min=0.,   max=2e-5, vary=0)
        self.p.add('particle_size_1',  108.676202248, min=90.0, max=200,  vary=0)
        self.p.add('particle_size_2',  30.0,          min=0.0,  max=200,  vary=0)
        self.p.add('packing_density_1',0.584971550956,min=0.0,  max=1.0,  vary=0)
        self.p.add('packing_density_2',0.2,           min=0.0,  max=1.0,  vary=0)
        self.p.add('oleic_spacer1',    15.0834474946, min=0.,   max=50,   vary=0)
        self.p.add('oleic_spacer2',    15.0834474946, min=0.,   max=50,   vary=0)
#        self.p.add('particle_sld_real',42.146*1e-6,   min=3e-5, max=5e-5, vary=0)
#        self.p.add('particle_sld_imag',3.121e-6,      min=-1e-5,max=1e-5, vary=0)
        # Bruker properties
        self.dlam_ov_lam = 0.01 #? 
        self.wavelength = 1.3414 # Cu-K-alpha 
        self.beamwidth = 0.2 # 200 microns
        self.samplelen = 10 # 10mm
        
        self.xsld = np.linspace(-50, 400, 450)
        self.particle_sld = 42.146*1e-6 - 1j*3.121e-6
        self.si_sld = 20.065*1E-6 - 1j*0.351e-6 # A-2
        
        self.get_model(self.p)

        
    def get_model(self, p):
        packing_density_1 = p["packing_density_1"]
        packing_density_2 = p["packing_density_2"]
        particle_size_1 = p["particle_size_1"]
        particle_size_2 = p["particle_size_2"]
        oleic_spacer1 = p["oleic_spacer1"]
        oleic_spacer2 = p["oleic_spacer2"]
        oleic_sld = p["oleic_sld"]
#        self.particle_sld = p['particle_sld_real'] -1j*p['particle_sld_imag']
        
        sld = np.array([self.si_sld,\
                        oleic_sld,
                        self.particle_sld*packing_density_1,\
                        oleic_sld,
                        self.particle_sld*packing_density_2,
                        0.])
        
        roughness = p["roughness"]*np.ones(len(sld))
        thickness = np.array([0,\
                        oleic_spacer1,
                        particle_size_1,\
                        oleic_spacer2,\
                        particle_size_2,\
                        10])
        
        q = reflectivity.math.angle_correct(self.x, p["th_offset"],\
                            self.wavelength)
                            
        self.ymodel = reflectivity.models.parrat(q, sld, roughness,\
                            thickness)
                            
        self.ymodel = reflectivity.math.footprint_correct(q, self.ymodel,\
                            self.samplelen, self.beamwidth, \
                            self.wavelength)
        
        
        sigQ = np.sqrt((self.dlam_ov_lam * q)**2 +\
                    (4.*np.pi/self.wavelength * p["dai"])**2)
        self.ymodel = reflectivity.math.resolution_smear(q, self.ymodel, sigQ)
        self.ymodel = p["I0"]*self.ymodel + p["Ibg"]
        
        self.ysld = reflectivity.models.roughsld_thick_layers(self.xsld,\
                                sld, roughness, thickness)
                                
        
if __name__ == "__main__":
    app = QApplication(sys.argv)
    aw = SliderFitApp(ReflGuiApp)
    aw.show()
    app.exec_()
