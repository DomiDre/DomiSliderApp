import matplotlib
matplotlib.use("Qt5Agg")
from PyQt5 import QtCore
import PyQt5.QtWidgets as pyqt5widget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas,\
        NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import warnings, lmfit, sys, os, datetime
import sas_methods
import matplotlib.pyplot as plt
import numpy as np

# remove some annoying deprecation warnings
warnings.filterwarnings("ignore", category=UserWarning, module='matplotlib')

class SliderFitApp(pyqt5widget.QMainWindow):
    def __init__(self, PlotClass):
        super().__init__()
        
        self.version = 0.6

        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        
#       Set up menubar
        self.file_menu = pyqt5widget.QMenu('&File', self)
        self.file_menu.addAction('&Quit', self.fileQuit, 'Ctrl+C')
        self.menuBar().addMenu(self.file_menu)

        self.help_menu = pyqt5widget.QMenu('&Help', self)
        self.help_menu.addAction('&About', self.about)
        self.menuBar().addSeparator()
        self.menuBar().addMenu(self.help_menu)
        
#        Set up mainwindow
        self.main_widget = pyqt5widget.QWidget(self)
        
        self.plot_window = PlotClass(self)

        mpl_toolbar = NavigationToolbar(self.plot_window, self)
        b_refresh = pyqt5widget.QPushButton("Refresh")

        slider_widgets = pyqt5widget.QWidget(self) 

        slider_layout = pyqt5widget.QGridLayout(slider_widgets)
        self.sliders = {}
        self.checkboxes = {}
        for i, parameter in enumerate(self.plot_window.p):
            
            slider_label = pyqt5widget.QLabel(parameter)
            slider_bar = pyqt5widget.QSlider(QtCore.Qt.Horizontal, self)
            checkbox = pyqt5widget.QCheckBox(self)
            
            cur_param = self.plot_window.p[parameter]

            slider_bar.setRange(0, 1000)
            slider_bar.setTickInterval(5)
            slider_bar.setSingleStep(1)
            slider_bar.setPageStep(10)
            
            curval = cur_param.value
            minval = cur_param.min
            maxval = cur_param.max
            
            if minval == -np.inf:
                if curval > 0:
                    minval = 0
                elif curval < 0:
                    minval = 10*curval
                else:
                    minval = -1
                cur_param.min = minval

            if maxval == np.inf:
                if curval > 0:
                    maxval = 10*curval
                elif curval < 0:
                    minval = 0
                else:
                    minval = 1
                cur_param.max = maxval

            delta = (maxval - minval)/1000.
            checkbox.setChecked(cur_param.vary)
            slider_value = int((curval-minval)/delta)
            slider_bar.setValue(slider_value)
            new_value = minval + slider_bar.value()*delta
            if new_value > 1e3 or new_value < 1e-3:
                prec = '{:.3e}'
            else:
                prec = '{:.3f}' 
            slider_bar.label = pyqt5widget.QLabel(prec.format(new_value))
            
            slider_bar.valueChanged.connect(self.slider_value_changed)
            slider_layout.addWidget(slider_label, i, 0)
            slider_layout.addWidget(slider_bar, i, 1)
            slider_layout.addWidget(slider_bar.label, i, 2)
            slider_layout.addWidget(checkbox, i, 3)
            self.sliders[parameter] = slider_bar
            self.checkboxes[parameter] = checkbox
        self.sliders_inverse = dict(zip(self.sliders.values(),self.sliders.keys()))
        
#        slider_layout.setColumnMinimumWidth(1,300)
#        slider_layout.setColumnMinimumWidth(2,100)

        
        button_widget = pyqt5widget.QWidget(self)
        button_layout = pyqt5widget.QGridLayout(button_widget)
        
        self.label_chi2 = pyqt5widget.QLabel("")
        self.but_globalfit = pyqt5widget.QPushButton("Global Fit (Differential Evolution)", self)
        self.but_globalfit.setToolTip("Fit parameters set on vary in code.")
        self.but_globalfit.clicked.connect(self.plot_window.fit_global)
        self.but_localfit = pyqt5widget.QPushButton("Local Fit (LM)", self)
        self.but_localfit.setToolTip("Fit parameters set on vary in code.")
        self.but_localfit.clicked.connect(self.plot_window.fit_local)
        self.but_fit_without_bounds = pyqt5widget.QPushButton("Local Fit Without Bounds (LM)", self)
        self.but_fit_without_bounds.setToolTip("Run fit algorithm ignoring bounds from slider. Removes bias.")
        self.but_fit_without_bounds.clicked.connect(self.plot_window.fit_local_wo_bounds)
        but_saveparascript = pyqt5widget.QPushButton("Save Parameters to Script", self)
        but_saveparascript.setToolTip("Overwrite parameters in script with set values.")
        but_saveparascript.clicked.connect(self.plot_window.save_para_to_script)
#        but_export = pyqt5widget.QPushButton("Export", self)
#        but_export.setToolTip("Export parameter values to file.")
#        but_export.clicked.connect(self.plot_window.export_params)
        but_exportmodel = pyqt5widget.QPushButton("Export Model", self)
        but_exportmodel.setToolTip("Export Model.")
        but_exportmodel.clicked.connect(self.plot_window.export_model)
        button_layout.addWidget(self.label_chi2, 0, 0)
        button_layout.addWidget(self.but_globalfit, 1, 0)
        button_layout.addWidget(self.but_localfit, 2, 0)
        button_layout.addWidget(self.but_fit_without_bounds, 3, 0)
        button_layout.addWidget(but_saveparascript, 1, 1)
        button_layout.addWidget(but_exportmodel, 2, 1)
#        button_layout.addWidget(but_export, 3, 0)
        
        layout = pyqt5widget.QGridLayout(self.main_widget)
        layout.addWidget(self.plot_window, 0, 0)
        layout.addWidget(mpl_toolbar, 1, 0)
        layout.addWidget(slider_widgets, 0, 1)
        layout.addWidget(button_widget, 1, 1)
        
        layout.setColumnMinimumWidth(0, 1024)
        layout.setRowMinimumHeight(0, 768)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)

        self.statusBar().showMessage("Domi Slider Fit App v" + str(self.version))
    
    def slider_value_changed(self, value):
        changed_slider = self.sender()
        cur_param = self.plot_window.p[self.sliders_inverse[changed_slider]]
        
        minval = cur_param.min
        maxval = cur_param.max
        delta = (maxval - minval)/1000.
        new_value = minval + changed_slider.value()*delta
        if new_value > 1e3 or new_value < 1e-3:
            prec = '{:.3e}'
        else:
            prec = '{:.3f}'
        changed_slider.label.setText(prec.format(new_value))
        cur_param.value = new_value
        self.plot_window.update_plot()
        

    def closeEvent(self, event):
        self.fileQuit()

    def fileQuit(self):
        self.close()

    def about(self):
        pyqt5widget.QMessageBox.about(self, "About",
              """
              Slider - Fitting App
              Copyright 2016 Dominixue Dresen

              Version """+str(self.version)+"""

              Program to estimate fit parameters as initial step.
              """)
              
class cPlotAndFit(FigureCanvas):
    def __init__(self, parent=None):
        self.parent = parent
        
        self.p = None # has to be set in init_data !
        self.x = None # has to be set in init_data !
        
        self.y = None # can be set in init_data
        self.sy = None # can be set in init_data
        self.ymodel = None # has to be set in init_data !
        
        self.data_path = None
        
        self.fit_result = None
        self.init_data()
        self.fig = Figure(figsize=(4, 4))#, dpi=100)
        self.define_plot_canvas()
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
                pyqt5widget.QSizePolicy.Expanding,
                pyqt5widget.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        self.fig.tight_layout()
        
    def init_data(self):
        sys.exit("Define init_data() for cPlotAndFit. Setting initial "+\
                 "parameters p as well as x and ymodel. Optionally y and sy")

    def get_model(self, p, x):
        sys.exit("Define self.ymodel=get_model(p,x) for cPlotAndFit")

    def define_plot_canvas(self):
        self.ax1 = self.fig.add_subplot(111)
        
        self.ax1.set_xlabel("$\mathit{x}$")
        self.ax1.set_ylabel("$\mathit{y}$")
        if self.x is not None and self.y is not None and self.sy is not None:
            self.ax1.errorbar(self.x, self.y, self.sy, marker='.',\
                    linestyle='None', color='#4dac26', label=self.data_path)

        self.model_plot, = self.ax1.plot(self.x, self.ymodel, marker='None',\
                linestyle='-', color='#ca0020', lw=1, label="Model")
    
    def update_plot(self):
        self.ymodel = self.get_model(self.p, self.x)
        self.model_plot.set_ydata(self.ymodel)
        self.update_chi2()
        self.draw()

    def update_vary_vals_of_params(self):
        for parameter in self.p:
            self.p[parameter].vary =\
                     self.parent.checkboxes[parameter].isChecked()

    def figure_of_merit(self, p):
        self.ymodel = self.get_model(p, self.x)
        return (self.ymodel-self.y)/self.sy
        
    def update_parent_sliders(self):
        sliders = self.parent.sliders
        for parameter in self.p:
            cur_param = self.p[parameter]

            slider_bar = sliders[parameter]

            minval = cur_param.min
            maxval = cur_param.max
            delta = (maxval - minval)/1000.

            slider_value = int((cur_param.value-minval)/delta)
            slider_bar.setValue(slider_value)

    def fit_local(self):
        self.update_vary_vals_of_params()
        print("Running Levenberg-Marquardt.")
        self.fit_result = lmfit.minimize(self.figure_of_merit, self.p)
        self.p = self.fit_result.params
        print(lmfit.fit_report(self.fit_result))
        self.update_parent_sliders()
    
    def fit_local_wo_bounds(self):
        self.update_vary_vals_of_params()
        p_wo_bounds = lmfit.Parameters()
        for parameter in self.p:
            p_wo_bounds.add(parameter, self.p[parameter].value,\
                            vary=self.p[parameter].vary)
        print("Running Levenberg-Marquardt without bounds in parameters.")
        self.fit_result = lmfit.minimize(self.figure_of_merit, p_wo_bounds)
        p_wo_bounds = self.fit_result.params
        for parameter in self.p:
            self.p[parameter].value = p_wo_bounds[parameter].value
        print(lmfit.fit_report(self.fit_result))
        self.update_parent_sliders()

    
    def fit_global(self):
        self.update_vary_vals_of_params()
        print("Running Differential Evolution.")
        self.fit_result = lmfit.minimize(self.figure_of_merit, self.p,\
                method="differential_evolution")
        self.p = self.fit_result.params
        print(lmfit.fit_report(self.fit_result))
        self.update_parent_sliders()
            
    def export_model(self):
        modelfile = "modelfile.dat"
        savefile = open(modelfile, "w")
        if self.fit_result is not None:
            savefile.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))

        savefile.write("#Fitrange: " + str(self.qmin) +\
                " .. " + str(self.qmax) + "\n")
        savefile.write("#x \t y \t sy \t ymodel\n")
        for ix, xval in enumerate(self.qz):
            savefile.write(str(xval) +"\t"+ str(self.y[ix])+"\t"+\
                    str(self.sy[ix]) + "\t" + str(self.ymodel[ix])+"\n")
        print("Wrote results to " + modelfile)
        savefile.close()
    
    def save_para_to_script(self):
        self.update_vary_vals_of_params()
        script_file_name = sys.argv[0]
        script_file = open(script_file_name, "r")
        
        script_file_string = ""
        for line in script_file:
            if "self.p.add" in line:
                line_beginning = line.split('self.p.add')[0]
                parameter_name = line.split('self.p.add("')[-1].split('"')[0]
                line_end = line.split(")")[-1]
                para_value = self.p[parameter_name].value
                para_min = self.p[parameter_name].min
                para_max = self.p[parameter_name].max
                para_vary = self.p[parameter_name].vary
                script_file_string += line_beginning+\
                                      'self.p.add("'+parameter_name+'", '+\
                                       str(para_value) +", "+\
                                       "min = " + str(para_min) +", "+\
                                       "max = " + str(para_max) +", "+\
                                       "vary = " + str(para_vary) + ")"+\
                                       line_end
            else:
                script_file_string += line
        script_file.close()
        script_file = open(script_file_name, "w")
        script_file.write(script_file_string)
        script_file.close()
        print("Updated script file: " + script_file_name)
        
    def update_chi2(self):
        self.parent.label_chi2.setText("chi2/ndof: " +\
                             "{:.3f}".format(self.chi2))

    
#    
#    def export_params(self):
#        save_file = open("parameters.dat", "w")
#        
#        save_file.write("#Parameter estimated and exported using sas gui version " +\
#                 str(self.version) + "\n")
#        save_file.write("#Data file was: " + self.exp_data_path + "\n")
#        if self.fit_result is not None:
#            save_file.write("#"+lmfit.fit_report(self.fit_result).replace("\n", "\n#"))
#        save_file.write("\n\n")
#        for parameter in self.p:
#            save_file.write(parameter + "\t" + str(self.p[parameter].value) + "\n")
#        save_file.write("\n\n")
#        save_file.write(sas_methods.get_params_list(self.p))
#        save_file.close()
#        print("Saved parameters to parameters.dat")


if __name__ == '__main__':
    
    p = lmfit.Parameters()
    p.add("test_param1", 10, min=0, max=100, vary=1)
    p.add("test_param2", 50, vary=0)
    app = pyqt5widget.QApplication(sys.argv)
    aw = SliderFitApp(cPlotAndFit)
    aw.setWindowTitle("Slider Fitting Application")
    aw.show()
    app.exec_()
    
