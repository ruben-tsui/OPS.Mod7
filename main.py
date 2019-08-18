# -*- coding: utf-8 -*-

import time, datetime, os
#os.environ["QT_AUTO_SCREEN_SCALE_FACTOR"] = "1"

import sys, copy
#import design
#from PyQt5.QtCore import Qt
from PyQt5 import QtCore, QtGui, QtWidgets, uic
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import QDialog, QApplication, QMainWindow, QPushButton, QVBoxLayout, QHBoxLayout, QSlider, QLineEdit, QLabel, QTableWidgetItem, QSizePolicy, QTableWidgetSelectionRange
from matplotlib.widgets import Slider, Button, RadioButtons

# Libraries for drawing figures
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.ticker import FuncFormatter
# Ensure using PyQt5 backend
import matplotlib as mpl
mpl.use('QT5Agg')
import matplotlib.pyplot as plt
#plt.style.use('ggplot')
#plt.style.use('seaborn')
#plt.tight_layout(pad = 1.25)
#plt.style.use('dark_background')

# Libraries for stochastic simulation
import numpy as np
from numpy import log, exp, sqrt
from scipy.stats import skew, kurtosis
from math import ceil, floor
from scipy.special import binom
from scipy.stats import norm
import bsm_base as bsm

#from Simulation import merton_jump_diffusion_simulate as mjd
from Simulation import geometric_brownian_motion_simulate as gbm

class Window(QMainWindow):
    def __init__(self):
        QtWidgets.QDialog.__init__(self)
        # load the .ui file created in Qt Creator directly
        print(f"screen width = {screen_width}")
        if screen_width >= 1920:
            Interface = uic.loadUi('mainwindow7.ui', self)
        else:
            Interface = uic.loadUi('mod7.ui', self)

        # default values for all parameters
        self.S0 = self.slider_S0.value() #100
        # call options 1 - 3: strike prices
        self.K_C1 = self.slider_K_C1.value()
        # put options 1 - 3: strike prices
        self.K_P1 = self.slider_K_P1.value()
        # call options - times to maturity
        self.T_C1 = self.slider_T_C1.value()*.25
        # put options - times to maturity
        self.T_P1 = self.slider_T_P1.value()*.25
        # share parameters
        self.r   = self.slider_r.value()/100
        self.q   = self.slider_q.value()/100
        self.σ   = self.slider_sigma.value()/100
        self.μ   = self.slider_mu.value()/100
        # share simulation parameters
        self.Δt  = 0.001
        self.seed = 12345
        self.simSteps = self.spin_numSteps.value() # no. of steps in simulated share price changes; default 5
        self.simCurrentStep = 0
        self.simY = None
        self.Ydiscrete = None
        self.circleCoords = []
        self.PL_TVOM = [[None for col in range(4)] for row in range(self.simSteps+1)] # initialize matrix (6x4) 
        self.PL      = [[None for col in range(4)] for row in range(self.simSteps+1)] # same as above, but no TVOM 
        #params = {'S0': log(self.S0), 'K': log(self.K_C), 'r': self.r, 'q': self.q, 'σ': self.σ, 'T': self.T_C, 'seed': self.seed}
        #self.params = params
        self.NEWSIM = True # Boolean flag to indicate whether a new simulation is necessary 

        self.button_Simulate.hide()

        self.slider_S0.valueChanged.connect(self.on_change_S0)

        self.slider_K_C1.valueChanged.connect(self.on_change_K_C1)
        self.slider_K_P1.valueChanged.connect(self.on_change_K_P1)
        self.slider_T_C1.valueChanged.connect(self.on_change_T_C1)
        self.slider_T_P1.valueChanged.connect(self.on_change_T_P1)

        self.slider_r.valueChanged.connect(self.on_change_r)
        self.slider_q.valueChanged.connect(self.on_change_q)
        self.slider_sigma.valueChanged.connect(self.on_change_sigma)
        self.slider_mu.valueChanged.connect(self.on_change_mu)

        self.button_Simulate.clicked.connect(self.on_click_ShowSimulation)
        self.button_NextStep.clicked.connect(self.on_click_NextStep)
        self.button_NewSim.clicked.connect(self.on_click_NewSim)

        self.spin_numSteps.valueChanged.connect(self.on_valueChange_numSteps)
        self.spin_numSteps.visible = False

        self.radio_ignore_TVOM.clicked.connect(self.on_click_ignore_TVOM)

        self.radio_Profit.clicked.connect(self.on_click_Profit)
        self.radio_Payoff.clicked.connect(self.on_click_Payoff)

        # initialize default parameter values

        # Run some initialization stuff
        # simulate a batch of prices

        # Initialize all parameter labels 
        self.label_S0.setText(f"S0 = {self.S0}") 
        self.label_K_C1.setText(f"K = {self.K_C1}") 
        self.label_K_P1.setText(f"K = {self.K_P1}")   
        self.label_T_C1.setText(f"T = {self.T_C1}")   
        self.label_T_P1.setText(f"T = {self.T_P1}")   
        self.label_q.setText(f"q = {self.q}")     
        self.label_r.setText(f"r = {self.r}")      
        self.label_sigma.setText(f"σ = {self.σ}") 
        self.label_mu.setText(f"μ = {self.μ}") 

        self.radio_CoveredCall.clicked.connect(self.on_click_CoveredCall)
        self.radio_CoveredCallReverse.clicked.connect(self.on_click_CoveredCallReverse)
        self.radio_ProtectivePut.clicked.connect(self.on_click_ProtectivePut)
        self.radio_ProtectivePutReverse.clicked.connect(self.on_click_ProtectivePutReverse)

        self.setGeometry(0, 0, 1366, 768)

        # initialize default strategy (Covered Call)
        self.on_click_CoveredCall()
        self.on_click_NewSim()

    def on_click_Payoff(self):
        self.showGraph()

    def on_click_Profit(self):
        self.showGraph()

    def on_click_CoveredCall(self):
        self.spin_C1.setValue(-1)
        self.spin_Share.setValue(1)
        self.spin_P1.setValue(0)
        self.on_click_NewSim()

    def on_click_CoveredCallReverse(self):
        self.spin_C1.setValue(1)
        self.spin_Share.setValue(-1)
        self.spin_P1.setValue(0)

    def on_click_ProtectivePut(self):
        self.spin_P1.setValue(1)
        self.spin_Share.setValue(1)
        self.spin_C1.setValue(0)
        self.on_click_NewSim()

    def on_click_ProtectivePutReverse(self):
        self.spin_P1.setValue(-1)
        self.spin_Share.setValue(-1)
        self.spin_C1.setValue(0)

    ## The following 2 functions are for selecting between the Call and Put option in the Groupbox
    def on_click_Call(self):
        if self.radio_Call.isChecked():
            self.groupBox_Put.setGeometry(1340, 270, 271, 221)  # reset Put groupbox's position (off-screen)
            self.groupBox_Call.setGeometry(830, 270, 271, 221)  # move Call groupbox to "on-screen" position
        self.on_click_NewSim()
        #self.showGraph()

    def on_click_Put(self):
        if self.radio_Put.isChecked():
            self.groupBox_Call.setGeometry(1340, 20, 271, 221)  # reset Call groupbox's position (off-screen)
            self.groupBox_Put.setGeometry(830, 270, 271, 221)   # move Put groupbox to "on-screen" position
        self.on_click_NewSim()
        #self.showGraph()
        
    # The following is run when the "Ignore time value of money" radio button is clicked
    def on_click_ignore_TVOM(self):
        if self.radio_ignore_TVOM.isChecked():
            #print(f"ignore checked!")
            self.PL_TVOM = copy.deepcopy(self.PL)
        else:
            #print(f"ignore not chekced!")
            pass
        self.UpdatePL_TVOM()
        #self.DisplayPL_Table()
        self.showGraph()

    # Show the Excel-like table
    def DisplayPL_Table(self, hori_labels):

        # Set up the Table Widget (QTableView) environment - headers, widths, etc.
        if self.simSteps <= 5:
            self.tablePL.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.tablePL.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.tablePL.setColumnCount(4)
        self.tablePL.setRowCount(self.simSteps+1) # was 6
        #hori_labels = ['1P Value', '1C Value', 'xC+yP Value', 'P/L']
        #hori_labels = ['1P\nValue', '1C\nValue', 'ΣC+ΣP\nValue', 'P/L']
        self.tablePL.setHorizontalHeaderLabels(hori_labels)
        font = QFont('Segoe UI', 9)
        font.setBold(True)
        self.tablePL.horizontalHeader().setFont(font)
        vert_labels = [f"t{str(j)}" for j in range(self.simSteps+1)]
        self.tablePL.setVerticalHeaderLabels(vert_labels)
        font = QFont('Segoe UI', 9)
        font.setBold(False)
        self.tablePL.verticalHeader().setFont(font)
        # Headers
        self.tablePL.resizeColumnsToContents()
        header = self.tablePL.horizontalHeader()
        header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeToContents)
        for col in range(4):    
            header.setSectionResizeMode(col, QtWidgets.QHeaderView.Stretch)

        self.tablePL.clearSelection()

        # copy from self.PL_TVOM to self.tablePL
        k = self.simCurrentStep
        for row in range(k+1):
            for col in range(4):
                item = QTableWidgetItem(f"{self.PL_TVOM[row][col]:.2f}")
                item.setTextAlignment(Qt.AlignRight)
                self.tablePL.setItem(row, col, item)
        

    def on_valueChange_numSteps(self, value):
        self.simSteps = value

    def on_click_NextStep(self, value):
        self.tablePL.clear()
        #print(f"length of circleCoords = {len(self.circleCoords)}")
        self.simCurrentStep += 1
        if self.simCurrentStep > self.simSteps:
            self.simCurrentStep = 0
            #self.circleCoords = []
        #print(f"Current step = {self.simCurrentStep}")
        self.showGraph()

    def on_change_S0(self, value):
        self.S0 = value
        self.label_S0.setText(f"S0 = {str(self.S0)}")
        self.on_click_NewSim()

    def on_change_K_C1(self, value):
        self.K_C1 = value
        self.label_K_C1.setText(f"K = {str(self.K_C1)}")
        self.on_click_NewSim()
        
    def on_change_K_P1(self, value):
        self.K_P1 = value
        self.label_K_P1.setText(f"K = {str(self.K_P1)}")
        self.on_click_NewSim()
        
    def on_change_T_C1(self, value): 
        self.T_C1 = value*.25
        self.label_T_C1.setText(f"T = {str(self.T_C1)}")
        self.on_click_NewSim()
        
    def on_change_T_P1(self, value): 
        self.T_P1 = value*.25
        self.label_T_P1.setText(f"T = {str(self.T_P1)}")
        self.on_click_NewSim()
        
    def on_change_r(self, value):
        self.r = value/100
        self.label_r.setText(f"r = {str(self.r)}")
        self.on_click_NewSim()
        
    def on_change_q(self, value):
        self.q = value/100
        self.label_q.setText(f"q = {str(self.q)}")
        self.on_click_NewSim()
        
    def on_change_sigma(self, value):
        self.σ = value/100
        self.label_sigma.setText(f"σ = {str(self.σ)}")
        self.on_click_NewSim()
        #self.recomputeCircle()
        
    def on_change_mu(self, value):
        self.μ = value/100
        self.label_mu.setText(f"μ = {str(self.μ)}")
        self.showGraph()

    def on_click_ShowSimulation(self):
        #print("Simulation button clicked")
        self.showGraph()

    def on_click_NewSim(self):
        if not self.NEWSIM: return None
        params = {'S0': self.S0, 'μ': self.μ, 'σ': self.σ, 'Δt': self.Δt, \
                  'T': self.T_C1, \
                  'N': 1, 'P': 1, 'seed': np.random.randint(32768)}
        self.simY = gbm(params)
        self.simCurrentStep = 0
        self.circleCoords = []
        Max_T = self.T_C1
        self.X = np.linspace(0.0, Max_T, int(1/self.Δt)+1)
        self.simSteps = self.spin_numSteps.value()

        self.PL_TVOM = [[None for col in range(4)] for row in range(self.simSteps+1)] # initialize matrix (6x4) 
        self.PL      = [[None for col in range(4)] for row in range(self.simSteps+1)] # same as above, but no TVOM 

        self.Xdiscrete = np.linspace(0.0, Max_T, self.simSteps+1)      ### CHANGES_REQUIRED
        self.Ydiscrete = [self.simY[int((len(self.X)-1)/self.simSteps * i)] for i in range(self.simSteps+1) ] # sample of simY at discrete time points
        S0   = self.S0
        K_C1  = self.K_C1
        K_P1  = self.K_P1
        r    = self.r
        q    = self.q
        σ    = self.σ
        T_C1  = self.T_C1
        T_P1  = self.T_P1
        #seed = np.random.randint(0,65535)
        # CALL
        α1 = self.spin_C1.value()
        # PUT
        β1 = self.spin_P1.value()
        # SHARE
        Share = self.spin_Share.value()

        for i in range(0, self.simSteps+1):
            if i == 0:
                current_spot = S0
            else:
                current_spot = self.Ydiscrete[i][0]
            spot_optval = α1 * bsm.bs_call_price(current_spot, K_C1, r, q, σ, T_C1*(1-(i)/(self.simSteps))) + \
                          β1 * bsm.bs_put_price(current_spot, K_P1, r, q, σ, T_P1*(1-(i)/(self.simSteps))) + \
                          Share * (current_spot - self.S0)
            self.circleCoords.append((current_spot, spot_optval))
        
        # clear P/L table widget
        self.tablePL.clear()
        self.showGraph()


    def UpdatePL_TVOM(self):
        '''
        This method recomputes the P/L table based on the current time step
        Results (no explicit output, only self.PL_TVOM is updated):
            An updated PL_TVOM matrix up to the current step as a side effect
        '''
        #print(f"entering UpdatePL_TVOM...")
        k = self.simCurrentStep
        r = self.r # current interest rate
        deltaT = self.T_C1 / self.simSteps
        self.PL_TVOM = copy.deepcopy(self.PL)
        if not self.radio_ignore_TVOM.isChecked(): # TVOM
            for row in range(k+1):
                for col in range(3):
                    self.PL_TVOM[row][col] = self.PL[row][col] * np.exp(r*(k - row)*deltaT)
            self.PL_TVOM[0][3] = 0.0   
            for row in range(1, k+1):
                # P/L column
                self.PL_TVOM[row][3] = self.PL_TVOM[row][2] - self.PL_TVOM[0][2]

    def showGraph(self):

        mpl = self.oMplCanvas.canvas
        [s1, s2, s3, s4, s5] = mpl.axes

        #print(f"showGraph::very:beginning:PL={self.PL}")

        S0   = self.S0
        K_C1 = self.K_C1
        K_P1 = self.K_P1
        r    = self.r
        q    = self.q
        σ    = self.σ
        T_C1 = self.T_C1
        T_P1 = self.T_P1

        k = self.simCurrentStep  # this is used a lot below 

        # Get the x, y units of Call and Put options, respectively
        # CALL
        α1 = self.spin_C1.value()
        # PUT
        β1 = self.spin_P1.value()
        #print(f"α={α}, β={β}")
        # SHARE
        Share = self.spin_Share.value()


        # Set minimum and maximum x-axis values for s1 (main plot)
        # get the y-coordinates of self.circleCoords (options values) for t0, t1, ..., t5
        circle_y_list = [self.circleCoords[i][0] for i in range(len(self.circleCoords))]
        s1_x_max = max(140, max(circle_y_list)+10)
        s1_x_min = min(60, min(circle_y_list)-10)

        # x-axis points for s1 (main) plot
        x      = np.linspace(s1_x_min, s1_x_max, 51)
        self.x = x
        s      = x # stock value per share
        # y-axis points for s3 (1-unit call) plot
        y_call1 = bsm.bs_call_price(x, K_C1, r, q, σ, T_C1*(1-k/(self.simSteps)))
        # y-axis points for s2 (1-unit put) plot
        y_put1  = bsm.bs_put_price(x, K_P1, r, q, σ, T_P1*(1-k/(self.simSteps)))
        # y-axis points for s1 (combined options)
        y      = α1 * y_call1 + β1 * y_put1 + Share * s
        # payoff for 1-unit call
        z_call1 = bsm.payoff_call(x, K_C1)
        # payoff for 1-unit put
        z_put1  = bsm.payoff_put(x, K_P1)
        # payoff for combined
        z       = α1 * z_call1 + β1 * z_put1 + Share * s 

        profit_option = 0
        if α1:
            # profit = payoff + call premium (short call)
            profit_option = α1 * z_call1 + (-α1) * bsm.bs_call_price(S0, K_C1, r, q, σ, T_C1)
        if β1:
            # profit = payoff - put premium (long put)
            profit_option = β1 * z_put1 + (-β1) * bsm.bs_put_price(S0, K_P1, r, q, σ, T_P1)
        profit_stock  = Share * (x - S0)
        profit_combined = profit_stock + profit_option

        # s2: 1-unit Put/Call options
        s2.clear()
        if α1: # non-zero
            #title = f'$C_{1}$'
            title = f'{α1}$C$'
            s2.plot(x, α1*y_call1, 'r', lw=0.6, label='value')
            s2.plot(x, α1*z_call1, 'b-.', lw=1, label='payoff')
            s2.set_title(title, fontsize=12, color='brown')
            # plot circle on s2
            put_circle_x = self.circleCoords[k][0]
            put_circle_y = α1*bsm.bs_put_price(put_circle_x, K_C1, r, q, σ, T_C1*(1-k/(self.simSteps)))
            #s2.plot(put_circle_x, put_circle_y, 'black', marker="o")
        #self.text_PutValue1Unit.setPlainText(f"{put_circle_y:.2f}")
        s2.grid(True)

        # s3: 1-unit Put/Call options
        s3.clear()
        if β1: # non-zero
            title = f'{β1}$P_{1}$'
            s3.plot(x, β1*y_put1, 'r', lw=0.6, label='value')
            s3.plot(x, β1*z_put1, 'b-.', lw=1, label='payoff')
            s3.set_title(title, fontsize=12, color='brown')
            put_circle_x = self.circleCoords[k][0]
            put_circle_y = β1*bsm.bs_put_price(put_circle_x, K_P1, r, q, σ, T_P1*(1-k/(self.simSteps)))
            s3.plot(put_circle_x, put_circle_y, 'black', marker="o")
        s3.grid(True)
        # plot circle on s3
        call_circle_x = self.circleCoords[k][0]
        call_circle_y = bsm.bs_call_price(call_circle_x, K_C1, r, q, σ, T_C1*(1-k/(self.simSteps)))
        #s3.plot(call_circle_x, call_circle_y, 'black', marker="o")
        #self.text_CallValue1Unit.setPlainText(f"{call_circle_y:.2f}")

        # Plot all (discretized continuous) simulated share prices
        s4.clear()
        s4.set_title(f'Simulated share prices', fontsize=12, color='brown')
        s4.yaxis.set_tick_params(labelright=False, labelleft=True)
        s4.grid(True)
        #print(f"size of X = {len(self.X)}")
        # set the upper limit of index for which the simulated share prices will be shown on s4
        upperLimit = int((len(self.X)-1) * k / self.simSteps)
        #print(f"upperLimit = {upperLimit}")
        s4.plot(self.X[:upperLimit], self.simY[:upperLimit], 'g', lw=0.75, label='$S_T$')
        s4.set_xlim(right=self.T_C1, left=0.0)
        s4.set_ylim(top=self.simY.max()+10, bottom=self.simY.min()-10)
        # plot simulated prices at discrete time points (t0, t1, ..., t5)
        s4.plot(self.Xdiscrete[:k+1], self.Ydiscrete[:k+1], 'black', lw=1.25, label='price')
        #legend = s4.legend(loc='upper left', shadow=False, fontsize='medium')
        #print(f"Ydiscrete = {self.Ydiscrete}")

        # Plot main graph (s1: top row, middle)
        s1.clear()
        #if int(α) - float(α) == 0.0: # α is an integer
        main_title = '' # was f'${α1}C + {β1}P$'
        if Share > 0:
            main_title = f"${Share}S$"
        elif Share < 0:
            main_title = f"$-{Share}S$"
        weights = [α1, β1]
        for j, w in enumerate(weights[:1]): # α
            if w > 0:
                #main_title += f"$+{w} C_{j+1}$ "
                main_title += f"$+{w} C$ "
            if w < 0:
                #main_title += f"${w} C_{j+1}$ "
                main_title += f"${w} C$ "
        for j, w in enumerate(weights[1:]): # β
            if w > 0:
                #main_title += f"$+{w} P_{j+1}$ "
                main_title += f"$+{w} P$ "
            if w < 0:
                #main_title += f"${w} P_{j+1}$ "
                main_title += f"${w} P$ "

        s1.set_title(main_title, fontsize=12, color='brown')
        s1.set_aspect(aspect='equal')
        s1.yaxis.set_tick_params(labelright=False, labelleft=True)
        s1.grid(True)

        # Do we want to display Profit or Payoff?
        payoff_selected = True if self.radio_Payoff.isChecked() else False
        profit_selected = True if self.radio_Profit.isChecked() else False

        if payoff_selected:
            # payoff 
            if α1:
                s1.plot(x, z, 'r', lw=1.5, label='Covered Call')
            if β1:
                s1.plot(x, z, 'r', lw=1.5, label='Protective Put')
            s1.set_xlim(right=s1_x_max, left=s1_x_min)
            # long stock line
            if α1:
                #s1.plot(x, Share * (x - self.S0), 'g--', lw=1.5, label='Long Stock')
                s1.plot(x, Share * (x - 0), 'g--', lw=1.5, label='Long Stock')
            if β1:
                #s1.plot(x, Share * (x - self.S0), 'g--', lw=1.5, label='Long Stock')
                s1.plot(x, Share * (x - 0), 'g--', lw=1.5, label='Long Stock')
            # call/put
            if α1:
                #s1.plot(x, α1 * (z_call1), 'b-.', lw=1, label='Short Call')
                s1.plot(x, α1 * z_call1 + β1 * z_put1 + S0, 'b-.', lw=1, label='Short Call')
            if β1:
                #s1.plot(x, β1 * (z_put1), 'b-.', lw=1, label='Long Put')
                s1.plot(x, α1 * z_call1 + β1 * z_put1 + S0, 'b-.', lw=1, label='Long Put')
            s1.legend(loc='upper left', shadow=True, fontsize='large')

        if profit_selected:
            # payoff 
            if α1:
                s1.plot(x, profit_combined, 'r', lw=1.5, label='Covered Call')
            if β1:
                s1.plot(x, profit_combined, 'r', lw=1.5, label='Protective Put')
            s1.set_xlim(right=s1_x_max, left=s1_x_min)
            # long stock line
            if α1:
                s1.plot(x, profit_stock, 'g--', lw=1.5, label='Long Stock')
            if β1:
                s1.plot(x, profit_stock, 'g--', lw=1.5, label='Long Stock')
            # call/put
            if α1:
                s1.plot(x, profit_option, 'b-.', lw=1, label='Short Call')
            if β1:
                s1.plot(x, profit_option, 'b-.', lw=1, label='Long Put')
            s1.legend(loc='upper left', shadow=True, fontsize='large')

        # option premium (1 unit)
        if profit_selected:
            premium_call = bsm.bs_call_price(S0, K_C1, r, q, σ, T_C1)
            premium_put  = bsm.bs_put_price(S0, K_P1, r, q, σ, T_P1)
        if payoff_selected:
            premium_call = 0
            premium_put  = 0
        premium_offset = 0 if payoff_selected else (-α1) * premium_call + (-β1) * premium_put

        # plot a number of option value curves (from 0 to 5)
        for i in range(0, k+1):
            y_call1 = 0.0; y_put1 = 0.0
            if α1:
                y_call1 = bsm.bs_call_price(x, K_C1, r, q, σ, T_C1*(1-i/(self.simSteps)))
            # y-axis points for s2 (1-unit put) plot
            if β1:
                y_put1  = bsm.bs_put_price(x, K_P1, r, q, σ, T_P1*(1-i/(self.simSteps)))
            # y-axis points for s1 (combined options)
            if payoff_selected:
                y       = α1 * y_call1 + β1 * y_put1 + Share * (x - 0)
            elif profit_selected:
                y       = α1 * y_call1 + β1 * y_put1 + Share * (x - self.S0)
            #y      = α1 * y_call + β1 * y_put
            # plot option value curve
            if α1:
                s1.plot(x, y + (-α1) * premium_call, \
                    'red', alpha=0.6, lw=0.4+0.05*i, label=f'value = {put_circle_y + call_circle_y}')
            if β1:
                s1.plot(x, y + (-β1) * premium_put, \
                    'red', alpha=0.6, lw=0.4+0.05*i, label=f'value = {put_circle_y + call_circle_y}')
            #s1.legend(loc='upper left', shadow=False, fontsize='medium')
            
        # Now plot the circles
        #self.recomputeCircle()
        
        for i, point in enumerate(self.circleCoords):
            if i == 0:
                if payoff_selected:
                    s1.plot(point[0], point[1] + premium_offset + S0, 'black', marker="o")
                elif profit_selected:
                    s1.plot(point[0], point[1] + premium_offset, 'black', marker="o")
            elif i > 0 and i <= k: # skip step 0, i.e. when there's only one curve
                #print(f"i={i}, {self.circleCoords}")
                p_x, p_y = self.circleCoords[i-1] # previous circle's coordinates
                c_x, c_y = self.circleCoords[i]   # current circle's coordinates
                #print(f"p_x = {p_x:.2f}, p_y = {p_y:.2f}, c_x = {c_x:.2f}, c_y = {c_y:.2f}")
                if payoff_selected:
                    s1.plot(point[0], point[1] + premium_offset + S0, 'black', marker="o")
                    s1.arrow(p_x, p_y + premium_offset + S0, c_x - p_x, c_y - p_y, \
                        color='purple', width=0.00025, head_width=1.0, length_includes_head=True)
                elif profit_selected:
                    s1.plot(point[0], point[1] + premium_offset, 'black', marker="o")
                    s1.arrow(p_x, p_y + premium_offset, c_x - p_x, c_y - p_y, \
                        color='purple', width=0.00025, head_width=1.0, length_includes_head=True)


        # Plot P&L graph
        s5.clear()
        s5.set_title(f'Profit/Loss', fontsize=12, color='brown')
        s5.yaxis.set_tick_params(labelright=False, labelleft=True)
        s5.grid(True)

        # default horizontal header labels
        hori_labels = ['1P\nProfit', '1C\nProfit', 'ΣC+ΣP\nProfit', 'P/L']

        number_basis = 'Profit' if profit_selected else 'Value'
        # 1C vlaue
        stepSize = int((1/self.Δt)/self.simSteps)
        if α1:
            #self.PL[k][0] = put_circle_y
            self.PL[k][0] = (-α1)*bsm.bs_call_price(self.simY[k*stepSize][0], K_C1, r, q, σ, T_C1*(1-k/(self.simSteps)))
            hori_labels[0] = f"{α1}C\n{number_basis}"
        # 1P value
        if β1:
            self.PL[k][0] = call_circle_y
            self.PL[k][0] = (-β1)*bsm.bs_put_price(self.simY[k*stepSize][0], K_P1, r, q, σ, T_P1*(1-k/(self.simSteps)))
            hori_labels[0] = f"{β1}P\n{number_basis}"
        # zS value

        self.PL[k][1] = self.simY[k*stepSize][0] - S0 if profit_selected else self.simY[k*stepSize][0]
        hori_labels[1] = f"{Share}S\n{number_basis}"
        # xC|yP + zS
        self.PL[k][2] = (-α1) * self.PL[k][0] + β1 * self.PL[k][0] + Share*self.PL[k][1]
        #hori_labels[2] = f"{hori_labels[0]} + {hori_labels[1]}"
        hori_labels[2] = f"O+S\n{number_basis}"

        hori_labels[3] = "P/L" if profit_selected else "Value"
        # P/L
        #print(f"k={k}")
        self.PL[k][3] = self.PL[k][2] - self.PL[0][2]
        self.UpdatePL_TVOM()
        # copy from self.PL_TVOM to self.tablePL

        self.DisplayPL_Table(hori_labels)
                
        #item.setTextAlignment(Qt.AlignVCenter)
        # Select a row
        self.tablePL.setRangeSelected(QTableWidgetSelectionRange(k, 0, k, 3), True)
        for row in range(self.simSteps+1):
            self.tablePL.setRowHeight(row, 20)

        # Finally, plot the P/L line
        # find min and max of option values    
        all_options_values = [self.circleCoords[i][1] for i in range(len(self.circleCoords))]
        #s5_x_max = max(all_options_values)+10
        #s5_x_min = min(all_options_values)-10

        xPL = self.Xdiscrete
        option_val_t0 = self.PL_TVOM[0][2]  # option value at t0
        yPL = [self.PL_TVOM[row][2] for row in range(k+1)]
        #s5.set_ylim(top=s5_x_max, bottom=s5_x_min)
        s5.plot(xPL[:k+1], yPL, 'r-.', lw=1.5, label='P&L')
        # plot horizontal line indicating the option value at t0
        s5.axhline(y=option_val_t0, color="black", alpha=1, linewidth=1)
        # annotate plot with P/L values 
        for row in range(0, k+1):
            s5.annotate(f"{self.PL_TVOM[row][3]:+.2f}", xy=[self.Xdiscrete[row]-0.1, self.PL_TVOM[row][3]+option_val_t0+2])
        for row in range(k+1):
            s5.plot(self.Xdiscrete[row], self.PL_TVOM[row][3]+option_val_t0, marker='o', color='black')

        #self.tablePL.setStyleSheet("QTableView::item:selected { color:red; background:yellow; font-weight: bold; } " + "QTableView::vertical:header {font-size: 9}")


        # Display entire canvas!
        mpl.fig.set_visible(True)
        mpl.draw()

 




# This creates the GUI window:
if __name__ == '__main__':

    import sys
    app = QtWidgets.QApplication(sys.argv)
    ### load appropriately sized UI based on screen resolution detected
    screen = app.primaryScreen()
    screen_width=screen.size().width()
    window = Window()
    w = window
    t = w.tablePL
    #window.showGraph()
    window.show()
    sys.exit(app.exec_())

