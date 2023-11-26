# -*- coding: utf-8 -*-
"""
Created on Fri Feb 15 19:16:23 2019

@author: Thalpy
"""

import PySimpleGUIedit as sg
import pyAFM_FC.mfp_force_curves as mfp
import pyAFM_FC.proc_force_sep as sep

#import winsound
import random
import glob
import os

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasAgg
#import matplotlib.backends.tkagg as tkagg
import tkinter as Tk

"""
TODO:
    #port autofit
    #join parts together
    #batch processing
    #testing
    #import parameters
"""

#def _RunCode():
    
def draw_figure(canvas, figure, loc=(0, 0)): 
    """ Draw a matplotlib figure onto a Tk canvas
    loc: location of top-left corner of figure on canvas in pixels.
    Inspired by matplotlib source: lib/matplotlib/backends/backend_tkagg.py
    """
    figure_canvas_agg = FigureCanvasAgg(figure)
    figure_canvas_agg.draw()
    figure_x, figure_y, figure_w, figure_h = figure.bbox.bounds
    figure_w, figure_h = int(figure_w), int(figure_h)
    photo = Tk.PhotoImage(master=canvas, width=figure_w, height=figure_h)
    canvas.create_image(loc[0] + figure_w/2, loc[1] + figure_h/2, image=photo)
    tkagg.blit(photo, figure_canvas_agg.get_renderer()._renderer, colormode=2)
    return photo

def fcPlotWindow():
    """
    Part of the matplotlib plotting function
    """
    #figure_x, figure_y, figure_w, figure_h = figbox
    
    plotLayout = [[sg.Text('Plot test', font='Any 18', key="_canvNam")],
               [sg.Canvas(size=(1000, 400), key='canvas')]]
    
    plotWindow = sg.Window('MFP force curve plots', force_toplevel=True).Layout(plotLayout).Finalize()
    
    return plotWindow

class forceGUI(object):
    """
    The main window function, sets up UI and handles inputs.
    """
    
    def __init__(self):

        openBox = False
        test = "howdy doody"
        iTTip = "The base folder/file used for processing, for a batch file, select a folder containing all subfolders with files in."
        oTTip = "The base folder/file used for outputting data/graphs into, for a batch file, select a folder where all outputs will be put into"
        
        #sg.FolderBrowse(disabled=True, visible=False, key = "_inpFold")
        
        
        
        self.layout =    [[sg.Text('ProcForceSep Alpha v0.2', text_color="#0ea7ee")],
                       
                      [sg.Text('Script to run:', tooltip = iTTip), 
                       sg.InputCombo(('MFP Extract', 'Force Analysis', 'Both'), size=(12, 1), key="_scriptChoice", enable_events = True), 
                       sg.InputCombo(('Approach', 'Retract', 'Both'), size=(10, 1), key = "_curveChoice"),
                       sg.InputCombo(('.pdf', '.jpg', '.png'), size=(10, 1), key = "_filetype"),
                       sg.Button('Run', key ="_run"), sg.Button('?')],
                      
                      [sg.Text('Input file/path:', tooltip = iTTip, text_color="#f92672")],
                      
                      [sg.Text('Input file:', tooltip = iTTip, key = "_InputPath"), sg.Input(tooltip = iTTip, key = "_inpTxt", do_not_clear=True), sg.FileBrowse(key = "_inpFile"), sg.Checkbox('Batch job', enable_events=True, key = "_Batch")],
                      [sg.Text('Output path:', tooltip = oTTip, key = "_OutputPath"), sg.Input(disabled = True, tooltip = oTTip, key = "_OutTxt", do_not_clear=True), sg.FolderBrowse(disabled = True, key = "_OutInp"), sg.Checkbox('Enable', enable_events=True, tooltip = "This is for testing", key = "_OutCheck")],
                      #[sg.Text(text=str(openBox), tooltip = "This is for testing")],  
                      
                      #[sg.Button("Browse", key = "_browTest", button_type=21)],
                      
                      [sg.Text('Fitting parameters:', text_color="#f92672")],
                      
                      [sg.Input("Kc", tooltip = "Spring constant of the cantilever", do_not_clear=True, size=(12,1), key = "_k_c"),  
                      sg.Input("Fitbin", do_not_clear=True, size=(12,1), key = "_fitbin", disabled = True),
                      sg.Input("Threshhold", do_not_clear=True, size=(12,1), key = "_threshhold", disabled = True),
                      sg.Checkbox('Load from file', tooltip = "Load parameters from file", enable_events=True, key = "_LoadParam",size=(9,1), disabled = True),
                      sg.Checkbox('MFP iGraphs', enable_events=True, size=(9,1), key ="_Graphs"), 
                      sg.Checkbox('Force iGraphs', visible = False, disabled = True, key = "_EGraphs")],
                       
                      [sg.Input("dfit_min", do_not_clear=True, size=(12,1), key = "_dfit_min", disabled = True),
                      sg.Input("dfit_max", do_not_clear=True, size=(12,1), key = "_dfit_max", disabled = True),
                      sg.Input("cfit_min", do_not_clear=True, size=(12,1), key = "_cfix_min", disabled = True),
                      sg.Input("cfit_max", do_not_clear=True, size=(12,1), key = "_cfit_max", disabled = True)],
                       
                      [sg.Checkbox('Advanced', tooltip = "This is for testing", enable_events=True, key = "_Advanced",size=(8,1), disabled = True),
                      sg.Input("Binsize", do_not_clear=True, size=(12,1), key='_Binsize', tooltip = "Binsize I'M NOT SURE WHAT THIS DOES", disabled = True, visible=True),
                      sg.Input("Squarerange", do_not_clear=True, size=(12,1), key='_Squarerange', disabled = True, visible=True),
                      sg.Input("h", do_not_clear=True, size=(12,1), key='_h', disabled = True, visible=True)],
                       
                      [sg.Checkbox('Auto dfit'), sg.Checkbox('Sound', key = "_Sound", default = True), sg.Checkbox('Individual Graphs', enable_events=True, key ="_Graphs"), sg.Checkbox('Extra graphs', visible = False, disabled = True, key = "_EGraphs")],  
                      #[sg.Input(do_not_clear=True)],
                      
                      [sg.ProgressBar(100, orientation='h', size=(40, 20), key='_BatchBar'), sg.Text('jCount / TotBatch', key = "_BatchTxt")],
                      [sg.ProgressBar(100, orientation='h', size=(40, 20), key='_InduBar'), sg.Text('numCount / TotCurve', key = "_InduBar")],
                      
                      [sg.Output(size=(80,10))]]    
                  
        layoutInfo = [[sg.Text('Readme Here!')]]
        sg.Input()
        
        self.window = sg.Window('ProcForceSep Alpha v0.1').Layout(self.layout)      
        self.infoBox = sg.Window('Information popup').Layout(layoutInfo)   
        openBox = False
        sg.ChangeLookAndFeel('Dark')
        
        paramList = ["_k_c", "_fitbin", "_threshhold", "_dfit_min", "_dfit_max", "_cfix_min", "_cfit_max"]
        
       
        
        print("If you want to print to terminal, follow with a self.window.Refresh()")
        
        #Check for pygame and disable sound if not present.
        event, values = self.window.Read()  
        try:
            import pygame
            #for sfx:
            pygame.init()
            pygame.mixer.init()
            mysound = pygame.mixer.Sound("bleep_high.ogg")
            pygameP = True
            
        
        except:
            pygameP = False
            self.window.FindElement("_Sound").Update(disabled=True, value = False)
            print("Please install pygame for sound")
        self.window.Refresh()
        #self.window.FindElement("_Binsize").Update(disabled=True, visible=False)
        #self.window.FindElement("_Squarerange").Update(disabled=True, visible=False)
        #self.window.FindElement("_h").Update(disabled=True, visible=False)
        
        """
        PROGRESS BAR CODE:
        
        [sg.ProgressBar(100, orientation='h', size=(20, 20), key='progressbar')]
        progress_bar = self.window.FindElement('progressbar')   
        progress_bar.UpdateBar(i + 1)      
        """
        
        while True:      
            event, values = self.window.Read()     
            #print(event)
            
            if event is None or event == 'Exit':      
                break 
            
            if event == "_run":
                """
                Main running fuction of the code
                """
                #print(self.window.FindElement("_inpTxt").Get())
                inpTxtDir = self.window.FindElement("_inpTxt").Get()
                #k_c = self.window.FindElement("_Graphs").Get()
                
               # print(inpTxtDir)
                
                if self.window.FindElement("_Batch").Get() == True:
                    
                    namelist = sorted(glob.glob(inpTxtDir+'/**/*.txt', recursive=True))
                    OutTxtDir = self.window.FindElement("_OutTxt").Get()
                    
                else:
                    
                    namelist = inpTxtDir
                    OutTxtDir = None
                    
                #print(namelist)
                """
                MFPextract script handler
                """
                if values["_scriptChoice"] == "MFP Extract" or values["_scriptChoice"] == "Both":
                    print(values["_k_c"], type(values["_k_c"]), values["_filetype"], type(values["_filetype"]), values["_Graphs"], type(values["_Graphs"]))
                    try:
                        values["_k_c"] = float(values["_k_c"])
                    except:
                        print("Error converting Kc")
                    
                    if isinstance(values["_k_c"], str):
                        print("Please input Kc.")
                    
                    else:            
                        """
                        Graphs enabled depends on igraphs
                        """
                        for item in namelist:
                            #if "Def" in item:
                            print("Processing file:", item)
                            MFP = mfp.MFP_curves()
                            totFiles, outdir = MFP.setup(item, OutTxtDir)
                            self.window.FindElement('_BatchBar').Update(totFiles)
                            batch_bar = self.window.FindElement('_BatchBar')   
                            self.window.FindElement("_OutTxt").Update(outdir)
                            self.window.Refresh()
                            plotWindow = fcPlotWindow()
                            for i in range (0, totFiles):
                                batch_bar.UpdateBar(i+1) 
                                print("Processing file", i, "of", totFiles, end="\r")
                                self.window.FindElement('_BatchTxt').Update((str(i)+"/"+str(totFiles)))
                                #self, k_c, filetype, i, iGraphs = True)
                                try:
                                    fig = MFP.split_curves(values["_k_c"], values["_filetype"], i, values["_Graphs"])
                                except:
                                    print(item+" failed to be processed, skipping")
                                    raise
                                #fcPlotWindow(fig, figbox)
                                draw_figure(plotWindow.FindElement('canvas').TKCanvas, fig)
                                self.window.Refresh()
                                plotWindow.Refresh()
                    
        
            if event == "_scriptChoice":
                #if self.window.FindElement("_scriptChoice").Get() == "MFP Extract":
                if values["_scriptChoice"] == "MFP Extract":
                    for item in paramList[1:]:
                        self.window.FindElement(item).Update(disabled = True)
                        self.window.FindElement("_Advanced").Update(disabled = True, value = False)
                        self.window.FindElement("_Binsize").Update(disabled=True)
                        self.window.FindElement("_Squarerange").Update(disabled=True)
                        self.window.FindElement("_LoadParam").Update(disabled=True)
                        self.window.FindElement("_h").Update(disabled=True)
                        
                        
                        
                else:
                    for item in paramList[1:]:
                        self.window.FindElement(item).Update(disabled = False)
                        self.window.FindElement("_Advanced").Update(disabled = False)
                        self.window.FindElement("_LoadParam").Update(disabled=False)
                        
                        
            
            if event == "?" or openBox == True:
                self.infoBox.Show()
                eventI, valuesI = self.infoBox.Read()
                openBox = True
                if eventI is None or eventI == 'Exit': 
                    self.infoBox.Close()
                    openBox = False
            
            if event == "_Graphs":
                if self.window.FindElement("_Graphs").Get() == True:
                    self.window.FindElement("_EGraphs").Update(disabled = False, visible = True) 
                else:
                    self.window.FindElement("_EGraphs").Update(disabled = True, visible = False, value = False)
                    
                    
            if event == "_OutCheck":
                if self.window.FindElement("_OutCheck").Get() == True:
                    self.window.FindElement("_OutInp").Update(disabled = False)
                    self.window.FindElement("_OutTxt").Update(disabled = False)
                    
                else:
                    self.window.FindElement("_OutInp").Update(disabled = True)
                    self.window.FindElement("_OutTxt").Update(disabled = True)
            
            if event == "_Batch":
                if self.window.FindElement("_Batch").Get() == True:
                    self.window.FindElement("_InputPath").Update(value = "Input path:")
                    self.window.FindElement("_inpTxt").Update(value = "")
                    self.window.FindElement("_inpFile").Update(button_type=1)
                    #self.window.FindElement("_inpFold").Update(disabled = False, visible = True)
                    #self.window.FindElement("_inpFile").Update(disabled = True, visible = False, value = "")
                else:
                    self.window.FindElement("_InputPath").Update(value = "Input file:")
                    self.window.FindElement("_inpTxt").Update(value = "")
                    self.window.FindElement("_inpFile").Update(button_type=21)
                    #self.window.FindElement("_inpFold").Update(disabled = True, visible = False, value = "")
                    #self.window.FindElement("_inpFile").Update(disabled = False, visible = True)
                    
            if event == "_Advanced":
                if self.window.FindElement("_Advanced").Get() == True:
                    #print("FalseBranch")
                    self.window.FindElement("_Binsize").Update(disabled=False)
                    self.window.FindElement("_Squarerange").Update(disabled=False)
                    self.window.FindElement("_h").Update(disabled=False)
                elif self.window.FindElement("_Advanced").Get() == False:
                    #print("TrueBranch")
                    self.window.FindElement("_Binsize").Update(disabled=True)
                    self.window.FindElement("_Squarerange").Update(disabled=True)
                    self.window.FindElement("_h").Update(disabled=True)
                #print("test:", self.window.FindElement("_Advanced").Get())
                #print(self.window.FindElement("_h").FindElement(disabled))
                #print(self.window.FindElement(disabled))
            
            if self.window.FindElement("_Sound").Get() == True and pygameP == True:
                #winsound.Beep(random.randint(37,3767), 200)
                
                

                mysound.play()

                
            #print("Event:",event,"\n ValueList:", values) #'Kc', 'Binsize', 'Threshhold', 'h', 'dfit_min', 'dfit_max', 
            #'cfit_min', 'cfit_max', 'Fitbin', 'Squarerange', 'h', False, False, False, 'Approach'
            
            
        self.window.Close()
        
        #def proc_force_sep(name, ThalpyCount, k_c = 0.144, binsize = 0.5, thresh = 0, dfit_min = 1000, dfit_max = 2000,
        #    cfit_min = 38, cfit_max = 60, fitbin = 5, squareRange = 20, h=3.444102299, auto_dfit = False, extra = False, approach = True):
     
forceGUI()

"""
layout = [[sg.Text('All graphic widgets in one self.window!', size=(30, 1), font=("Helvetica", 25), text_color='blue')], #  
   [sg.Text('Here is some text.... and a place to enter text')],      
   [sg.InputText()],      
   [sg.Checkbox('My first checkbox!'), sg.Checkbox('My second checkbox!', default=True)],      
   [sg.Radio('My first Radio!     ', "RADIO1", default=True), sg.Radio('My second Radio!', "RADIO1")],      
   [sg.Multiline(default_text='This is the default Text shoulsd you decide not to type anything',)],      
[sg.InputCombo(['Combobox 1', 'Combobox 2'], size=(20, 3)),      
 sg.Slider(range=(1, 100), orientation='h', size=(35, 20), default_value=85)],      
[sg.Listbox(values=['Listbox 1', 'Listbox 2', 'Listbox 3'], size=(30, 6)),      
 sg.Slider(range=(1, 100), orientation='v', size=(10, 20), default_value=25),      
 sg.Slider(range=(1, 100), orientation='v', size=(10, 20), default_value=75),      
 sg.Slider(range=(1, 100), orientation='v', size=(10, 20), default_value=10)],      
[sg.Text('_'  * 100, size=(70, 1))],      
[sg.Text('Choose Source and Destination Folders', size=(35, 1))],      
[sg.Text('Source Folder', size=(15, 1), auto_size_text=False, justification='right'), sg.InputText('Source'),      
 sg.FolderBrowse()],      
[sg.Text('Destination Folder', size=(15, 1), auto_size_text=False, justification='right'), sg.InputText('Dest'),      
 sg.FolderBrowse()],      
[sg.Submit(), sg.Cancel(), sg.Button('Customized', button_color=('white', 'green'))]]

event, values  = sg.Window('Everything bagel', auto_size_text=True, default_element_size=(40, 1)).Layout(layout).Read()
"""
