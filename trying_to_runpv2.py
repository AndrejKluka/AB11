# -*- coding: utf-8 -*-
"""
Created on Mon May 14 16:28:12 2018

@author: timvd
"""

import time
from paraview.simple import *

frames=1
frame_names=[]
if frames>0:
    for i in range(frames):
        frame_names.append('uvwp_00{:03}.h5' .format(i+1))

times=1   
if frames!=1:
    times=frames

for frame in range(times) : 
    #import name
    ##################################
    colouring = 'Vorticity z'#raw_input("What type of colouring to you want? Normal, Vorticityx, Vorticityy of Vorticityz? ")
    Isosurface = .2 #float(raw_input("What is the isosurface value?"))
    ##################################
    
    
    
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()
    
    # create a new 'XML Rectilinear Grid Reader'
    validation_Q_l295of96Qvtr = XMLRectilinearGridReader(FileName=["C:\\Users\\timvd\\Documents\\GitHub\\AB11\\calculated data\\validation_Q_l2-Q.vtr"])
    validation_Q_l295of96Qvtr.PointArrayStatus = ['Q', 'Vorticity normal', 'Vorticity x', 'Vorticity y', 'Vorticity z']
    
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1168, 772]
    
    # show data in view
    validation_Q_l295of96QvtrDisplay = Show(validation_Q_l295of96Qvtr, renderView1)
    
    
    # trace defaults for the display properties.
    validation_Q_l295of96QvtrDisplay.Representation = 'Outline'
    validation_Q_l295of96QvtrDisplay.ColorArrayName = [None, '']
    validation_Q_l295of96QvtrDisplay.OSPRayScaleArray = 'Q'
    validation_Q_l295of96QvtrDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    validation_Q_l295of96QvtrDisplay.SelectOrientationVectors = 'None'
    validation_Q_l295of96QvtrDisplay.ScaleFactor = 9.4
    validation_Q_l295of96QvtrDisplay.SelectScaleArray = 'None'
    validation_Q_l295of96QvtrDisplay.GlyphType = 'Arrow'
    validation_Q_l295of96QvtrDisplay.GlyphTableIndexArray = 'None'
    validation_Q_l295of96QvtrDisplay.GaussianRadius = 0.47000000000000003
    validation_Q_l295of96QvtrDisplay.SetScaleArray = ['POINTS', 'Q']
    validation_Q_l295of96QvtrDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    validation_Q_l295of96QvtrDisplay.OpacityArray = ['POINTS', 'Q']
    validation_Q_l295of96QvtrDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    validation_Q_l295of96QvtrDisplay.DataAxesGrid = 'GridAxesRepresentation'
    validation_Q_l295of96QvtrDisplay.SelectionCellLabelFontFile = ''
    validation_Q_l295of96QvtrDisplay.SelectionPointLabelFontFile = ''
    validation_Q_l295of96QvtrDisplay.PolarAxes = 'PolarAxesRepresentation'
    
    
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    validation_Q_l295of96QvtrDisplay.ScaleTransferFunction.Points = [-2.5586026646007163, 0.0, 0.5, 0.0, 19.220113274066406, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    validation_Q_l295of96QvtrDisplay.OpacityTransferFunction.Points = [-2.5586026646007163, 0.0, 0.5, 0.0, 19.220113274066406, 1.0, 0.5, 0.0]
    
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    validation_Q_l295of96QvtrDisplay.DataAxesGrid.XTitleFontFile = ''
    validation_Q_l295of96QvtrDisplay.DataAxesGrid.YTitleFontFile = ''
    validation_Q_l295of96QvtrDisplay.DataAxesGrid.ZTitleFontFile = ''
    validation_Q_l295of96QvtrDisplay.DataAxesGrid.XLabelFontFile = ''
    validation_Q_l295of96QvtrDisplay.DataAxesGrid.YLabelFontFile = ''
    validation_Q_l295of96QvtrDisplay.DataAxesGrid.ZLabelFontFile = ''
    
    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    validation_Q_l295of96QvtrDisplay.PolarAxes.PolarAxisTitleFontFile = ''
    validation_Q_l295of96QvtrDisplay.PolarAxes.PolarAxisLabelFontFile = ''
    validation_Q_l295of96QvtrDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
    validation_Q_l295of96QvtrDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # create a new 'Contour'
    contour1 = Contour(Input=validation_Q_l295of96Qvtr)
    contour1.ContourBy = ['POINTS', 'Q']
    contour1.Isosurfaces = [Isosurface]
    contour1.PointMergeMethod = 'Uniform Binning'
    
    # Properties modified on contour1
    #contour1.Isosurfaces = [Isosurface]
    
    # show data in view
    contour1Display = Show(contour1, renderView1)
    
    # trace defaults for the display properties.
    contour1Display.Representation = 'Surface'
    contour1Display.ColorArrayName = [None, '']
    contour1Display.OSPRayScaleArray = 'Normals'
    contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour1Display.SelectOrientationVectors = 'None'
    contour1Display.ScaleFactor = 0.9620841979980469
    contour1Display.SelectScaleArray = 'None'
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = 'None'
    contour1Display.GaussianRadius = 0.04810420989990234
    contour1Display.SetScaleArray = ['POINTS', 'Normals']
    contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour1Display.OpacityArray = ['POINTS', 'Normals']
    contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour1Display.DataAxesGrid = 'GridAxesRepresentation'
    contour1Display.SelectionCellLabelFontFile = ''
    contour1Display.SelectionPointLabelFontFile = ''
    contour1Display.PolarAxes = 'PolarAxesRepresentation'
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour1Display.ScaleTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [-1.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
    
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    contour1Display.DataAxesGrid.XTitleFontFile = ''
    contour1Display.DataAxesGrid.YTitleFontFile = ''
    contour1Display.DataAxesGrid.ZTitleFontFile = ''
    contour1Display.DataAxesGrid.XLabelFontFile = ''
    contour1Display.DataAxesGrid.YLabelFontFile = ''
    contour1Display.DataAxesGrid.ZLabelFontFile = ''
    
    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    contour1Display.PolarAxes.PolarAxisTitleFontFile = ''
    contour1Display.PolarAxes.PolarAxisLabelFontFile = ''
    contour1Display.PolarAxes.LastRadialAxisTextFontFile = ''
    contour1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # set scalar coloring
    ColorBy(contour1Display, ('POINTS', 'Vorticity z'))
    
    # rescale color and/or opacity maps used to include current data range
    contour1Display.RescaleTransferFunctionToDataRange(True, False)
    
    # show color bar/color legend
    contour1Display.SetScalarBarVisibility(renderView1, True)
    
    # get color transfer function/color map for 'Vorticityx'
    vorticityxLUT = GetColorTransferFunction('Vorticityz')
    vorticityxLUT.RGBPoints = [.807555, 0.231373, 0.298039, 0.752941, 4.440892098500626e-16, 0.865003, 0.865003, 0.865003, 3.01771, 0.705882, 0.0156863, 0.14902]
    vorticityxLUT.ScalarRangeInitialized = 1.0
    
    #### saving camera placements for all active views
    
    # current camera placement for renderView1
    renderView1.CameraPosition = [73.65260386450542, -266.17953498135154, 58.7226440337288]
    renderView1.CameraFocalPoint = [47.00000000000004, 46.999999999999964, 47.000000000000014]
    renderView1.CameraViewUp = [-0.05299728857042391, 0.03284799252690742, 0.9980542554346109]
    renderView1.CameraParallelScale = 81.40638795573723
    
    #### uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).
    #paraview.simple.SaveScreenshot(jpg)
    ####################################------------------------------------------------------------------------
    addon=frame_names[frame]
    #Render()
    #Interact()
    SaveScreenshot('C:/Users/timvd/Documents/GitHub/AB11/photos/'+ addon + '.jpg')#, ImageResolution=[1172, 826])    #Saves the picture
    ####################################------------------------------------------------------------------------