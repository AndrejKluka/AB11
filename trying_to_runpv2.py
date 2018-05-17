# -*- coding: utf-8 -*-
"""
Created on Mon May 14 16:28:12 2018

@author: timvd
"""

import time
from paraview.simple import *
stop11 = time.clock()

frames=20
frame_names=[]
if frames>0:
    for i in range(frames):
        frame_names.append('uvwp_00{:03}.h5-Q.vtr' .format(i+1))

times=1   
if frames!=1:
    times=frames

for frame in range(times) : 
# =============================================================================
#     #import name
#     ##################################
#     colouring = 'Vorticity z'#raw_input("What type of colouring to you want? Normal, Vorticityx, Vorticityy of Vorticityz? ")
#     Isosurface = .2 #float(raw_input("What is the isosurface value?"))
#     ##################################
# =============================================================================
    uvwp_00001h5Qvtr = XMLRectilinearGridReader(FileName=['C:\\Users\\timvd\\Documents\\GitHub\\AB11\\calculated data\\{}'.format(frame_names[frame])])
    uvwp_00001h5Qvtr.PointArrayStatus = ['Q', 'Vorticity normal', 'Vorticity z']
    
    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1172, 826]
    
    # show data in view
    uvwp_00001h5QvtrDisplay = Show(uvwp_00001h5Qvtr, renderView1)
    
    # trace defaults for the display properties.
    uvwp_00001h5QvtrDisplay.Representation = 'Outline'
    uvwp_00001h5QvtrDisplay.ColorArrayName = [None, '']
    uvwp_00001h5QvtrDisplay.OSPRayScaleArray = 'Q'
    uvwp_00001h5QvtrDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    uvwp_00001h5QvtrDisplay.SelectOrientationVectors = 'None'
    uvwp_00001h5QvtrDisplay.ScaleFactor = 19.1
    uvwp_00001h5QvtrDisplay.SelectScaleArray = 'None'
    uvwp_00001h5QvtrDisplay.GlyphType = 'Arrow'
    uvwp_00001h5QvtrDisplay.GlyphTableIndexArray = 'None'
    uvwp_00001h5QvtrDisplay.GaussianRadius = 0.9550000000000001
    uvwp_00001h5QvtrDisplay.SetScaleArray = ['POINTS', 'Q']
    uvwp_00001h5QvtrDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    uvwp_00001h5QvtrDisplay.OpacityArray = ['POINTS', 'Q']
    uvwp_00001h5QvtrDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    uvwp_00001h5QvtrDisplay.DataAxesGrid = 'GridAxesRepresentation'
    uvwp_00001h5QvtrDisplay.SelectionCellLabelFontFile = ''
    uvwp_00001h5QvtrDisplay.SelectionPointLabelFontFile = ''
    uvwp_00001h5QvtrDisplay.PolarAxes = 'PolarAxesRepresentation'
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    uvwp_00001h5QvtrDisplay.ScaleTransferFunction.Points = [-3266.90625, 0.0, 0.5, 0.0, 1098.2564697265625, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    uvwp_00001h5QvtrDisplay.OpacityTransferFunction.Points = [-3266.90625, 0.0, 0.5, 0.0, 1098.2564697265625, 1.0, 0.5, 0.0]
    
    # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
    uvwp_00001h5QvtrDisplay.DataAxesGrid.XTitleFontFile = ''
    uvwp_00001h5QvtrDisplay.DataAxesGrid.YTitleFontFile = ''
    uvwp_00001h5QvtrDisplay.DataAxesGrid.ZTitleFontFile = ''
    uvwp_00001h5QvtrDisplay.DataAxesGrid.XLabelFontFile = ''
    uvwp_00001h5QvtrDisplay.DataAxesGrid.YLabelFontFile = ''
    uvwp_00001h5QvtrDisplay.DataAxesGrid.ZLabelFontFile = ''
    
    # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
    uvwp_00001h5QvtrDisplay.PolarAxes.PolarAxisTitleFontFile = ''
    uvwp_00001h5QvtrDisplay.PolarAxes.PolarAxisLabelFontFile = ''
    uvwp_00001h5QvtrDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
    uvwp_00001h5QvtrDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
    
    # reset view to fit data
    renderView1.ResetCamera()
    
    # update the view to ensure updated data information
    renderView1.Update()
    
    # create a new 'Contour'
    contour1 = Contour(Input=uvwp_00001h5Qvtr)
    contour1.ContourBy = ['POINTS', 'Q']
    contour1.Isosurfaces = [-1084.3248901367188]
    contour1.PointMergeMethod = 'Uniform Binning'
    
    # Properties modified on contour1
    contour1.Isosurfaces = [100.0]
    
    # show data in view
    contour1Display = Show(contour1, renderView1)
    
    # trace defaults for the display properties.
    contour1Display.Representation = 'Surface'
    contour1Display.ColorArrayName = [None, '']
    contour1Display.OSPRayScaleArray = 'Normals'
    contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    contour1Display.SelectOrientationVectors = 'None'
    contour1Display.ScaleFactor = 19.06903533935547
    contour1Display.SelectScaleArray = 'None'
    contour1Display.GlyphType = 'Arrow'
    contour1Display.GlyphTableIndexArray = 'None'
    contour1Display.GaussianRadius = 0.9534517669677735
    contour1Display.SetScaleArray = ['POINTS', 'Normals']
    contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
    contour1Display.OpacityArray = ['POINTS', 'Normals']
    contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
    contour1Display.DataAxesGrid = 'GridAxesRepresentation'
    contour1Display.SelectionCellLabelFontFile = ''
    contour1Display.SelectionPointLabelFontFile = ''
    contour1Display.PolarAxes = 'PolarAxesRepresentation'
    
    # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    contour1Display.ScaleTransferFunction.Points = [-0.9999080300331116, 0.0, 0.5, 0.0, 0.9999626278877258, 1.0, 0.5, 0.0]
    
    # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    contour1Display.OpacityTransferFunction.Points = [-0.9999080300331116, 0.0, 0.5, 0.0, 0.9999626278877258, 1.0, 0.5, 0.0]
    
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
    
    # get color transfer function/color map for 'Vorticityz'
    vorticityzLUT = GetColorTransferFunction('Vorticityz')
    vorticityzLUT.RGBPoints = [-0.9996499419212341, 0.231373, 0.298039, 0.752941, -0.0003439188003540039, 0.865003, 0.865003, 0.865003, 0.9989621043205261, 0.705882, 0.0156863, 0.14902]
    vorticityzLUT.ScalarRangeInitialized = 1.0
    
    #### saving camera placements for all active views
    
    # current camera placement for renderView1
    renderView1.CameraPosition = [95.5, 95.5, 734.5984560556852]
    renderView1.CameraFocalPoint = [95.5, 95.5, 95.5]
    renderView1.CameraParallelScale = 165.41085212282778
    
    #### uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).
    ####################################------------------------------------------------------------------------
    #Interact()
    SaveScreenshot('C:/Users/timvd/Documents/GitHub/AB11/photos/'+ frame_names[frame] + '.jpg')#, ImageResolution=[1172, 826])    #Saves the picture
    ####################################------------------------------------------------------------------------
#for loading the right file validation_Q_l295of96Qvtr = XMLRectilinearGridReader(FileName=["C:\\Users\\timvd\\Documents\\GitHub\\AB11\\calculated data\\{}".format(frame_names[frame])])
print ('\n',int((time.clock()-stop11)*10000)/10000.,'sec  calculations done')
time.sleep(10)