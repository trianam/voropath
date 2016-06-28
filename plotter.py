import numpy as np
import scipy as sp
import scipy.interpolate
import scipy.spatial
import pickle
import vtk
import vtk.util.colors
import math
import warnings
warnings.filterwarnings("error")

class Plotter:
    COLOR_BG = vtk.util.colors.light_grey
    COLOR_BG_PLOT = vtk.util.colors.ghost_white
    #vtk.util.colors.ivory
    COLOR_OBSTACLE = vtk.util.colors.banana
    COLOR_SITES = vtk.util.colors.cobalt
    COLOR_PATH = vtk.util.colors.brick
    COLOR_CONTROL_POINTS = vtk.util.colors.tomato
    COLOR_CONTROL_POLIG = vtk.util.colors.mint
    COLOR_GRAPH = vtk.util.colors.sepia
    COLOR_PLOT_CURV = vtk.util.colors.blue
    COLOR_PLOT_TORS = vtk.util.colors.red
    COLOR_LABELS = vtk.util.colors.blue
    COLOR_LENGTH = vtk.util.colors.red

    _DEFAULT_LINE_THICKNESS = 0.0005
    _DEFAULT_POINT_THICKNESS = 0.002
    _DEFAULT_BSPLINE_THICKNESS = 0.001

    class KeyPressInteractorStyle(vtk.vtkInteractorStyleUnicam):
        _screenshotFile = "/tmp/screenshot.png"
        _cameraFile = "/tmp/cameraData.dat"
        _cameraFile2 = "/tmp/cameraData2.dat"
        def __init__(self, parent=None):
            self.AddObserver("KeyPressEvent", self._keyPressEvent)
            self.AddObserver("RightButtonPressEvent", self._mousePressEvent)
            #super(KeyPressInteractorStyle, self).__init__()

        def SetCamera(self, camera):
            self._camera = camera

        def SetRenderer(self, renderer):
            self._renderer = renderer

        def SetRenderWindow(self, renderWindow):
            self._renderWindow = renderWindow

        def _keyPressEvent(self, obj, event):
            if obj.GetInteractor().GetKeySym() == "l":
                print("Scene screenshot in "+self._screenshotFile)
                w2if = vtk.vtkWindowToImageFilter()
                w2if.SetInput(self._renderWindow)
                w2if.Update()
 
                writer = vtk.vtkPNGWriter()
                writer.SetFileName(self._screenshotFile)
                writer.SetInputData(w2if.GetOutput())
                writer.Write()

            elif obj.GetInteractor().GetKeySym() == "c":
                print("Save camera data in "+self._cameraFile)
                record = {}
                record['position'] = self._camera.GetPosition()
                record['focalPoint'] = self._camera.GetFocalPoint()
                record['viewAngle'] = self._camera.GetViewAngle()
                record['viewUp'] = self._camera.GetViewUp()
                record['clippingRange'] = self._camera.GetClippingRange()

                with open(self._cameraFile, 'wb') as f:
                    pickle.dump(record, f)

            elif obj.GetInteractor().GetKeySym() == "v":
                print("Restore camera data from "+self._cameraFile)

                with open(self._cameraFile, 'rb') as f:
                    record = pickle.load(f)

                    self._camera.SetPosition(record['position'])
                    self._camera.SetFocalPoint(record['focalPoint'])
                    self._camera.SetViewAngle(record['viewAngle'])
                    self._camera.SetViewUp(record['viewUp'])
                    self._camera.SetClippingRange(record['clippingRange'])

                    self._renderWindow.Render()

            elif obj.GetInteractor().GetKeySym() == "b":
                print("Restore camera data from "+self._cameraFile2)

                with open(self._cameraFile2, 'rb') as f:
                    record = pickle.load(f)

                    self._camera.SetPosition(record['position'])
                    self._camera.SetFocalPoint(record['focalPoint'])
                    self._camera.SetViewAngle(record['viewAngle'])
                    self._camera.SetViewUp(record['viewUp'])
                    self._camera.SetClippingRange(record['clippingRange'])

                    self._renderWindow.Render()


            self.OnKeyPress()

        def _mousePressEvent(self, obj, event):
            clickPos = obj.GetInteractor().GetEventPosition()
            picker =vtk.vtkPropPicker()
            picker.Pick(clickPos[0], clickPos[1], 0, self._renderer)
            pos = picker.GetPickPosition()
            print(pos)


    class KeyPressContextInteractorStyle(vtk.vtkContextInteractorStyle):
        _screenshotFile = "/tmp/screenshot.png"
        def __init__(self, parent=None):
            self.AddObserver("KeyPressEvent",self._keyPressEvent)

        def SetRenderWindow(self, renderWindow):
            self._renderWindow = renderWindow

        def _keyPressEvent(self, obj, event):
            if obj.GetInteractor().GetKeySym() == "l":
                print("Plot screenshot in "+self._screenshotFile)
                w2if = vtk.vtkWindowToImageFilter()
                w2if.SetInput(self._renderWindow)
                w2if.Update()
 
                writer = vtk.vtkPNGWriter()
                writer.SetFileName(self._screenshotFile)
                writer.SetInputData(w2if.GetOutput())
                writer.Write()


            
    def __init__(self):
        self._rendererScene = vtk.vtkRenderer()
        self._rendererScene.SetBackground(self.COLOR_BG)

        self._renderWindowScene = vtk.vtkRenderWindow()
        self._renderWindowScene.AddRenderer(self._rendererScene)
        self._renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        self._renderWindowInteractor.SetRenderWindow(self._renderWindowScene)
        #self._interactorStyle = vtk.vtkInteractorStyleUnicam()
        self._interactorStyle = self.KeyPressInteractorStyle()
        self._interactorStyle.SetCamera(self._rendererScene.GetActiveCamera())
        self._interactorStyle.SetRenderer(self._rendererScene)
        self._interactorStyle.SetRenderWindow(self._renderWindowScene)

        self._contextViewPlotCurv = vtk.vtkContextView()
        self._contextViewPlotCurv.GetRenderer().SetBackground(self.COLOR_BG_PLOT)

        self._contextInteractorStyleCurv = self.KeyPressContextInteractorStyle()
        self._contextInteractorStyleCurv.SetRenderWindow(self._contextViewPlotCurv.GetRenderWindow())
        
        self._chartXYCurv = vtk.vtkChartXY()
        self._contextViewPlotCurv.GetScene().AddItem(self._chartXYCurv)
        self._chartXYCurv.SetShowLegend(True)
        self._chartXYCurv.GetAxis(vtk.vtkAxis.LEFT).SetTitle("")
        self._chartXYCurv.GetAxis(vtk.vtkAxis.BOTTOM).SetTitle("")

        self._contextViewPlotTors = vtk.vtkContextView()
        self._contextViewPlotTors.GetRenderer().SetBackground(self.COLOR_BG_PLOT)

        self._contextInteractorStyleTors = self.KeyPressContextInteractorStyle()
        self._contextInteractorStyleTors.SetRenderWindow(self._contextViewPlotTors.GetRenderWindow())
        
        self._chartXYTors = vtk.vtkChartXY()
        self._contextViewPlotTors.GetScene().AddItem(self._chartXYTors)
        self._chartXYTors.SetShowLegend(True)
        self._chartXYTors.GetAxis(vtk.vtkAxis.LEFT).SetTitle("")
        self._chartXYTors.GetAxis(vtk.vtkAxis.BOTTOM).SetTitle("")

        self._textActor = vtk.vtkTextActor()
        self._textActor.GetTextProperty().SetColor(self.COLOR_LENGTH)

        self._addedBSpline = False
        
    def draw(self):
        self._renderWindowInteractor.Initialize()
        self._renderWindowInteractor.SetInteractorStyle(self._interactorStyle)

        axes = vtk.vtkAxesActor()
        widget = vtk.vtkOrientationMarkerWidget()
        widget.SetOutlineColor(0.9300, 0.5700, 0.1300)
        widget.SetOrientationMarker(axes)
        widget.SetInteractor(self._renderWindowInteractor)
        #widget.SetViewport(0.0, 0.0, 0.1, 0.1)
        widget.SetViewport(0.0, 0.0, 0.2, 0.4)
        widget.SetEnabled(True)
        widget.InteractiveOn()

        textWidget = vtk.vtkTextWidget()
 
        textRepresentation = vtk.vtkTextRepresentation()
        textRepresentation.GetPositionCoordinate().SetValue(.0,.0 )
        textRepresentation.GetPosition2Coordinate().SetValue(.3,.04 )
        textWidget.SetRepresentation(textRepresentation)
 
        textWidget.SetInteractor(self._renderWindowInteractor)
        textWidget.SetTextActor(self._textActor)
        textWidget.SelectableOff()
        textWidget.On()

        self._rendererScene.ResetCamera()
        camPos = self._rendererScene.GetActiveCamera().GetPosition()
        self._rendererScene.GetActiveCamera().SetPosition((camPos[2],camPos[1],camPos[0]))
        self._rendererScene.GetActiveCamera().SetViewUp((0.0,0.0,1.0))
        self._rendererScene.GetActiveCamera().Zoom(1.4)

        self._renderWindowScene.Render()

        if self._addedBSpline:
            self._contextViewPlotCurv.GetRenderWindow().SetMultiSamples(0)
            self._contextViewPlotCurv.GetInteractor().Initialize()
            self._contextViewPlotCurv.GetInteractor().SetInteractorStyle(self._contextInteractorStyleCurv)
            #self._contextViewPlotCurv.GetInteractor().Start()

            self._contextViewPlotTors.GetRenderWindow().SetMultiSamples(0)
            self._contextViewPlotTors.GetInteractor().Initialize()
            self._contextViewPlotTors.GetInteractor().SetInteractorStyle(self._contextInteractorStyleTors)
            self._contextViewPlotTors.GetInteractor().Start()
        else:
            self._renderWindowInteractor.Start()
            

    def addTetrahedron(self, vertexes, color):
        vtkPoints = vtk.vtkPoints()
        vtkPoints.InsertNextPoint(vertexes[0][0], vertexes[0][1], vertexes[0][2])
        vtkPoints.InsertNextPoint(vertexes[1][0], vertexes[1][1], vertexes[1][2])
        vtkPoints.InsertNextPoint(vertexes[2][0], vertexes[2][1], vertexes[2][2])
        vtkPoints.InsertNextPoint(vertexes[3][0], vertexes[3][1], vertexes[3][2])

        unstructuredGrid = vtk.vtkUnstructuredGrid()
        unstructuredGrid.SetPoints(vtkPoints)
        unstructuredGrid.InsertNextCell(vtk.VTK_TETRA, 4, range(4))

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(unstructuredGrid)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        self._rendererScene.AddActor(actor)

    def addTriangles(self, triangles, color):
        vtkPoints = vtk.vtkPoints()
        idPoint = 0
        allIdsTriangle = []

        for triangle in triangles:
            idsTriangle = []

            for point in triangle:
                vtkPoints.InsertNextPoint(point[0], point[1], point[2])
                idsTriangle.append(idPoint)
                idPoint += 1

            allIdsTriangle.append(idsTriangle)

        unstructuredGrid = vtk.vtkUnstructuredGrid()
        unstructuredGrid.SetPoints(vtkPoints)
        for idsTriangle in allIdsTriangle:
            unstructuredGrid.InsertNextCell(vtk.VTK_TRIANGLE, 3, idsTriangle)

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(unstructuredGrid)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        self._rendererScene.AddActor(actor)

    def addPolyLine(self, points, color, thick=False, thickness=_DEFAULT_LINE_THICKNESS):
        vtkPoints = vtk.vtkPoints()
        for point in points:
            vtkPoints.InsertNextPoint(point[0], point[1], point[2])

        if thick:
            cellArray = vtk.vtkCellArray()
            cellArray.InsertNextCell(len(points))
            for i in range(len(points)):
                cellArray.InsertCellPoint(i)

            polyData = vtk.vtkPolyData()
            polyData.SetPoints(vtkPoints)
            polyData.SetLines(cellArray)

            tubeFilter = vtk.vtkTubeFilter()
            tubeFilter.SetNumberOfSides(8)
            tubeFilter.SetInputData(polyData)
            tubeFilter.SetRadius(thickness)
            tubeFilter.Update()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(tubeFilter.GetOutputPort())

        else:
            unstructuredGrid = vtk.vtkUnstructuredGrid()
            unstructuredGrid.SetPoints(vtkPoints)
            for i in range(1, len(points)):
                unstructuredGrid.InsertNextCell(vtk.VTK_LINE, 2, [i-1, i])

            mapper = vtk.vtkDataSetMapper()
            mapper.SetInputData(unstructuredGrid)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        self._rendererScene.AddActor(actor)

    def addPoints(self, points, color, thick=False, thickness=_DEFAULT_POINT_THICKNESS):
        vtkPoints = vtk.vtkPoints()
        for point in points:
            vtkPoints.InsertNextPoint(point[0], point[1], point[2])

        pointsPolyData = vtk.vtkPolyData()
        pointsPolyData.SetPoints(vtkPoints)

        if thick:
            sphereSource = vtk.vtkSphereSource()
            sphereSource.SetRadius(thickness)
            
            glyph3D = vtk.vtkGlyph3D()
            glyph3D.SetSourceConnection(sphereSource.GetOutputPort())
            glyph3D.SetInputData(pointsPolyData)
            glyph3D.Update()
 
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(glyph3D.GetOutputPort())
        else:
            vertexFilter = vtk.vtkVertexGlyphFilter()
            vertexFilter.SetInputData(pointsPolyData)
            vertexFilter.Update()                          

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(vertexFilter.GetOutput())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        self._rendererScene.AddActor(actor)

    def addBSpline(self, path, degree, color, thick=False, thickness=_DEFAULT_BSPLINE_THICKNESS):
        self._addedBSpline = True

        tau, u, spline, splineD1, splineD2, splineD3, curv, tors, arcLength, polLength = path.splinePoints()

        self._textActor.SetInput("Length: "+str(arcLength))

        numIntervals = len(tau)-1

        curvPlotActor = vtk.vtkXYPlotActor()
        curvPlotActor.SetTitle("Curvature")
        curvPlotActor.SetXTitle("")
        curvPlotActor.SetYTitle("")
        curvPlotActor.SetXValuesToIndex()
        
        torsPlotActor = vtk.vtkXYPlotActor()
        torsPlotActor.SetTitle("Torsion")
        torsPlotActor.SetXTitle("")
        torsPlotActor.SetYTitle("")
        torsPlotActor.SetXValuesToIndex()

        uArrays = []
        curvArrays = []
        torsArrays = []
        for i in range(numIntervals):
            uArrays.append(vtk.vtkDoubleArray())
            uArrays[i].SetName("t")
            
            curvArrays.append(vtk.vtkDoubleArray())
            curvArrays[i].SetName("Curvature")

            torsArrays.append(vtk.vtkDoubleArray())
            torsArrays[i].SetName("Torsion")
        
        curvTorsArray = vtk.vtkDoubleArray()

        #minCurv = minTors = minNd1Xd2 = float("inf")
        #maxCurv = maxTors = float("-inf")
        
        for i in range(len(u)):
            for j in range(numIntervals):
                if u[i] >= tau[j] and u[i] < tau[j+1]:
                    break
            
            uArrays[j].InsertNextValue(u[i])
            curvArrays[j].InsertNextValue(curv[i])
            torsArrays[j].InsertNextValue(tors[i])
            
            curvTorsArray.InsertNextValue(curv[i])# + abs(tors[i]))

        #print("minCurv: {:e}; maxCurv: {:e}; minTors: {:e}; maxTors: {:e}; minNd1Xd2: {:e}".format(minCurv, maxCurv, minTors, maxTors, minNd1Xd2))

        for inter in range(numIntervals):
            plotTable = vtk.vtkTable()
            plotTable.AddColumn(uArrays[inter])
            plotTable.AddColumn(curvArrays[inter])
            plotTable.AddColumn(torsArrays[inter])
        
            points = self._chartXYCurv.AddPlot(vtk.vtkChart.LINE)
            points.SetInputData(plotTable, 0, 1)
            points.SetColor(self.COLOR_PLOT_CURV[0], self.COLOR_PLOT_CURV[1], self.COLOR_PLOT_CURV[2])
            points.SetWidth(1.0)
            if inter > 0:
                points.SetLegendVisibility(False)
            
            points = self._chartXYTors.AddPlot(vtk.vtkChart.LINE)
            points.SetInputData(plotTable, 0, 2)
            points.SetColor(self.COLOR_PLOT_TORS[0], self.COLOR_PLOT_TORS[1], self.COLOR_PLOT_TORS[2])
            points.SetWidth(1.0)
            if inter > 0:
                points.SetLegendVisibility(False)

            
        vtkPoints = vtk.vtkPoints()
        for point in spline:
            vtkPoints.InsertNextPoint(point[0], point[1], point[2])

        polyDataLabelP = vtk.vtkPolyData()
        polyDataLabelP.SetPoints(vtkPoints)
        
        labels = vtk.vtkStringArray()
        labels.SetNumberOfValues(len(spline))
        labels.SetName("labels")
        for i in range(len(spline)):
            if i == 0:
                labels.SetValue(i, "S")
            elif i == len(spline)-1:
                labels.SetValue(i, "E")
            else:
                labels.SetValue(i, "")

        polyDataLabelP.GetPointData().AddArray(labels)

        sizes = vtk.vtkIntArray()
        sizes.SetNumberOfValues(len(spline))
        sizes.SetName("sizes")
        for i in range(len(spline)):
            if i == 0 or i == len(spline)-1:
                sizes.SetValue(i, 10)
            else:
                sizes.SetValue(i,1)

        polyDataLabelP.GetPointData().AddArray(sizes)
        
        pointMapper = vtk.vtkPolyDataMapper()
        pointMapper.SetInputData(polyDataLabelP)
 
        pointActor = vtk.vtkActor()
        pointActor.SetMapper(pointMapper)

        pointSetToLabelHierarchyFilter = vtk.vtkPointSetToLabelHierarchy()
        pointSetToLabelHierarchyFilter.SetInputData(polyDataLabelP)
        pointSetToLabelHierarchyFilter.SetLabelArrayName("labels")
        pointSetToLabelHierarchyFilter.SetPriorityArrayName("sizes")
        pointSetToLabelHierarchyFilter.GetTextProperty().SetColor(self.COLOR_LABELS)
        pointSetToLabelHierarchyFilter.GetTextProperty().SetFontSize(15)
        pointSetToLabelHierarchyFilter.GetTextProperty().SetBold(True)
        pointSetToLabelHierarchyFilter.Update()
        
        labelMapper = vtk.vtkLabelPlacementMapper()
        labelMapper.SetInputConnection(pointSetToLabelHierarchyFilter.GetOutputPort())
        labelActor = vtk.vtkActor2D()
        labelActor.SetMapper(labelMapper)

        self._rendererScene.AddActor(labelActor)

        if thick:
            cellArray = vtk.vtkCellArray()
            cellArray.InsertNextCell(len(spline))
            for i in range(len(spline)):
                cellArray.InsertCellPoint(i)

            polyData = vtk.vtkPolyData()
            polyData.SetPoints(vtkPoints)
            polyData.SetLines(cellArray)

            polyData.GetPointData().SetScalars(curvTorsArray)
            
            tubeFilter = vtk.vtkTubeFilter()
            tubeFilter.SetNumberOfSides(8)
            tubeFilter.SetInputData(polyData)
            tubeFilter.SetRadius(thickness)
            tubeFilter.Update()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(tubeFilter.GetOutputPort())

        else:
            unstructuredGrid = vtk.vtkUnstructuredGrid()
            unstructuredGrid.SetPoints(vtkPoints)
            for i in range(1, len(spline)):
                unstructuredGrid.InsertNextCell(vtk.VTK_LINE, 2, [i-1, i])

            unstructuredGrid.GetPointData().SetScalars(curvArray)
            
            mapper = vtk.vtkDataSetMapper()
            mapper.SetInputData(unstructuredGrid)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        self._rendererScene.AddActor(actor)


        #self.addPolyLine(list(zip(out[0], out[1], out[2])), color, thick, thickness)

    def addBSplineDEPRECATED(self, controlPolygon, degree, color, thick=False, thickness=_DEFAULT_BSPLINE_THICKNESS):
        x = controlPolygon[:,0]
        y = controlPolygon[:,1]
        z = controlPolygon[:,2]
    
        polLen = 0.
        for i in range(1, len(controlPolygon)):
            polLen += sp.spatial.distance.euclidean(controlPolygon[i-1], controlPolygon[i])

        t = range(len(controlPolygon))
        ipl_t = np.linspace(0.0, len(controlPolygon) - 1, max(polLen*100,100))

        x_tup = sp.interpolate.splrep(t, x, k = degree)
        y_tup = sp.interpolate.splrep(t, y, k = degree)
        z_tup = sp.interpolate.splrep(t, z, k = degree)

        x_list = list(x_tup)
        xl = x.tolist()
        x_list[1] = xl + [0.0, 0.0, 0.0, 0.0]

        y_list = list(y_tup)
        yl = y.tolist()
        y_list[1] = yl + [0.0, 0.0, 0.0, 0.0]

        z_list = list(z_tup)
        zl = z.tolist()
        z_list[1] = zl + [0.0, 0.0, 0.0, 0.0]

        x_i = sp.interpolate.splev(ipl_t, x_list)
        y_i = sp.interpolate.splev(ipl_t, y_list)
        z_i = sp.interpolate.splev(ipl_t, z_list)

        self.addPolyLine(list(zip(x_i, y_i, z_i)), color, thick, thickness)

    def addGraph(self, graph, color):
        vtkPoints = vtk.vtkPoints()
        vtkId = 0
        graph2VtkId = {}
        
        for node in graph.nodes():
            vtkPoints.InsertNextPoint(graph.node[node]['coord'][0], graph.node[node]['coord'][1], graph.node[node]['coord'][2])
            graph2VtkId[node] = vtkId
            vtkId += 1
            
        unstructuredGrid = vtk.vtkUnstructuredGrid()
        unstructuredGrid.SetPoints(vtkPoints)

        for edge in graph.edges():
            unstructuredGrid.InsertNextCell(vtk.VTK_LINE, 2, [graph2VtkId[edge[0]], graph2VtkId[edge[1]]])

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(unstructuredGrid)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        self._rendererScene.AddActor(actor)

    
