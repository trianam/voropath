import numpy as np
import scipy as sp
import scipy.interpolate
import vtk
import vtk.util.colors

class Plotter:
    COLOR_BG = vtk.util.colors.ivory
    COLOR_OBSTACLE = vtk.util.colors.banana
    COLOR_SITES = vtk.util.colors.cobalt
    COLOR_PATH = vtk.util.colors.brick
    COLOR_CONTROL_POLIG = vtk.util.colors.mint
    COLOR_GRAPH = vtk.util.colors.sepia

    def __init__(self):
        self._renderer = vtk.vtkRenderer()
        self._renderer.SetBackground(self.COLOR_BG)

        self._renderWindow = vtk.vtkRenderWindow()
        self._renderWindow.AddRenderer(self._renderer)
        self._renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        self._renderWindowInteractor.SetRenderWindow(self._renderWindow)
        self._interactorStyle = vtk.vtkInteractorStyleUnicam()

    def draw(self):
        self._renderWindowInteractor.Initialize()
        self._renderWindowInteractor.SetInteractorStyle(self._interactorStyle)

        axes = vtk.vtkAxesActor()
        widget = vtk.vtkOrientationMarkerWidget()
        widget.SetOutlineColor(0.9300, 0.5700, 0.1300)
        widget.SetOrientationMarker(axes)
        widget.SetInteractor(self._renderWindowInteractor)
        widget.SetViewport(0.0, 0.0, 0.1, 0.1)
        widget.SetEnabled(True)
        widget.InteractiveOn()

        self._renderer.ResetCamera()
        camPos = self._renderer.GetActiveCamera().GetPosition()
        self._renderer.GetActiveCamera().SetPosition((camPos[2],camPos[1],camPos[0]))
        self._renderer.GetActiveCamera().SetViewUp((0.0,0.0,1.0))

        self._renderWindow.Render()
        self._renderWindowInteractor.Start()

    def addTetrahedron(self, vertexes, color):
        vtkPoints = vtk.vtkPoints()
        vtkPoints.InsertNextPoint(vertexes[0][0], vertexes[0][1], vertexes[0][2])
        vtkPoints.InsertNextPoint(vertexes[1][0], vertexes[1][1], vertexes[1][2])
        vtkPoints.InsertNextPoint(vertexes[2][0], vertexes[2][1], vertexes[2][2])
        vtkPoints.InsertNextPoint(vertexes[3][0], vertexes[3][1], vertexes[3][2])

        unstructuredGrid = vtk.vtkUnstructuredGrid()
        unstructuredGrid.SetPoints(vtkPoints)
        unstructuredGrid.InsertNextCell(vtk.VTK_TETRA, 4, [0,1,2,3])

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(unstructuredGrid)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        self._renderer.AddActor(actor)

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

        self._renderer.AddActor(actor)

    def addPolyLine(self, points, color):
        vtkPoints = vtk.vtkPoints()
        for point in points:
            vtkPoints.InsertNextPoint(point[0], point[1], point[2])

        unstructuredGrid = vtk.vtkUnstructuredGrid()
        unstructuredGrid.SetPoints(vtkPoints)
        for i in range(1, len(points)):
            unstructuredGrid.InsertNextCell(vtk.VTK_LINE, 2, [i-1, i])

        mapper = vtk.vtkDataSetMapper()
        mapper.SetInputData(unstructuredGrid)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        self._renderer.AddActor(actor)

    def addPoints(self, points, color):
        vtkPoints = vtk.vtkPoints()
        for point in points:
            vtkPoints.InsertNextPoint(point[0], point[1], point[2])

        pointsPolyData = vtk.vtkPolyData()
        pointsPolyData.SetPoints(vtkPoints)

        vertexFilter = vtk.vtkVertexGlyphFilter()
        vertexFilter.SetInputData(pointsPolyData)
        vertexFilter.Update()                          

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(vertexFilter.GetOutput())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color)

        self._renderer.AddActor(actor)

    def addBSpline(self, controlPolygon, degree, color):
        x = controlPolygon[:,0]
        y = controlPolygon[:,1]
        z = controlPolygon[:,2]

        t = range(len(controlPolygon))
        ipl_t = np.linspace(0.0, len(controlPolygon) - 1, 1000)
        #TODO: find a better way to substitute 1000 above
        x_tup = sp.interpolate.splrep(t, x, k = degree + 1)
        y_tup = sp.interpolate.splrep(t, y, k = degree + 1)
        z_tup = sp.interpolate.splrep(t, z, k = degree + 1)

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

        self.addPoints(zip(x_i, y_i, z_i), color)

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

        self._renderer.AddActor(actor)
