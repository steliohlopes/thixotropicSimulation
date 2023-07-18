import meshio
from dolfin import *

class Mesh:
    def __init__(self,meshPath,meshFile):
        self.meshPath=meshPath
        self.meshFile=meshFile
        self.msh = meshio.read(meshPath+meshFile+".msh")

    def readDomains(self):
    # Read .msh File
        fid = open(self.meshPath+self.meshFile+'.msh', 'r')
        # Initialize variables
        found = 0
        finished = 0
        physicalNames = {}
        # Loop througn .msh lines
        for line in fid:
            if '$EndPhysicalNames' in line:
                finished == 1
                break
            elif '$PhysicalNames' in line:
                found = 1
            elif found==1 and finished == 0:
                word=line.split()
                if len(word)==3:
                    physicalNames[word[2][1:len(word[2])-1]] = int(word[1])

        return physicalNames
    
    def msh2hdmf3D(self):
        for key in self.msh.cell_data_dict["gmsh:physical"].keys():
            if key == "triangle":
                triangle_data = self.msh.cell_data_dict["gmsh:physical"][key]
            elif key == "tetra":
                tetra_data = self.msh.cell_data_dict["gmsh:physical"][key]
        for cell in self.msh.cells:
            if cell.type == "tetra":
                tetra_cells = cell.data
            elif cell.type == "triangle":
                triangle_cells = cell.data
        tetra_mesh = meshio.Mesh(points=self.msh.points, cells={"tetra": tetra_cells},
                                cell_data={"name_to_read":[tetra_data]})
        triangle_mesh =meshio.Mesh(points=self.msh.points,
                                cells=[("triangle", triangle_cells)],
                                cell_data={"name_to_read":[triangle_data]})

        meshio.write(self.meshPath+"mesh.xdmf", tetra_mesh)
        meshio.write(self.meshPath+"mf.xdmf", triangle_mesh)

    def msh2hdmf2D(self):
        for key in self.msh.cell_data_dict["gmsh:physical"].keys():
            if key == "line":
                line_data = self.msh.cell_data_dict["gmsh:physical"][key]
            elif key == "triangle":
                triangle_data = self.msh.cell_data_dict["gmsh:physical"][key]

        for cell in self.msh.cells:
            if cell.type == "triangle":
                triangle_cells = cell.data
            elif cell.type == "line":
                line_cells = cell.data

        triangle_mesh = meshio.Mesh(
            points=self.msh.points[:,:2],
            cells={"triangle": triangle_cells},
            cell_data={"name_to_read": [triangle_data]})

        line_mesh = meshio.Mesh(
            points=self.msh.points[:,:2],
            cells=[("line", line_cells)],
            cell_data={"name_to_read": [line_data]})

        meshio.write(self.meshPath+"mesh.xdmf", triangle_mesh)
        meshio.write(self.meshPath+"mf.xdmf", line_mesh)

    def meshObject3D(self):
        self.msh2hdmf3D()
        self.meshObj = Mesh()

        with XDMFFile(self.meshPath+"mesh.xdmf") as infile:
            infile.read(self.meshObj)
        mvc = MeshValueCollection("size_t", self.meshObj, 2)
        with XDMFFile(self.meshPath+"mf.xdmf") as infile:
            infile.read(mvc, "name_to_read")
        mf = cpp.mesh.MeshFunctionSizet(self.meshObj, mvc)

        mvc2 = MeshValueCollection("size_t", self.meshObj, 3)
        with XDMFFile(self.meshPath+"mesh.xdmf") as infile:
            infile.read(mvc2, "name_to_read")
        cf = cpp.mesh.MeshFunctionSizet(self.meshObj, mvc2)

        # Get Element Shape: Triangle, etc...
        elementShape = self.meshObj.ufl_cell()  
        
        # Define any measure associated with domain and subdomains
        dx = Measure('dx', domain=meshObj, subdomain_data=cf)
        ds = Measure('ds', domain=meshObj, subdomain_data=mf)



    
