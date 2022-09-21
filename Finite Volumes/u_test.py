from mesh import Mesh


nx = 3
ny = 3
lx = 3
ly = 3

mesh = Mesh(nx, ny, lx, ly)

print(mesh.elements)

print(mesh.borders)