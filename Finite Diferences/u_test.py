from meshFd import Mesh


nx = 5
ny = 5
lx = 1
ly = 1

mesh = Mesh(nx-1, ny-1, lx, ly)

#print(mesh.nodes)

print(mesh.borders)