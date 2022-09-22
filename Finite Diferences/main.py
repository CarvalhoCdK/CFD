import numpy as np

from scipy.sparse import csr_matrix
from scipy import linalg

from meshFd import Mesh


# Geometria
nx = 5
ny = 5
lx = 1
ly = 1

# Propriedades
k = 1

# Cronstução da malha
mesh = Mesh(nx-1, ny-1, lx, ly)

nodes = mesh.nodes
borders = mesh.borders

# Condições de fronteira
# Tipo 1: temperatura mantida em T1 = 500k
# bc_1

T1 = 500

# Tipo 2: fronteira convectiva com h = 10 W/(m2 K) e Tc = 300 K
# bc_2
h = 10
Tc = 300

hk = h * mesh.deltax / k

def boundary2(hk, Tc):
    '''
    Calcula os termos da equação na forma:
    Tp = aw*Tw + ae*Te + an*Tn + bc*Tc
    [aw, ae, an, bc]
    '''
    an = 1 / (hk + 2)
    ae = an / 2
    aw = ae
    bc = hk / (hk + 2)

    a = np.array([aw, ae, an, bc*Tc])

    return a


# Vetor de temperatura nodais e matrix de coeficientes
# [Tp] = [A][T] = [B]

dof = (nx) * (ny)
u_dof = np.ones(dof, dtype=bool) # Guarda os nós com temperatura por bc_1

Temperatures = np.zeros(dof)
B = np.zeros(dof)
A = np.zeros((dof, dof))

for i, node in enumerate(mesh.borders):
       
    bc_1 = node['N'] == -1 or node['W'] == -1 or node['E'] == -1

    bc_2 = node['S'] == -1

    if bc_1:
        # print(f'Node: {i} : {bc_1}')
        Temperatures[i] = T1
        A[i,i] = 1
        u_dof[i] = False

    elif bc_2:
        # W, E, N, B

        a = boundary2(hk, Tc)

        A[i, node['W']] = a[0]
        A[i, node['E']] = a[1]
        A[i, node['N']] = a[2]

        B[i] = a[3]

    else:
        ai = 1 / 4
        A[i, node['W']] = ai
        A[i, node['E']] = ai
        A[i, node['S']] = ai
        A[i, node['N']] = ai

  
# Simplifica as equações conforme bc_1
Auu = A#[u_dof, :][:, u_dof]
Buu = B#[u_dof]
Tuu = Temperatures#[u_dof]

# Solver:



it = 1
max_it = 1e3

diff = T1
tol = 1e-2

print('Iniciando solver')

while diff >= tol:

    Tuu0 = np.copy(Tuu)
    print(f'Número da iteração: {it}')
    
    for i, tp in enumerate(Tuu):

        Tuu[i] = np.dot(Auu[0,:], Tuu) + Buu[i]

    diff = np.sqrt(np.dot(Tuu-Tuu0, Tuu-Tuu0))
    it += 1

    if it > max_it:
        print('Excedido limite de iterações')
        break


print(Tuu)

# Set temperature boundary

# Convective boundary
