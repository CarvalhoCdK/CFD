import numpy as np

from scipy.sparse import csr_matrix

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

def boundary2(hk):
    '''
    Calcula os termos da equação na forma:
    Tp = aw*Tw + ae*Te + an*Tn + bc*Tc
    [aw, ae, an, bc]
    '''
    an = 1 / (hk + 2)
    ae = an / 2
    aw = ae
    bc = hk / (hk + 2)

    A = np.array([aw, ae, an, bc])

    return A


# Vetor de temperatura nodais e matrix de coeficientes
# [A][T] = [B]

dof = (nx) * (ny)
set_dof = [] # Guarda os nós com temperatura por bc_1

Temperatures = np.zeros(dof)
B = np.zeros(dof)
A = np.zeros((dof, dof))

for i, node in enumerate(mesh.borders):
       
    bc_1 = node['N'] == -1 or node['W'] == -1 or node['E'] == -1

    bc_2 = node['S'] == -1

    if bc_1:
        # print(f'Node: {i} : {bc_1}')
        Temperatures[i] = T1
        set_dof.append(i)

    elif bc_2:
        # W, E, N, B

        a = boundary2(hk)

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






# Set temperature boundary

# Convective boundary
