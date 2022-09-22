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
    ap*tp+ aw*Tw + ae*Te + an*Tn + bc*Tc = 0
    [ap, aw, ae, an, bc]
    '''
    ap = -2 * (hk + 2)

    an = 2
    ae = 1
    aw = 1
    bc = 2*hk*Tc

    a = np.array([ap, aw, ae, an, bc])

    return a


# Vetor de temperatura nodais e matrix de coeficientes
# [Tp] = [A][T] = [B]

dof = (nx) * (ny)
u_dof = np.ones(dof, dtype=bool) # Guarda os nós com temperatura por bc_1

Temperatures = np.ones(dof)*T1
B = np.zeros(dof)
A = np.zeros((dof, dof))

# Constrói os coeficientes adequando as condições de contorno
for i, node in enumerate(mesh.borders):
       
    bc_1 = node['N'] == -1 or node['W'] == -1 or node['E'] == -1

    bc_2 = node['S'] == -1

    if bc_1:
        # Fronteira de temperatura prescrita

        for key, value in node.items():
            if value != -1:
                
                B[value] += T1

        Temperatures[i] = T1
        #A[i,i] = 1
        u_dof[i] = False

    elif bc_2:
        # Fronteira convectiva

        a = boundary2(hk, Tc)

        A[i, i] = a[0]
        A[i, node['W']] = a[1]
        A[i, node['E']] = a[2]
        A[i, node['N']] = a[3]

        B[i] += a[4]

    else:
        # Nó interior
        
        A[i, i] = -4
        A[i, node['W']] = 1
        A[i, node['E']] = 1
        A[i, node['S']] = 1
        A[i, node['N']] = 1

  
# Simplifica as equações conforme bc_1
Auu = A[u_dof, :][:, u_dof]
Buu = B[u_dof]
Tuu = Temperatures[u_dof]

# Solver:



it = 1
max_it = 1e3

diff = T1
tol = 1e-3

print('Iniciando solver')


while diff >= tol:

    Tuu0 = np.copy(Tuu)
    print(f'Número da iteração: {it}')
    
    for i, tp in enumerate(Tuu):


        App = np.append(Auu[i,:i], Auu[i, i+1:])
        Tpp = np.append(Tuu[:i], Tuu[i+1:])
        Bpp = np.append(Buu[:i], Buu[i+1:])

        Tuu[i] = -(np.dot(App, Tpp) + Buu[i]) / Auu[i, i]

    diff = np.sqrt(np.dot(Tuu-Tuu0, Tuu-Tuu0))
    it += 1

    if it > max_it:
        print('Excedido limite de iterações')
        break


# Plots
Tplots = np.ones(dof)*T1
j = 0

for i, tp in enumerate(u_dof):
    if tp:
        Tplots[i] = Tuu[j]
        j += 1


import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6,6))

X = np.linspace(0, 0.4, 142)
Y = np.linspace(0, 0.45, 17767)

Xmesh, Ymesh = np.meshgrid(X,Y)

T = np.random.randint(300, 1000, Xmesh.shape)

plt.pcolormesh(X, Y, T)
plt.colorbar()
plt.show()
