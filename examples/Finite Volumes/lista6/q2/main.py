import numpy as np
import matplotlib.pyplot as plt


from mesh import Mesh
from interpolations import WUDS
from solver import gauss_seidel

## PARAMETROS

## MALHA
nx = 50
ny = 50
a = 1
b = 1
tol = 1e-4


model = {
    'u' : 0.0,  # Velocidade em x
    'v' : 0.0,  # Velocidade em y
    'rho' : 1.0,  # Densidade
    'gamma' : 1.0,
    'a' : a,
    'b' : b,
    'deltax' : a/nx,  # Número de elementos em x
    'deltay' : b/ny,   # Número de elementos em y
    'sp' : 0.0
    }

## Condições de contorno
tw = 0
te = 0
tn = 0

def ts(x,a):
    t = np.sin(np.pi*x/a)
    return t

## Expande malha para volumes fictícios
nxe = nx + 2
nye = ny + 2
# Nova origem dos eixos x,y
xe0 = model['deltax']
ye0 = model['deltay']

a = a + 2*xe0
b = b + 2*ye0
# Cronstução da malha expandida
mesh = Mesh(nxe, nye, a, b)

volumes = mesh.volumes
neighbours = mesh.neighbours


interp = WUDS(model)

# Anota posições
dof = (nxe) * (nye)
u_dof = np.ones(dof, dtype=bool)
real_t = np.zeros(dof, dtype=bool)

B = np.zeros(dof)
A = np.zeros((dof, dof))

# TDMA
am = np.zeros(dof)
bm = np.zeros(dof)
cm = np.zeros(dof)
dm = np.zeros(dof)

for i, nb in enumerate(neighbours):
    # [Ap, Aw, Ae, As, An, B]
    # [0,  1,  2,  3,  4,  5]
    # w, e, s, n
    nb = np.array(nb, dtype="int")

    if np.sum(nb < 0) > 1:
        print(f'{i} : Corner')

        u_dof[i] = False
        # Não há volumes fictícios nos cantos
    
    # Fronteira W
    elif nb[0] == -1:
        #print(f'{i} : W')
        #print(interp.boundary('W', 0))
        aa = interp.boundary('W', tw)
        A[i, i] = aa[0]
        A[i, nb[1]] = -aa[2]
        B[i] = aa[-1]

    # Fronteira E
    elif nb[1] == -1:
        #print(f'{i} : E')
        #print(interp.boundary('E', 0))
        aa = interp.boundary('E', te)
        A[i, i] = aa[0]
        A[i, nb[0]] = -aa[1]
        B[i] = aa[-1]
    
    # Fronteira N
    elif nb[3] == -1:
        #print(f'{i} : N')
        #print(interp.boundary('N', 0))
        aa = interp.boundary('N', tn)
        A[i, i] = aa[0]
        A[i, nb[2]] = -aa[3]
        B[i] = aa[-1]

    # Fronteira S
    elif nb[2] == -1:
        #print(f'{i} : S')
        #print(interp.boundary('S', ts(x, a)))
        x = volumes[i][0]
        aa = interp.boundary('S', ts(x, a))
        A[i, i] = aa[0]
        A[i, nb[3]] = -aa[4]
        B[i] = aa[-1]

    # Interno
    else:
        #print(f'{i} : Interno')
        #print(interp.internal())
        real_t[i] = True
        a = interp.internal()
        A[i, i] = a[0]
        A[i, nb[0]] = -a[1]
        A[i, nb[1]] = -a[2]
        A[i, nb[2]] = -a[3]
        A[i, nb[3]] = -a[4]

# Exclui cantos e separa os volumes não fictícios            
Auu = A[u_dof, :][:, u_dof]
Buu = B[u_dof]

real_tuu = real_t[u_dof]

# Solução do sistema
#T = np.linalg.solve(Auu, Buu)

print('Starting solver')
T, itj = gauss_seidel(Auu, Buu, tol)

## PLOTS
x = volumes[real_t][:,0]
y = volumes[real_t][:,1]
z = T[real_tuu]

z_min, z_max = np.abs(z).min(), np.abs(z).max()

x = x.reshape((nx, ny))-xe0
y = y.reshape((nx, ny))-ye0
z = z.reshape((nx, ny))

fig, ax = plt.subplots()#figsize=(6, 6))
c = ax.pcolormesh(x, y, z, cmap='coolwarm',  vmin=z_min, vmax=z_max)
ax.axis([x.min(), x.max(), y.min(), y.max()])
fig.colorbar(c, ax=ax, label='T')
ax.set_xlabel('x', fontsize=14)  
ax.set_ylabel('y', fontsize=14)
ax.grid()


## Solução analítica : u = v = 0
a = model['a']
b = model['b']

t_anl = np.sinh(np.pi*y/a) * np.sin(np.pi*x/b) / np.sinh(np.pi*b/a)
z1 = t_anl.reshape((nx, ny))

fig, ax = plt.subplots()#figsize=(6, 6))
c = ax.pcolormesh(x, y, np.flip(z1), cmap='coolwarm',  vmin=z_min, vmax=z_max)
ax.axis([x.min(), x.max(), y.min(), y.max()])
fig.colorbar(c, ax=ax, label='T')
ax.set_xlabel('x', fontsize=14)  
ax.set_ylabel('y', fontsize=14)
ax.grid()
