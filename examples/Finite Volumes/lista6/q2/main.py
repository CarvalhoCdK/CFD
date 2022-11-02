import numpy as np
import matplotlib.pyplot as plt


from mesh import Mesh
from interpolations import WUDS
from plots import plot_field


## PARAMETROS

## MALHA
nx = 50
ny = 50
a = 1
b = 1


model = {
    'u' : 10.0,  # Velocidade em x
    'v' : 5.0,  # Velocidade em y
    'rho' : 1.0,  # Densidade
    'gamma' : 1.0,
    'deltax' : 1.0/nx,  # Número de elementos em x
    'deltay' : 1.0/ny   # Número de elementos em y
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

dof = (nxe) * (nye)
u_dof = np.ones(dof, dtype=bool)
real_t = np.zeros(dof, dtype=bool)

B = np.zeros(dof)
A = np.zeros((dof, dof))

for i, nb in enumerate(neighbours):
    # [Ap, Aw, Ae, As, An, B]
    # [0,  1,  2,  3,  4,  5]
    # w, e, s, n
    #print(nb)
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

        x = volumes[i][0]
        #print(interp.boundary('S', ts(x, a)))

        aa = interp.boundary('S', ts(x, a))
        A[i, i] = aa[0]
        A[i, nb[3]] = -aa[4]
        B[i] = aa[-1]

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


Auu = A[u_dof, :][:, u_dof]
Buu = B[u_dof]

real_tuu = real_t[u_dof]


t = np.linalg.solve(Auu, Buu)
#print(t)

x = volumes[real_t][:,0]
y = volumes[real_t][:,1]
z = t[real_tuu]

z_min, z_max = np.abs(z).min(), np.abs(z).max()

x = x.reshape((nx, ny))
y = y.reshape((nx, ny))
z = z.reshape((nx, ny))

fig, ax = plt.subplots(figsize=(6, 6))
ax.pcolormesh(x, y, z, cmap='coolwarm',  vmin=z_min, vmax=z_max)
ax.axis([x.min(), x.max(), y.min(), y.max()])

#plot_field(nx, ny, a, b, t[real_tuu])