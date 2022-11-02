import numpy as np
from time import perf_counter


from meshFd import Mesh
from plots import plot_field


# Geometria
nx = 3
ny = 3
lx = 1
ly = 1

# Propriedades
k = 1

# Cronstução da malha
mesh = Mesh(nx-1, ny-1, lx, ly)

nodes = mesh.nodes
borders = mesh.borders

# Propriedades do solver

max_it = 1e3
tol = 0.01#1e-3

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


# Resolvendo Iterativamente
diff = T1
it = 1

print('Iniciando solver')
start_time = perf_counter()

while diff >= tol:

    #print(f'Número da iteração: {it}')
    
    Tuu0 = np.copy(Tuu)
        
    for i, tp in enumerate(Tuu):


        #App = np.append(Auu[i,:i], Auu[i, i+1:])
        #Tpp = np.append(Tuu[:i], Tuu[i+1:])
        #Bpp = np.append(Buu[:i], Buu[i+1:])

        #Tuu[i] = -(np.dot(App, Tpp) + Buu[i]) / Auu[i, i]

        pfront = np.dot(Auu[i,:i], Tuu[:i])
        pback = np.dot(Auu[i, i+1:], Tuu[i+1:])

        Tuu[i] = -(pfront + pback + Buu[i]) / Auu[i, i]

    diff = np.max(np.abs(Tuu-Tuu0))# np.sqrt(np.dot(Tuu-Tuu0, Tuu-Tuu0))
    it += 1

    if it > max_it:
        print('Excedido limite de iterações')
        break

end_time = perf_counter()
print(f'Tempo de execução: {end_time - start_time}')
print(f'Número de iterações: {it}')
print(f'Erro: {diff}')
print(f'Temperatura no centro da face sul: {Tuu[int(np.floor(nx/2))]}')

# Plots
Tplot = np.ones(dof)*T1
j = 0

for i, tp in enumerate(u_dof):
    if tp:
        Tplot[i] = Tuu[j]
        j += 1

shading = 'nearest'
plot_field(nx, ny, lx, ly, Tplot)

print('EoF')

