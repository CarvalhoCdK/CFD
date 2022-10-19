import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter

from mesh import Mesh
from solver import jacobi
from solver import tdma



def lista5_1(n, tol, model):
    #n = 100

    L = model['L']
    D = model['D']

    # Propriedades


    #tol = 0.001
    max_it = 1e6


    def dirichlet(model):
        """
        Aplica condição de contorno de temperatura prescrita.
        Equação do volume tem forma:
        ap Tp = ae Te + aw Tw + B + a0 T0
        """
        H = model['H']
        K = model['K']
        Alpha = model['Alpha']

        Tb = model['Tb']
        Tamb = model['Tamb']
        
        dx = model['dx']
        dt = model['dt']

        hh = 4*H / (K*D)

        a0 = dx / (Alpha*dt)
        ap = a0 + 3 / dx + hh*dx
        ae = 1 / dx
        b = hh*Tamb*dx + 2*Tb / dx
        
        return np.array([ap, 0, ae, a0, b])


    def newmann(model):
        """
        Aplica a condição de contorno de fluxo prescrito.
        Equação do volume tem forma:
        ap Tp = ae Te + aw Tw + B + a0 T0
        """
        
        H = model['H']
        K = model['K']
        Alpha = model['Alpha']

        Tb = model['Tb']
        Tamb = model['Tamb']

        dx = model['dx']
        dt = model['dt']
        hh = 4*H / (K*D)

        a0 = dx / (Alpha*dt)
        ap = a0 + 1 / dx + hh*dx
        aw = 1 / dx
        b = hh*Tamb*dx
        
        return np.array([ap, aw, 0, a0, b])


    def internal_volume(model):
        """
        Implementa um volume interno genérico.
        Equação do volume tem forma:
        ap Tp = ae Te + aw Tw + B + a0 T0
        """
        
        H = model['H']
        K = model['K']
        Alpha = model['Alpha']
        
        Tamb = model['Tamb']

        dx = model['dx']
        dt = model['dt']
        hh = 4*H / (K*D)

        a0 = dx / (Alpha*dt)
        ap = a0 + 2 / dx + hh*dx
        aw = ae = 1 / dx
        b = hh*Tamb*dx
        
        return np.array([ap, aw, ae, a0, b])

    # Constroi matrizes da forma:
    # [A][T] = [b]+[A0][T0] 

    A = np.zeros((n,n))
    A0 = np.zeros(n)
    b = np.zeros(n)

    for i in range(n):

        if i == 0:
            # Primeiro volume
            a = dirichlet(model)
            A[i, i+1] = -a[2]

            #print(f'{i} : Dirichlet')

        elif i == n-1:
            # Ultimo volume
            a = newmann(model)
            A[i, i-1] = -a[1]

            #print(f'{i} : Newmann')

        else:
            a = internal_volume(model)
            A[i, i-1] = -a[1]
            A[i, i+1] = -a[2]
            #print(f'{i} : Interno')
        
        A[i, i] = a[0]
        
        A0[i] = a[3]
        b[i] = a[4]

    #print(A)

    # SOLVER
    diff = model['Tb']
    it = 1

    # Condição inicial
    T = model['T0'] * np.ones(n)
    #T00 = np.copy(T)
    #B = b +  np.multiply(A0, T00)

    print('Iniciando solver')
    start_time = perf_counter()
    while diff >= tol and it < max_it:

        #print(it)
        #print(T)

        T00 = np.copy(T)

        B = b + np.multiply(A0, T00)

        T = tdma(A, B, T00)
        #T, itj = jacobi(A, B, T00)
        #print(itj)

        diff = np.max(np.abs(T - T00))
        it += 1

    end_time = perf_counter()
    print(f'Solução {n} volumes')
    print(f'Tempo de execução: {end_time - start_time}')
    print(f'Número de iterações: {it}')
    #print(T)

    return T


n = 100
L = 0.05
D = 0.01

model = {'K' : 10.0,
        'H': 5.0,
        'Alpha': 1e-6,
        'Tb' : 373,
        'Tamb' :293,
        'T0' : 293,
        'dx' : L/(n),
        'dt' : 0.5,
        'L' : L,
        'D' : D}

###############################################################################
## PLOT VOLUMES DIFERENTES
###############################################################################
L = model['L']

## 10 Volumes
n = 10
tol = 1e-3

model['dx'] = L/n
x10 = np.linspace(0,L,n)
T10 = lista5_1(n, tol, model)

## 20 Volumes
n = 20
tol = 1e-3

model['dx'] = L/n
x20 = np.linspace(0,L,n)
T20 = lista5_1(n, tol, model)

## 40 Volumes
n = 40
tol = 1e-3

model['dx'] = L/n
x40 = np.linspace(0,L,n)
T40 = lista5_1(n, tol, model)

## 80 Volumes
n = 80
tol = 1e-3

model['dx'] = L/n
x80 = np.linspace(0,L,n)
T80 = lista5_1(n, tol, model)

## 160 Volumes
n = 160
tol = 1e-3

model['dx'] = L/n
x160 = np.linspace(0,L,n)
T160 = lista5_1(n, tol, model)

## 100 Volumes
#n = 100
#tol = 1e-3

#model['dx'] = L/n
#x100 = np.linspace(0,L,n)
#T100 = lista5_1(n, tol, model)

## Analítico
def analitic(n, model):

    H = model['H']
    K = model['K']
    
    Tamb = model['Tamb']
    Tb = model['Tb']

    hh = 4*H / (K*D)
    m = hh**(1/2)
    x = np.linspace(0,L,n)

    c = np.cosh(m * (L-x)) / np.cosh(m*L)

    T = c * (Tb - Tamb) + Tamb

    return x, T
    #x = np.linspace(0,L,n)

xanl, Tanl = analitic(100, model)

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(x10, T10, label='10 Volumes', linestyle='dashed', marker='.')
ax.plot(x20, T20, label='20 Volumes', linestyle='dashed', marker='+')
ax.plot(x40, T40, label='40 Volumes', linestyle='dotted')#, marker='.')
ax.plot(x80, T80, label='80 Volumes', linestyle='dashdot')#, marker='.')  
ax.plot(x160, T160, label='160 Volumes', linestyle='dashed')#, marker='.')
#ax.plot(xanl, Tanl, label='Solução analítica')
ax.set_xlabel('x [m]', fontsize=14)  
ax.set_ylabel('Temperatura [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=12)
#ax.set_title('Soluções para diferentes discretizações.') 
ax.legend()  
ax.grid()
fig.show()

###############################################################################
## PLOT TOLERÂNCIAS DIFERENTES
###############################################################################
L = model['L']
n = 80
x = np.linspace(0,L,n)

##
tol = 1e-3

model['dx'] = L/n
T3 = lista5_1(n, tol, model)

## 
tol = 1e-4

model['dx'] = L/n
T4 = lista5_1(n, tol, model)

## 
tol = 1e-5

model['dx'] = L/n
T5 = lista5_1(n, tol, model)

## 
tol = 1e-6

model['dx'] = L/n
T6 = lista5_1(n, tol, model)

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(x, T3, label='Tolerância : 1e-3', linestyle='dotted')#, marker='.')
ax.plot(x, T4, label='Tolerância : 1e-4', linestyle='dashdot')#, marker='.')
ax.plot(x, T5, label='Tolerância : 1e-5', linestyle='dashed')#, marker='.')
ax.plot(x, T6, label='Tolerância : 1e-6', linestyle='dashed')#, marker='.')  
#ax.plot(xanl, Tanl, label='Solução analítica')
ax.set_xlabel('x [m]', fontsize=14)  
ax.set_ylabel('Temperatura [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=12)
#ax.set_title('Soluções para diferentes tolerâncias considerando 80 volumes finitos') 
ax.legend()  
ax.grid()
fig.show()


## 160 Volumes
fig, ax = plt.subplots(figsize=(10, 6))
#ax.plot(x, T3, label='Tolerância : 1e-3', linestyle='dotted')#, marker='.')
#ax.plot(x, T4, label='Solução Numérica', linestyle='dashed')#, marker='.')
#ax.plot(x, T5, label='Tolerância : 1e-5', linestyle='dashed')#, marker='.')
ax.plot(x, T6, label='Tolerância : 1e-6', linestyle='dashed')#, marker='.')  
ax.plot(xanl, Tanl, label='Solução Analítica')
ax.set_xlabel('x [m]', fontsize=14)  
ax.set_ylabel('Temperatura [K]', fontsize=14)
ax.tick_params(axis='both', which='major', labelsize=12)
#ax.set_title('Soluções para diferentes tolerâncias considerando 80 volumes finitos') 
ax.legend()  
ax.grid()
fig.show()
