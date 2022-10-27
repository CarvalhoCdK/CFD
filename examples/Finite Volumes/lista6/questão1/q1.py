import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from time import perf_counter
from scipy.sparse.linalg import isolve

from interpolations import cds, uds, wuds
from solver import tdma
from plots import plot1, compare


def lista6_1(n, tol, model, solver, interpolation):
 
    interpolation_map = {'CDS': cds, 'UDS' : uds, 'WUDS' : wuds}
    solver_map = {'TDMA' : tdma}
    
    A = np.zeros((n,n))
    B = np.zeros(n)
    # ap, aw, ae, b
    # Primeiro volume
    a = interpolation_map[interpolation](model, 'Left')
    A[0, 1] = -a[2]
    A[0, 0] = a[0]
    B[0] = a[3]

    # Volumes internos
    for i in range(1,n-1):

        a = interpolation_map[interpolation](model, 'Internal')
        A[i, i-1] = -a[1]
        A[i, i+1] = -a[2]
        A[i, i] = a[0]
        B[i] = a[3]
        
    # Último volume
    a = interpolation_map[interpolation](model, 'Right')
    A[n-1,n-1] = a[0]
    A[n-1, n-2] = -a[1]
    B[n-1] = a[3]
    
    print(f'Iniciando solver: {solver}')
    T = solver_map[solver](A, B)

    # Debug print
    u = model['U']
    print(f'U: {u}')
    #print(f'{interpolation}')
    #print(f'A: {A}')
    #print(f'B: {B}')
    #print(f'T : {T}')

    return T # np.append(np.append(t0, T), tf)


# Solução analítica
def analitic(x, model):

    if model['U'] > 0:
        pe = model['Rho']*model['U']*x/model['Gamma']
        pl = model['Rho']*model['U']*model['L']/model['Gamma']

        r = (np.exp(pe) - 1) / (np.exp(pl) - 1)

        t = (model['Phi_1'] - model['Phi_0'])*r + model['Phi_0']
    
    else:
        r = x/model['L']
        t = (model['Phi_1'] - model['Phi_0'])*r + model['Phi_0']

    return t


# Parâmetros
n = 60
L = 1
tol = 1e-6

model = {'U' : 0.0,
        'Gamma' : 1.0,
        'Rho' : 1.0,
        'dx' : 1.0/n,
        'Phi_0' : 0.0,
        'Phi_1' : 1.0,
        'L' : 1.0}


dx = model['dx']

TA = pd.DataFrame()
TB = pd.DataFrame()
TC = pd.DataFrame()

x = np.linspace(dx/2, 1-dx/2, n)

x_anl = np.linspace(dx/2, 1-dx/2, 10*n)
x_anl = np.append(np.append(0, x_anl), L)
t_anl = analitic(x_anl, model)

## u = 0
model['U'] = 0
velocidade = 'u = 0'

TA[velocidade] = lista6_1(n, tol, model, 'TDMA', 'CDS')
TB[velocidade] = lista6_1(n, tol, model, 'TDMA', 'UDS')
TC[velocidade] = lista6_1(n, tol, model, 'TDMA', 'WUDS')
t_anl = analitic(x_anl, model)

ax, fig = compare(x, TA[velocidade], TB[velocidade], TC[velocidade], x_anl, t_anl)
ax.set_title(velocidade) 
#fig.show()


## u = 10
model['U'] = 10
velocidade = 'u = 10'

TA[velocidade] = lista6_1(n, tol, model, 'TDMA', 'CDS')
TB[velocidade] = lista6_1(n, tol, model, 'TDMA', 'UDS')
TC[velocidade] = lista6_1(n, tol, model, 'TDMA', 'WUDS')
t_anl = analitic(x_anl, model)

ax, fig = compare(x, TA[velocidade], TB[velocidade], TC[velocidade], x_anl, t_anl)
ax.set_title(velocidade) 
fig.show()

## u = 50
model['U'] = 50
velocidade = 'u = 50'

TA[velocidade] = lista6_1(n, tol, model, 'TDMA', 'CDS')
TB[velocidade] = lista6_1(n, tol, model, 'TDMA', 'UDS')
TC[velocidade] = lista6_1(n, tol, model, 'TDMA', 'WUDS')
t_anl = analitic(x_anl, model)

ax, fig = compare(x, TA[velocidade], TB[velocidade], TC[velocidade], x_anl, t_anl)
ax.set_title(velocidade) 
fig.show()


## u = 100
model['U'] = 100
velocidade = 'u = 100'

TA[velocidade] = lista6_1(n, tol, model, 'TDMA', 'CDS')
TB[velocidade] = lista6_1(n, tol, model, 'TDMA', 'UDS')
TC[velocidade] = lista6_1(n, tol, model, 'TDMA', 'WUDS')
t_anl = analitic(x_anl, model)

ax, fig = compare(x, TA[velocidade], TB[velocidade], TC[velocidade], x_anl, t_anl)
ax.set_title(velocidade) 
fig.show()

