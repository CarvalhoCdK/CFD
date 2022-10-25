import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from time import perf_counter
from scipy.sparse.linalg import isolve

from interpolations import cds, uds
from solver import tdma
from plots import plot1, compare


def lista6_1(n, tol, model, solver, interpolation):
 
    interpolation_map = {'CDS': cds, 'UDS' : uds}#, 'wuds' : wuds}
    solver_map = {'TDMA' : tdma}
    
    A = np.zeros((n,n))
    B = np.zeros(n)
    # ap, aw, ae, b
    # Primeiro volume
    a = interpolation_map[interpolation](model, 'Left')
    A[0, 1] = -a[2]
    A[0, 0] = a[0]

    # Volumes internos
    for i in range(1,n-1):

        a = interpolation_map[interpolation](model, 'Internal')
        A[i, i-1] = -a[1]
        A[i, i+1] = -a[2]
        A[i, i] = a[0]
        
    # Último volume
    a = interpolation_map[interpolation](model, 'Right')
    A[n-1,n-1] = a[0]
    A[n-1, n-2] = -a[1]
    
    B[n-1] = a[3]
    
    print(f'Iniciando solver: {solver}')
    T = solver_map[solver](A, B)

    return T

# Parâmetros
n = 100
L = 1
tol = 1e-3

model = {'U' : 0.0,
        'Gamma' : 1.0,
        'Rho' : 1.0,
        'dx' : 1.0/n,
        'Phi_0' : 0.0,
        'Phi_1' : 1.0,
        'L' : 1.0}

t1 = 'u = 0'
t2 = 'u = 10'
t3 = 'u = 50'
t4 = 'u = 100'

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


dx = model['dx']

x = np.linspace(dx, 1-dx, n)


## u = 0
TA = pd.DataFrame()
TA[t1] = lista6_1(n, tol, model, 'TDMA', 'CDS')

t_anl = analitic(x, model)

ax, fig = compare(x, TA[t1], t_anl)
fig.show()

## u = 10
model['U'] = 10
TA[t2] = lista6_1(n, tol, model, 'TDMA', 'CDS')
t_anl = analitic(x, model)

ax, fig = compare(x, TA[t2], t_anl)
fig.show()

## u = 50
model['U'] = 50
TA[t3] = lista6_1(n, tol, model, 'TDMA', 'CDS')

## u = 100
model['U'] = 100
TA[t4] = lista6_1(n, tol, model, 'TDMA', 'CDS')


ax, fig = plot1(x, TA, t1,t2,t3,t4)
fig.show()


np.set_printoptions(precision=2)
#print(f'TDMA : {T_A}')

#print(T_tdma[0:10])
