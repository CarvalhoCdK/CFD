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
    max_it = 1e4


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