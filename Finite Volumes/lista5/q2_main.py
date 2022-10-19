import numpy as np
import matplotlib.pyplot as plt
from time import perf_counter
from scipy.sparse.linalg import isolve

from mesh import Mesh
from solver import jacobi
from solver import tdma
from solver import sor
from solver import gauss_seidel


def lista5_2(n, tol, model, solver):
    #n = 100

    L = model['L']

    def dirichlet(model):
        """
        Aplica condição de contorno de temperatura prescrita.
        Equação do volume tem forma:
        ap Tp = ae Te + aw Tw + B + a0 T0
        """
        
        q = model['q']
        Tb = model['Tb']
        K = model['K']

        dx = model['dx']
                
        ap = 3 / dx
        aw = 1 / dx
        ae = 0
        b = 2 * Tb / dx + q * dx / K
        
        return np.array([ap, aw, ae, b])
       
    def newmann(model):
        """
        Aplica a condição de contorno de fluxo prescrito.
        Equação do volume tem forma:
        ap Tp = ae Te + aw Tw + B + a0 T0
        """
        
        q = model['q']
        q1 = model['q1']
        K = model['K']

        dx = model['dx']
              
        ap = 1 / dx
        aw = 1 / dx
        ae = 1 / dx
        b = 1 / K * (q1 + q * dx)
        
        return np.array([ap, aw, ae, b])


    def internal_volume(model):
        """
        Implementa um volume interno genérico.
        Equação do volume tem forma:
        ap Tp = ae Te + aw Tw + B + a0 T0
        """
        
        q = model['q']
        K = model['K']
    
        dx = model['dx']
    
        ap = 2 / dx
        aw = ae = 1 / dx
        b = q * dx / K
        
        return np.array([ap, aw, ae, b])

    # Constroi matrizes da forma:
    # [A][T] = [b]+[A0][T0]
    # [ap, aw, ae, b]

    A = np.zeros((n,n))
    B = np.zeros(n)

    for i in range(n):

        if i == 0:
            # Primeiro volume
            a = newmann(model)
            A[i, i+1] = -a[2]

            #print(f'{i} : Dirichlet')

        elif i == n-1:
            # Ultimo volume
            a = dirichlet(model)
            A[i, i-1] = -a[1]

            #print(f'{i} : Newmann')

        else:
            a = internal_volume(model)
            A[i, i-1] = -a[1]
            A[i, i+1] = -a[2]
            #print(f'{i} : Interno')
        
        A[i, i] = a[0]
        
        B[i] = a[3]

    #print(A)

    # SOLVER

    # Condição inicial
    T = model['T0'] * np.ones(n)


    if solver == 'TDMA':
        print('Iniciando solver: TDMA')
        T00 = np.copy(T)

        start_time = perf_counter()

        T = tdma(A, B, T00)
    
        end_time = perf_counter()
        print(f'Solução {n} volumes')
        print(f'Tempo de execução: {end_time - start_time} \n')
    
    
    elif solver == 'Jacobi':
        print('Iniciando solver: Jacobi')
        T00 = np.copy(T)

        start_time = perf_counter()
        
        T, itj = jacobi(A, B, tol, T00)

        end_time = perf_counter()
        print(f'Solução {n} volumes')
        print(f'Número de iterações: {itj}')
        print(f'Tempo de execução: {end_time - start_time} \n')
    

    elif solver == 'Gauss-Seidel':
        print('Iniciando solver: Gauss-Seidel')
        T00 = np.copy(T)

        start_time = perf_counter()
        
        T, itj = gauss_seidel(A, B, tol, T00)

        end_time = perf_counter()
        print(f'Solução {n} volumes')
        print(f'Número de iterações: {itj}')
        print(f'Tempo de execução: {end_time - start_time} \n')


    elif solver == 'S.O.R':
        print('Iniciando solver: S.O.R')
        T00 = np.copy(T)
        w = 1.5

        start_time = perf_counter()
        
        T, itj = sor(A, B, tol, w, T00)

        end_time = perf_counter()
        print(f'Solução {n} volumes')
        print(f'Número de iterações: {itj}')
        print(f'Tempo de execução: {end_time - start_time} \n')


    elif solver == 'Gradientes Conjugados':

        print('Iniciando solver: Gradientes Conjugados')
        T00 = np.copy(T)

        start_time = perf_counter()
        
        T, itj = isolve.cg (A, B, tol=1e-3)

        end_time = perf_counter()
        print(f'Solução {n} volumes')
        print(f'Número de iterações: {itj}')
        print(f'Tempo de execução: {end_time - start_time} \n')

    return T


n = 100
L = 3
tol = 1e-3

model = {'K' : 1.0,
        'Tb' : 283.15,
        'q' : 7.0,
        'q1' : 10.0,
        'dx' : L/(n),
        'L' : L,
        'T0' : 100}

T_tdma = lista5_2(n, tol, model, 'TDMA')
T_jacobi = lista5_2(n, tol, model, 'Jacobi')
T_gauss = lista5_2(n, tol, model, 'Gauss-Seidel')
T_sor = lista5_2(n, tol, model, 'S.O.R')
T_cg = lista5_2(n, tol, model, 'Gradientes Conjugados')

np.set_printoptions(precision=2)
#print(f'TDMA : {T_tdma}')
#print(f'Jacobi : {T_jacobi}')
#print(f'Gauss-Seidel : {T_gauss}')
#print(f'S.O.R : {T_sor}')
print(T_sor-T_tdma)
