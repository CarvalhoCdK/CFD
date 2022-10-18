import numpy as np
from time import perf_counter

from solver import tdma
from solver import jacobi
from solver import gauss_seidel

from mesh import Mesh


def mesh_test():
    """
    """
    nx = 3
    ny = 3
    lx = 3
    ly = 3

    mesh = Mesh(nx, ny, lx, ly)

    print(mesh.elements)

    print(mesh.borders)


## TEST DATASETS: 1D AND 2D
def test_data_2D():
    """
    Test problem setup
    Incropera 4.3
    """
    A = np.array([[-4, 1, 1, 0, 0, 0, 0, 0],
              [2, -4, 0, 1, 0, 0, 0, 0],
              [1, 0, -4, 1, 1, 0, 0, 0],
              [0, 1, 2, -4, 0, 1, 0, 0],
              [0, 0, 1, 0, -4, 1, 1, 0],
              [0, 0, 0, 1, 2, -4, 0, 1],
              [0, 0, 0, 0, 2, 0, -9, 1],
              [0, 0, 0, 0, 0, 2, 2, -9]])

    B = np.array([-1000,-500,-500,0,-500,0,-2000,-1500])

    direct_solution = np.linalg.solve(A, B)

    return {'A' : A, 'B' : B, 'Solução direta' : direct_solution}


def test_data_1D():
    """
    """
    A = -np.array([[-2, 1, 0, 0],
                  [1, -2, 1, 0],
                  [0, 1, -2, 1],
                  [0, 0, 1, -2]])

    B = -np.array([-100, 0, 0, -1000])

    direct_solution = np.linalg.solve(A, B)

    return {'A' : A, 'B' : B, 'Solução direta' : direct_solution}


## SOLVER TESTS
def test_tdma():
    """
    """
    data = test_data_1D()
    A = data['A']
    B = data['B']
    t0 = data['Solução direta']

    start_time = perf_counter()

    solution = tdma(A, B)

    end_time = perf_counter()

    #print(f'Solução direta : {t0}')
    print(f'TDMA T : {solution}')
    print(f'Tempo de execução: {end_time - start_time}')
    

def test_jacobi():
    """
    """
    data = test_data_1D()
    A = data['A']
    B = data['B']
    t0 = data['Solução direta']
    t_initial = np.ones(A.shape[0])

    
    start_time = perf_counter()

    solution, it = jacobi(A, B, t_initial)

    end_time = perf_counter()

    #print(f'Solução direta : {t0}')
    print(f'Jacobi T : {solution}')
    print(f'# de iterações : {it}')
    print(f'Tempo de execução: {end_time - start_time}')


def test_gauss_seidel():
    """
    """
    data = test_data_1D()
    A = data['A']
    B = data['B']
    t0 = data['Solução direta']

    start_time = perf_counter()

    solution, it = gauss_seidel(A, B)

    end_time = perf_counter()

    #print(f'Solução direta : {t0}')
    print(f'Gauss-Seidel T : {solution}')
    print(f'# de iterações : {it}')
    print(f'Tempo de execução: {end_time - start_time}')


## RUN TESTS
np.set_printoptions(precision=2)
test_jacobi()
test_gauss_seidel()
test_tdma()


# def test(func): 
#     # storing the function in a variable 
    
#     data = test_data_1D()
#     A = data['A']
#     B = data['B']
#     t0 = data['Solução direta']

    
#     start_time = perf_counter()

#     solution, it = func(A, B)

#     end_time = perf_counter()

#     #print(f'Solução direta : {t0}')
#     print(f'Jacobi T : {solution}')
#     print(f'# de iterações : {it}')
#     print(f'Tempo de execução: {end_time - start_time}')

# test(jacobi)
# test(gauss_seidel)
# test(tdma)