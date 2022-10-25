import numpy as np
from time import perf_counter

from solver import tdma



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
    data = test_data_2D()
    A = data['A']
    B = data['B']
    t0 = data['Solução direta']
    print(A.shape)
    nelx = 2
    nely = 4

    start_time = perf_counter()
    print(f'nelx:{nelx}')
    print(f'nely:{nely}')
    print(f'n:{nelx*nely}')

    solution = tdma(A, B, np.ones(nelx*nely), nelx, nely)

    end_time = perf_counter()

    #print(f'Solução direta : {t0}')
    print(f'TDMA T : {solution}')
    print(f'Tempo de execução: {end_time - start_time}')

test_tdma()