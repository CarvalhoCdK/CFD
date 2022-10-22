import numpy as np

from numba import jit, njit


def tdma(A, B, T0, nelx: int, nely: int, tol=1e-4)->np.ndarray:
    """
    """
    n = int(nelx*nely-1)
    print(f'nelx:{nelx}')
    print(f'nely:{nely}')
    print(f'n out:{nelx*nely}')

    def _get_coefficients(A, n):

        a = np.zeros(n+1)
        b = np.zeros(n+1)
        c = np.zeros(n+1)
        print(f'n in:{n}')
        # 1
        a[0] = A[0,0]
        b[0] = A[0,1]
        # 5
        a[n] = A[n,n]
        c[n] = A[n,n-1]
        # 2, 3, 4
        for i in np.arange(1, n):
            a[i] = A[i,i]   
            b[i] = A[i,i+1] 
            c[i] = A[i,i-1]

        return a, b, c

    a, b, c = _get_coefficients(A, n)
    print(A)
    print(f'a: {a}')
    print(f'b: {b}')
    print(f'c: {c}')

    def iterate():
        p = np.zeros(n)
        q = np.zeros(n)
        p[0] = -b[0] / a[0]
        q[0] = d[0] / a[0]


    #while error > tol:

    #    solution = iterate()
    #    error = 1

    #return solution