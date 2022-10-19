import numpy as np


def tdma(A, B, t_initial=1.0):
    """
    1D version
    [A][t] = [B]
    """

    n = A.shape[0]
    t = np.ones(n) * t_initial

    ## Rearange A, B terms in a, b, c, d
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)
    d = np.zeros(n)

    # Boundary volumes terms
    a[0] = A[0,0]
    b[0] = A[0,1]
    d[0] = B[0]

    a[n-1] = A[n-1,n-1]
    c[n-1] = A[n-1,n-2]
    d[n-1] = B[n-1]

    # Internal volumes terms
    for i in np.arange(1,n-1):#range(1, n-1):
      a[i] = A[i,i]   # Ap
      b[i] = A[i,i+1] # Ae
      c[i] = A[i,i-1] # Aw
      d[i] = B[i]

    ## Foward loop for p and q
    p = np.zeros(n)
    q = np.zeros(n)
    p[0] = -b[0] / a[0]
    q[0] = d[0] / a[0]

    for i in np.arange(1,n):#range(1, n):
      p[i] = -b[i] / (a[i] + c[i]*p[i-1])
      q[i] = (d[i] - c[i]*q[i-1]) / (a[i] + c[i]*p[i-1])

   
    ## Backward loop for t
    t[n-1] = q[n-1]

    for i in np.arange(n-2,-1,-1):#range(n-2,-1,-1):
      t[i] = p[i]*t[i+1] + q[i]

    return t


def gauss_seidel(A, B, tol, t_initial=1.0):
    """
    Linear system solver using Jacobi's method.
    [A][t] = [B]
    """
    n = A.shape[0]
    t = np.ones(n) * t_initial

    diff = t_initial[0]
    #tol = 1e-3#t_initial[0] * 0.01
    max_it = 1e6
    it = 0

    while diff >= tol:

        t0 = np.copy(t)

        for i in range(n):

          #App = np.append(A[i,:i], A[i, i+1:])
          #Tpp = np.append(t[:i], t[i+1:])
            
          #t[i] = (-np.dot(App, Tpp) + B[i]) / A[i, i]

          pfront = np.dot(A[i,:i], t[:i])
          pback = np.dot(A[i, i+1:], t[i+1:])

          t[i] = (-pfront - pback + B[i]) / A[i, i]

        diff = np.max(np.abs(t-t0))

        it += 1
        if it > max_it:
            print('Excedido limite de iterações')
            break

    return t, it


def jacobi(A, B, tol, t_initial):
    """
    Linear system solver using Jacobi's method.
    [A][t] = [B]
    """
    n = A.shape[0]
    t = np.copy(t_initial, subok=True)

    diff = t_initial[0]
    #tol = t_initial[0] * 1e-6
    max_it = 1e6
    it = 0

    while abs(diff) >= tol:

        t0 = np.copy(t)

        for i in range(n):

          #App = np.append(A[i,:i], A[i, i+1:])
          #Tpp = np.append(t0[:i], t0[i+1:])
            
          #t[i] = (-np.dot(App, Tpp) + B[i]) / A[i, i]

          pfront = np.dot(A[i,:i], t0[:i])
          pback = np.dot(A[i, i+1:], t0[i+1:])

          t[i] = (-pfront - pback + B[i]) / A[i, i]

        diff = np.max(np.abs(t-t0))

        it += 1
        #print(f'it : {it}')
        #print(f'Iguais: {t == t0}')
        if it > max_it:
            print('Excedido limite de iterações')
            break

    return t, it