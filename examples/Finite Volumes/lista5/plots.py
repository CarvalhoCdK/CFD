def plots_volumes(L, D, model):
    ## PLOTS
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
    ax.plot(x10, T10, label='10 Volumes', linestyle='dashed')#, marker='.')
    ax.plot(x20, T20, label='20 Volumes', linestyle='dashed')#, marker='.')
    ax.plot(x40, T40, label='40 Volumes', linestyle='dashed')#, marker='.')
    ax.plot(x80, T80, label='80 Volumes', linestyle='dashed')#, marker='.')  
    ax.plot(x160, T160, label='160 Volumes', linestyle='dashed')#, marker='.')
    ax.plot(xanl, Tanl, label='Solução analítica')
    ax.set_xlabel('x [m]')  
    ax.set_ylabel('Temperatura [K]')  
    #ax.set_title("Simple Plot") 
    ax.legend()  
    ax.grid()
    fig.show()