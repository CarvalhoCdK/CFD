import numpy as np


import numpy as np


def cds(model, type):
        """
        Condições de contorno criadas apenas para temperatura prescrita,
        num problema unidimensional
        """
        Gamma = model['Gamma']
        Rho = model['Rho']
        U = model['U']
        dx = model['dx']
        phi0 = model['Phi_0']
        phi1 = model['Phi_1']
                
        d = Gamma / dx # Auxiliar
        f = Rho*U  # auxiliar

        b = 0

        if type == 'Internal':
            ap = 2*d
            aw = d + f/2
            ae = d - f/2

        elif type == 'Left':
            ap = 3*d + f/2
            aw = 0
            ae = d - f/2
            b = 2*d*phi0 + f*phi0

        else: # 'Right'
            ap = 3*d - f/2
            aw = d + f/2
            ae = 0
            b = 2*d*phi1 - f*phi1
                
        return np.array([ap, aw, ae, b])
    

def uds(model, type):
    """
    """
    Gamma = model['Gamma']
    Rho = model['Rho']
    U = model['U']
    dx = model['dx']
            
    gdx = Gamma / dx  # auxiliar
    ru = Rho*U  # auxiliar
    b = 0
    if type == 'Internal':
        ap = ru + 2*gdx
        aw = ru + gdx
        ae = gdx

    elif type == 'Left':
        ap = ru + gdx
        aw = 0
        ae = gdx

    else: # 'Right'
        ap = ru + gdx
        aw = ru + gdx
        ae = 0
        
    return np.array([ap, aw, ae, b])

'''
def wuds(model):
    """
    """
    Gamma = model['Gamma']
    Rho = model['Rho']
    U = model['U']
    dx = model['dx']

    gdx = Gamma / dx  # auxiliar
    ru = Rho*U  # auxiliar

    Pe = ru*1 / gdx
    Pe2 = Pe**2

    alfa = Pe**2 / (10 + 2*Pe2)
    beta = (1 + 0.005*Pe2) / (1 + 0.05*Pe2)
            
    aw = ru*(0.5 + alfa) + beta*gdx
    ae = -ru*(0.5 - alfa) + beta*gdx
    ap = aw + ae
    
    return np.array([ap, aw, ae])
'''