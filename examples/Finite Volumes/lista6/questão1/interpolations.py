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
            ap = 2*d + f
            aw = d + f
            ae = d

        elif type == 'Left':
            ap = 3*d + f
            aw = 0
            ae = d
            b = 2*d*phi0 + f*phi0

        else: # 'Right'
            ap = 3*d + f
            aw = d + f
            ae = 0
            b = 2*d*phi1
                
        return np.array([ap, aw, ae, b])


def wuds(model, type):
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

        Pe = Rho*U*dx / Gamma
        alfa = Pe**2 / (10 + 2*Pe*2)
        beta = (1 + 0.005*Pe**2) / (1 + 0.05*Pe**2)

        # Debug print
        print(f'alfa: {alfa}')
        print(f'beta: {beta}')

        b = 0

        if type == 'Internal':
            ap = 2*beta*d + 2*alfa*f
            aw = beta*d + (0.5 + alfa)*f
            ae = beta*d - (0.5 - alfa)*f

        elif type == 'Left':
            ap = (beta + 2)*d + (0.5 + alfa)*f
            aw = 0
            ae = beta*d - (0.5 - alfa)*f
            b = 2*d*phi0 + f*phi0

        elif type == 'Right': # 'Right'
            ap = (beta + 2)*d - (0.5 - alfa)*f
            aw = beta*d + (0.5 + alfa)*f
            ae = 0
            b = 2*d*phi1 - f*phi1
                
        return np.array([ap, aw, ae, b])