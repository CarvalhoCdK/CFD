import numpy as np


class WUDS(object):
    """
    """
    def __init__(self, model) -> None:
        self.model = model        
        

    def internal(self):
        """
        """
        u = self.model['u']
        v = self.model['v']
        rho = self.model['rho']
        gamma = self.model['gamma']
        deltax = self.model['deltax']
        deltay = self.model['deltay']

        fx = rho*u*deltay
        fy = rho*v*deltax
        dx = gamma*deltay / deltax
        dy = gamma*deltax / deltay
        # x :>
        Pe = rho*u*dx / gamma
        alfax = Pe**2 / (10 + 2*Pe*2)
        betax = (1 + 0.005*Pe**2) / (1 + 0.05*Pe**2)
        # y :>
        Pe = rho*v*dx / gamma
        alfay = Pe**2 / (10 + 2*Pe*2)
        betay = (1 + 0.005*Pe**2) / (1 + 0.05*Pe**2)

        Aw =  (0.5 + alfax)*fx + betax*dx
        Ae = -(0.5 - alfax)*fx + betax*dx
        As =  (0.5 + alfay)*fy + betay*dy
        An = -(0.5 + alfay)*fy + betay*dy
        Ap = Aw + Ae + As + An

        return np.array([Ap, Aw, Ae, As, An, 0])


    def boundary(self, face, t):
        """
        """
        u = self.model['u']
        v = self.model['v']
        rho = self.model['rho']
        gamma = self.model['gamma']
        deltax = self.model['deltax']
        deltay = self.model['deltay']

        dx = gamma*deltay / deltax
        dy = gamma*deltax / deltay

        # x :>
        Pe = rho*u*dx / gamma
        alfax = Pe**2 / (10 + 2*Pe*2)
        # y :>
        Pe = rho*v*dx / gamma
        alfay = Pe**2 / (10 + 2*Pe*2)

        A =  np.zeros(6)
        A[-1] = t  # B

        # [Ap, Aw, Ae, As, An, B]
        # [0,  1,  2,  3,  4,  5]
        if face =='W':
            A[0] = 0.5 + alfax
            A[2] = -0.5 + alfax  # Ae

        if face =='E':
            A[0] = 0.5 + alfax  # Ap
            A[1] = -0.5 + alfax  # Aw

        if face =='S':
            A[0] = 0.5 + alfax  # Ap
            A[4] = -0.5 + alfax  # An

        if face =='N':
            A[0] = 0.5 + alfax  # Ap
            A[3] = -0.5 + alfax  # As


        return A

