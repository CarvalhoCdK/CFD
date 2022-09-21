import numpy as np


class Mesh(object):
    """
    """

    def __init__(self, nx, ny, lx, ly):
        '''
        nx : int
            Number of nodes in 'x' direction
        ny : int
            Number of nodes in 'y' direction
        lx : float
            Length in 'x' direction
        ly : float
            Length in 'y' direction
        '''
        self.nx = nx
        self.ny = ny
        self.lx = lx
        self.ly = ly

        self.deltax = lx / nx
        self.deltay = ly / ny

        self.nodes = self.map()

        self.borders = self.border()


    def map(self):
        '''
        '''

        # Get grid nodes coordinates:
        nx = self.nx
        ny = self.ny
        lx = self.lx
        ly = self.ly

        x = np.linspace(0, lx, nx+1)
        y = np.linspace(0, ly, ny+1)
        grid = np.meshgrid(x, y)
        xv = grid[0].ravel()
        yv = grid[1].ravel()

        nodes = np.column_stack((xv, yv))

        return nodes
    
        
    def border(self):
        '''
             n
             |
         w---p---e
             |
             s
        
        identifica os nós vizinhos por seus números 
        na ordem [w, e, s, n]. O valor -1 indica que
        o nó é adjacente à respectiva fronteira.

        '''
        
        lx = self.lx
        ly = self.ly
        nx = self.nx+1
        ny = self.ny+1
        deltax = self.deltax
        deltay = self.deltay

        margin = 1e-6 * min(deltax, deltay)

        borders = np.empty([nx * ny, 4])
        for i, el in enumerate(self.nodes):
                  
            xn = el[0]
            yn = el[1]

            border = np.array([i - 1, i + 1, i - nx, i + nx])
 
            wb = abs(xn) <= margin
            eb = abs(xn - lx) <= margin
            sb = abs(yn) <= margin
            nb = abs(yn - ly) <= margin

            bc = np.array([wb, eb, sb, nb])

            border[bc] = -1


            borders[i, :] = border

        return borders
        