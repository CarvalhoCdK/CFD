from msilib.schema import Class
import numpy as np


class Mesh(object):
    """
    A class used to create and manipulate retangular meshes.

    y
    |
    6---7---8
    | 2 | 3 |
    3---4---5
    | 0 | 1 |
    0---1---2 --> x
    ...

    Attributes
    ----------
    elements : list
    
    Methods
    -------
    map()
        Build the elements with respective center coordinates,
        for finite volumes method.
    """

    def __init__(self, nx, ny, lx, ly):
        '''
        nx : int
            Number of elements in 'x' direction
        ny : int
            Number of elements in 'y' direction
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

        self.elements = self.map()

        self.borders = self.border()


    def map(self):
        '''
        xxxx
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

        # 
        elements = np.empty([nx * ny, 2])
        for el in range(nx * ny):
            line = el // nx

            n0 = el + line
            n1 = n0 + 1
            n2 = n1 + nx
            n3 = n2 + 1

            xP = (xv[n1] + xv[n0]) / 2
            yP = (yv[n2] + yv[n0]) / 2

            elements[el,:] = [xP, yP]

        return elements

        
    def border(self):
        '''
        '''
        lx = self.lx
        ly = self.ly
        nx = self.nx
        ny = self.ny
        deltax = self.deltax
        deltay = self.deltay

        margin = 1e-6 * lx

        borders = np.empty([nx * ny, 4])
        for i, el in enumerate(self.elements):
                  
            xP = el[0]
            yP = el[1]

            border = np.array([i - 1, i + 1, i - nx, i + nx])
 
            wb = abs(xP - deltax/2) <= margin
            eb = abs(xP + deltax/2 - lx) <= margin
            sb = abs(yP - deltay/2) <= margin
            nb = abs(yP + deltay/2 - ly) <= margin

            bc = np.array([wb, eb, sb, nb])

            border[bc] = -1


            borders[i, :] = border

        return borders


class Mesh1D(object):
    """
    """

    def __init__(self, nx, ny, lx, ly):
        pass

