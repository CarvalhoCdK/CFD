from copyreg import add_extension
import numpy as np
######

from scipy.sparse import csr_matrix
  
row = np.array([0, 0, 1, 1, 2, 1])
col = np.array([0, 1, 2, 0, 2, 2])
  
# taking data
data = np.array([1, 4, 5, 8, 9, 6])
  
# creating sparse matrix
sparseMatrix = csr_matrix((data, (row, col)), 
                          shape = (3, 3)).toarray()
######
class Assembler(object):
    '''
    '''

    def __init__(self, properties, mesh):
        '''
        '''
        self.k = 0
        self.cp = 0
        self.h = 0

        self.mesh = mesh


    def build_element(self):
        '''
        '''
        k = self.k
        cp = self.cp

        deltax = self.mesh.deltax
        deltay = self.mesh.deltay

        border = self.mesh.border

        
        Aw = k / cp / deltay
        Ae = k / cp / deltay
        As = k / cp / deltax
        An = k / cp / deltax




