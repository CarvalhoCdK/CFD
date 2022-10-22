import matplotlib.pyplot as plt
import numpy as np


def plot_field(nx, ny,lx ,ly, values):
    '''
    '''

    x = np.linspace(0, lx, nx+1)
    y = np.linspace(0, ly, ny+1)
    x, y = np.meshgrid(x, y)

    z = values.reshape((nx, ny))
    z_min, z_max = np.abs(z).min(), np.abs(z).max()

    fig, ax = plt.subplots()

    c = ax.pcolormesh(x, y, z, cmap='coolwarm', vmin=z_min, vmax=z_max)
    ax.set_title('pcolormesh')
    # set the limits of the plot to the limits of the data
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    fig.colorbar(c, ax=ax)

    plt.show()