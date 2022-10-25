import matplotlib.pyplot as plt
import numpy as np


def compare(x, ta, tb):

    """
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    ax.plot(x, ta, label='ta')
    ax.plot(x, tb, label='tb', linestyle='dashed')#, marker='+')  

    ax.set_xlabel('x', fontsize=14)  
    ax.set_ylabel('Temperatura', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend()  
    ax.grid()
    #fig.show()

    return ax, fig


def plot1(x, TA, t1,t2,t3,t4):

    """
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.plot(x, TA[t1], label=t1)
    ax.plot(x, TA[t2], label=t2, linestyle='dashed')#, marker='+')
    ax.plot(x, TA[t3], label=t3, linestyle='dotted')#, marker='.')
    ax.plot(x, TA[t4], label=t4, linestyle='dashdot')#, marker='.')  

    ax.set_xlabel('x', fontsize=14)  
    ax.set_ylabel('Temperatura', fontsize=14)
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.legend()  
    ax.grid()
    #fig.show()

    return ax, fig


def plot4():
    """
    """
    fig, axs = plt.subplots(2, 2)
    axs[0, 0].plot(x, y)
    axs[0, 0].set_title('Axis [0, 0]')

    axs[0, 1].plot(x, y, 'tab:orange')
    axs[0, 1].set_title('Axis [0, 1]')

    axs[1, 0].plot(x, -y, 'tab:green')
    axs[1, 0].set_title('Axis [1, 0]')

    axs[1, 1].plot(x, -y, 'tab:red')
    axs[1, 1].set_title('Axis [1, 1]')

    for ax in axs.flat:
        ax.set(xlabel='x-label', ylabel='y-label')

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    for ax in axs.flat:
        ax.label_outer()