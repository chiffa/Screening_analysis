__author__ = 'ank'

from matplotlib import pyplot as plt

debug = 0

def rs(matrix, name):
    plt.title(name)
    plt.imshow(matrix, interpolation='nearest')
    plt.colorbar()
    plt.show()
    plt.clf()

def debug_wrapper(funct):

    def check_matrix(*args,**kwargs):
        result = funct(*args, **kwargs)
        if debug:
            rs(result, funct.__name__)
        return result

    return check_matrix